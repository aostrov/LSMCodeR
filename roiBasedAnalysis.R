stimParListRDS <- file.path(LSMCodeRConfig$srcdir,"objects","stimParListForParallel.RDS")

cleanLSMAnalysis <- function(){
  rm(matList)
  rm(stimulusParametersList)
}

if (file.exists(stimParListRDS) & exists("stimulusParametersList")) {
  print("stimulusParametersList exists and is already loaded. ")
  onDiskParameters <- readRDS(stimParListRDS)
  newParams <- setdiff(names(onDiskParameters),names(stimulusParametersList))
  if (length(newParams) > 0) {
    stimulusParametersList = onDiskParameters
    print("Using the on disk stimulusParametersList")
  } else {
    print("Using the in memory stimulusParametersList")
  }
} else if (file.exists(stimParListRDS) & !exists("stimulusParametersList")){
  print("There is no parameter list in memory. Loading stimulusParametersList from disk.")
  stimulusParametersList <- readRDS(file=stimParListRDS)
} else {
  print("There is no stimulusParametersList on disk or in memory.")
  print("Starting stimulusParametersList from scratch.")
  source(file.path(LSMCodeRConfig$srcdir,"createStimParamterList.R"))
}



######################################
##### more parallel, less clever #####
######################################


# for (myFile in sample(dir(imageDir,patt="[A-Z]{4}-[[:graph:]]",full=T,rec=TRUE)[1:2],
#                       length(dir(imageDir,patt="[A-Z]{4}-[[:graph:]]",full=T,rec=TRUE)[1:2]))) {

myFiles <- dir(imageDir,patt="[A-Z]{4}-[[:graph:]]",full=T,rec=TRUE)


for (myFile in myFiles) {
  rdsFile <- file.path(LSMCodeRConfig$srcdir,"objects",basename(myFile))
  lockFile <- sub(".mat",".lock",rdsFile)
  matFileCode <- substring(basename(myFile),1,4)
  if (file.exists(rdsFile) | file.exists(lockFile)) {
    print(paste(myFile,"exists or is being worked on. Skipping"))
  } else {
    print(paste("Processing:",myFile))
    cat("some text here",file=lockFile)
    
    print("Making ROI list")
    matFileROIListByZ <- list()
    for (tectumROIforZ in 1:nrow(tempDF <- tectumROIs[as.character(tectumROIs$matfile)==matFileCode,])) {
      print(tectumROIforZ)
      matFileROIListByZ[[paste("z",tempDF[tectumROIforZ,"z"],sep="_")]] = getROIs(roiEdgeLength = roiEdgeLength,
                                                                                  x=tempDF[tectumROIforZ,"x"],
                                                                                  y=tempDF[tectumROIforZ,"y"],
                                                                                  w=tempDF[tectumROIforZ,"w"],
                                                                                  h=tempDF[tectumROIforZ,"h"])
    } # attr(roiList, "ROI_info" ) = attrVector
    
    statisticsList <- list()
    
    print("Extracting raw data and useful statistics")
    
    # open the connection to the .mat file
    file.h5 <- H5File$new(myFile, mode = "r")
    # save a new variable that contains the image data
    # this step can be avoided, and might improve performance on very large datasets
    imageDataSlice<-file.h5[["imagedata"]]
    file.h5$close()

    # data is in the form:
    # [stacks,(channels),slices,rows, columns]
    imageDataSlice.dims <- imageDataSlice$dims
    names(imageDataSlice.dims) <- c('t','c','z','y','x')
    print(imageDataSlice.dims)
    
    for (stimulation in names(stimulusParametersList[[matFileCode]])){ # for a matfile get the raw data as Z plane and ROI
      # foreach(stimulation=names(currentStimulusParameters)) %dopar% {
      cat(".",sep="\n")
      print(stimulusParametersList[[matFileCode]][[stimulation]]$start)
      rawDataForMatFileByROIs <- mapply( # outer lapply breaks down roiList by z-plane
        function(roiList,z) 
        {
          lapply( # inner lapply gets the raw statistics for each ROI
            roiList,
            getROIsRawDataFromHDF5.lapply,
            z,
            hdf5Image.mat = imageDataSlice,
            frame.start = stimulusParametersList[[matFileCode]][[stimulation]]$start, # 750,
            frame.end = stimulusParametersList[[matFileCode]][[stimulation]]$end,
            offset = offset
          )
        },
        roiList=matFileROIListByZ,z=as.integer(sub("z_","",names(matFileROIListByZ))),SIMPLIFY = F
      )
      
      z_planes=attr(stimulusParametersList[[matFileCode]],"imageDimensions")[['z']]
      # save the relevant data from each ROI.
      # full data still exists as .mat files
      # but it might be nice to have the background
      # and signal areas saved externally
      test <- lapply(rawDataForMatFileByROIs, function(x) {
        lapply(x,
               function(y) {
                 background <- apply(y[,,
                                       c((protocolList$sabineProtocolSimple$background.start/z_planes):
                                           (protocolList$sabineProtocolSimple$background.end/z_planes))],3,mean)
                 attr(background,"start") <- protocolList$sabineProtocolSimple$background.start/z_planes
                 attr(background,"end") <- protocolList$sabineProtocolSimple$background.end/z_planes
                 
                 signal <- apply(y[,,c((protocolList$sabineProtocolSimple$signal.start/z_planes):
                                         (protocolList$sabineProtocolSimple$signal.end/z_planes))],3,mean)
                 attr(signal,"start") <- protocolList$sabineProtocolSimple$signal.start/z_planes
                 attr(signal,"end") <- protocolList$sabineProtocolSimple$signal.end/z_planes
                 return(list(background=background,signal=signal))
               })
      })
      saveRDS(test,file=file.path(LSMCodeRConfig$srcdir,"toDiscardEventually",paste(stimulation,".RDS",sep="")))
      
      
      statisticsList[[stimulation]]<- lapply(
        rawDataForMatFileByROIs,getUsefulStatisticsByROI,matFileROIListByZ,
        analysisWindow=c((protocolList$sabineProtocolSimple$signal.start/z_planes):
                           (protocolList$sabineProtocolSimple$signal.end/z_planes)),
        backgroundWindow=c((protocolList$sabineProtocolSimple$background.start/z_planes):
                             (protocolList$sabineProtocolSimple$background.start/z_planes))
      )
      
      attr(statisticsList[[stimulation]],"ROI_Location") <- c(frame.start=stimulusParametersList[[matFileCode]][[stimulation]]$start,
                                                              frame.end = stimulusParametersList[[matFileCode]][[stimulation]]$end)
      attr(statisticsList[[stimulation]],"StimulusBlock") <- c(stimulus=stimulusParametersList[[matFileCode]][[stimulation]]$stimulus,
                                                               block=stimulusParametersList[[matFileCode]][[stimulation]]$block)
      attr(statisticsList[[stimulation]],"matFile") <- myFile
      
    }
    matList <- list()

    matList[[substr(stimulation,1,4)]] <- statisticsList
    saveRDS(matList,file=rdsFile,compress = TRUE)
    unlink(lockFile)

    # clean up
    file.h5$close()
    imageDataSlice$close()
    
  } 
}
