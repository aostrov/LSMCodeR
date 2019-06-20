stimParListRDS <- file.path(LSMCodeRConfig$srcdir,"objects","stimParListForParallel.RDS")
stimulusParametersList <- readRDS(file=stimParListRDS)
# source(file.path(LSMCodeRConfig$scriptsUtility,"createStimParameterList.R"))

myFiles <- dir(imageDir,patt="^[A-Z]{4}-[[:graph:]]",full=T,rec=TRUE)
outPutFiles <- dir(file.path(LSMCodeRConfig$srcdir,"toDiscardEventually"),pattern = "^[A-Z]{4}",full=T)

for (myFile in myFiles[length(myFiles)]) {
  rdsFile <- file.path(LSMCodeRConfig$srcdir,"toDiscardEventually",basename(myFile))
  # lockFile <- sub(".mat",".lock",rdsFile)
  matFileCode <- substring(basename(myFile),1,4)
  if (FALSE) { # file.exists(rdsFile) | file.exists(lockFile)
    print(paste(myFile,"exists or is being worked on. Skipping"))
  } else {
    print(paste("Processing:",myFile))
    # cat("some text here",file=lockFile)
    
    print("Making ROI list")
    matFileROIListByZ <- list()
    for (tectumROIforZ in 1:nrow(tempDF <- tectumROIs[as.character(tectumROIs$matfile)==matFileCode,])) {
      print(tectumROIforZ)
      if (any(is.na(tempDF[1,]))){
        print("Check to see if you've updated and reloaded tectumROIs.csv")
        unlink(lockFile)
      }
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
      if (file.exists(file.path(LSMCodeRConfig$srcdir,"toDiscardEventually",paste(stimulation,".RDS",sep="")))){
        print(paste("Skipping",stimulation))
      } else {
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
        background.start  <- protocolList$sabineProtocolSimple$background.start/z_planes
        background.end    <- protocolList$sabineProtocolSimple$background.end/z_planes
        signal.start      <- protocolList$sabineProtocolSimple$signal.start/z_planes
        signal.end        <- protocolList$sabineProtocolSimple$signal.end/z_planes
        # save the relevant data from each ROI.
        # full data still exists as .mat files
        # but it might be nice to have the background
        # and signal areas saved externally
        rawData <- lapply(rawDataForMatFileByROIs, function(x) {
          lapply(x,
                 function(y) {
                   background <- apply(y[,,c((background.start):(background.end))],3,mean)
                   attr(background,"start") <- background.start
                   attr(background,"end") <- background.end
                   
                   signal <- apply(y[,,c((signal.start):(signal.end))],3,mean)
                   attr(signal,"start") <- signal.start
                   attr(signal,"end") <- signal.end
                   
                   sgolaySignal <- sgolayfilt(signal,p=2)
                   attr(sgolaySignal,"start") <- signal.start
                   attr(sgolaySignal,"end") <- signal.end
                   
                   return(list(background=background,
                               signal=signal,
                               sgolaySignal=sgolaySignal))
                 })
        })
        
        
        saveRDS(rawData,file=file.path(LSMCodeRConfig$srcdir,"toDiscardEventually",paste(stimulation,".RDS",sep="")))
      }
      
      # statisticsList[[stimulation]]<- lapply(
      #   rawDataForMatFileByROIs,getUsefulStatisticsByROI,matFileROIListByZ,
      #   analysisWindow=c((signal.start):(signal.end)),
      #   backgroundWindow=c((background.start):(background.end))
      #   
      # )
      # 
      # attr(statisticsList[[stimulation]],"ROI_Location") <- c(frame.start=stimulusParametersList[[matFileCode]][[stimulation]]$start,
      #                                                         frame.end = stimulusParametersList[[matFileCode]][[stimulation]]$end)
      # attr(statisticsList[[stimulation]],"StimulusBlock") <- c(stimulus=stimulusParametersList[[matFileCode]][[stimulation]]$stimulus,
      #                                                          block=stimulusParametersList[[matFileCode]][[stimulation]]$block)
      # attr(statisticsList[[stimulation]],"matFile") <- myFile
      # 
    }
    # matList <- list()
    # 
    # matList[[substr(stimulation,1,4)]] <- statisticsList
    # saveRDS(matList,file=rdsFile,compress = TRUE)
    # unlink(lockFile)
    # 
    # clean up
    file.h5$close()
    imageDataSlice$close()
    
  } 
}



ABCA <- dir(imageDir,patt="ABCA",rec=T,full=T)
file.h5 <- H5File$new(ABCA, mode = "r")
imageDataSlice<-file.h5[["imagedata"]]

processSingleStimulus.lapply(myList = stimulusParametersList$ABCA[1:4],
                             outputType = "raw",writeNRRD = T,downSampleImage = T,resizeFactor = 2,
                             image = imageDataSlice,offsetValue = 399,
                             outDir = "F:/",z=16)

lapply(stimulusParametersList$ABCA[1:4],processSingleStimulus.lapply,outputType = "raw",
       writeNRRD = T,downSampleImage = T,resizeFactor = 2,
       image = imageDataSlice,offsetValue = 399,
       outDir = "F:/",z=16)

