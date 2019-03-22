tectumROIs <- read.csv(file.path(LSMCodeRConfig$srcdir,"stuff","tectumROI.csv"))
anatomyFiles <- dir(imageDir,patt="^[A-Z]{3}-")
physiologyFilesSP <- dir(imageDir,patt="[A-Z]{4}-[[:graph:]]*SP",full=T,rec=TRUE)
multiplaneFiles <- setdiff(
  dir(imageDir,patt="[A-Z]{4}-[[:graph:]]",full=T,rec=TRUE),
  physiologyFilesSP
  )
# matFile <- "/Volumes/home/hdf5/20181204-gcamp7F-7d-SabineBars-1planeSP.mat"
# matFile <- file.path(imageDir,"AAFA-gen-A-laser3-SabineSimple","AAFA-gen-A-laser3-SabineSimple.mat")
# file.h5 <- H5File$new(file.path(matFile), mode = "r")
# imageDataSlice<-file.h5[["imagedata"]]
pixelSize=0.8
roiEdgeLength <- 26
matListRDS <- file.path(LSMCodeRConfig$srcdir,"objects","matList1.RDS")
stimParListRDS <- file.path(LSMCodeRConfig$srcdir,"objects","stimParList.RDS")

cleanLSMAnalysis <- function(){
  rm(matList)
  rm(stimulusParametersList)
}

if (file.exists(matListRDS) & exists("matList")) {
  print("matList exists and is already loaded. ")
  print("Using the current matList.")
  print("")
} else if (file.exists(matListRDS) & !exists("matList")){
  print("Loading matList from disk.")
  matList <- readRDS(file=matListRDS)
} else {
  print("There is no matList on disk or in memory.")
  print("Starting matList from scratch.")
  print("Go get a coffee, this will take a while.")
  matList <- list()
}

if (file.exists(stimParListRDS) & exists("stimulusParametersList")) {
  print("stimulusParametersList exists and is already loaded. ")
  print("Using the current stimulusParametersList")
  print("")
} else if (file.exists(stimParListRDS) & !exists("stimulusParametersList")){
  print("Loading stimulusParametersList from disk.")
  stimulusParametersList <- readRDS(file=stimParListRDS)
} else {
  print("There is no stimulusParametersList on disk or in memory.")
  print("Starting stimulusParametersList from scratch.")
  stimulusParametersList <- list()
}



# set up file handling
for (myFile in dir(imageDir,patt="[A-Z]{4}-[[:graph:]]",full=T,rec=TRUE)) {
  # foreach(myFile=physiologyFilesSP, .packages = "hdf5r") %dopar% {
  matFileCode <- substring(basename(myFile),1,4)
  animal <- list()
  print(paste("Starting to parse",matFileCode))
  if (is.null(matList[[matFileCode]])){
    # get the transition frames for the mat file
    source(file.path(LSMCodeRConfig$srcdir,"parseImageForGreenFlashAndLSMandStimLogs.R"))
    print("Passed parseImageForGreenFlash...")
    currentStimulusParameters <- makeTrial(myFile)
    # set up ROIs
    
    print("setting up ROIs")
    matFileROIListByZ <- list()
    for (tectumROIforZ in 1:nrow(tempDF <- tectumROIs[as.character(tectumROIs$matfile)==matFileCode,])) {
      print(tectumROIforZ)
      matFileROIListByZ[[paste("z",tempDF[tectumROIforZ,"z"],sep="_")]] = getROIs(roiEdgeLength = roiEdgeLength,
              x=tempDF[tectumROIforZ,"x"],
              y=tempDF[tectumROIforZ,"y"],
              w=tempDF[tectumROIforZ,"w"],
              h=tempDF[tectumROIforZ,"h"])
    }
    print("finished setting up ROIs")
    animal[[basename(myFile)]] <- currentStimulusParameters
    statisticsList <- list()
    for (stimulation in names(currentStimulusParameters)){ # for a matfile get the raw data as Z plane and ROI
      # foreach(stimulation=names(currentStimulusParameters)) %dopar% {
      print(currentStimulusParameters[[stimulation]]$start)
      rawDataForMatFileByROIs <- mapply( # outer lapply breaks down roiList by z-plane
        function(roiList,z) 
        {
        lapply( # inner lapply gets the raw statistics for each ROI
          roiList,
          getROIsRawDataFromHDF5.lapply,
          hdf5Image.mat = imageDataSlice,
          frame.start = currentStimulusParameters[[stimulation]]$start, # 750,
          frame.end = currentStimulusParameters[[stimulation]]$end,
          offset = 399
          )
        },
        matFileROIListByZ,sub("z_","",names(matFileROIListByZ))
      )
      
      saveRDS(rawDataForMatFileByROIs,
              file=file.path(LSMCodeRConfig$srcdir,"objects",paste(stimulation,".RDS",sep="")),
              compress = T)
      # statisticsList[[stimulation]]<- getUsefulStatisticsByROI(rawDataForMatFileByROIs,matFileROIListByZ,
      #                                                          analysisWindow=c(150:750),backgroundWindow=c(0:100))
      statisticsList[[stimulation]]<- lapply(
                                        rawDataForMatFileByROIs,getUsefulStatisticsByROI,matFileROIListByZ,
                                        analysisWindow=c((900/length(rawDataForMatFileByROIs)):(1500/length(rawDataForMatFileByROIs))),
                                        backgroundWindow=c((700/length(rawDataForMatFileByROIs)):(900/length(rawDataForMatFileByROIs)))
                                      )
      
      attr(statisticsList[[stimulation]],"ROI_Location") <- c(frame.start=currentStimulusParameters[[stimulation]]$start,
                                                              frame.end = currentStimulusParameters[[stimulation]]$end)
      attr(statisticsList[[stimulation]],"StimulusBlock") <- c(stimulus=currentStimulusParameters[[stimulation]]$stimulus,
                                                               block=currentStimulusParameters[[stimulation]]$block)
      attr(statisticsList[[stimulation]],"matFile") <- myFile
      
    }
    matList[[substr(stimulation,1,4)]] <- statisticsList
    saveRDS(matList,file=matListRDS,compress = TRUE)
    # clean up
    file.h5$close()
    imageDataSlice$close()
  } else {
    print(paste(matFileCode, "already exists in matList, skipping"))
  }
  stimulusParametersList[[matFileCode]] <- animal
  saveRDS(stimulusParametersList,file=stimParListRDS,compress = TRUE)
  count
}

saveRDS(matList,file=matListRDS,compress = TRUE)

analysisDF <- c()
for (k in 1:length(matList)){
  analysisDF.animals <- c()
  animal <- names(matList[k])
  print(paste("animal:",animal))
  for (j in 1: length(matList[[k]])){
    analysisDF.stimulus <- c()
    stimulus <- substr(names(matList[[k]][j]),nchar(names(matList[[k]][j]))-2,nchar(names(matList[[k]][j])))
    print(paste("stimulus:",stimulus))
    for (i in 1:length(matList[[k]][[j]])) {
      analysisDF.zplane <- c()
      z_plane <- names(matList[[k]][[j]][i])
      analysisDF.zplane <- matList[[k]][[j]][[i]][,c("snr.mean","background.mean","dff.mean")]
      analysisDF.zplane$z_plane <- as.factor(z_plane)
      analysisDF.stimulus <- rbind(analysisDF.stimulus,analysisDF.zplane)
    }
    analysisDF.stimulus$stimulus <- as.factor(stimulus)
    analysisDF.animals <- rbind(analysisDF.animals,analysisDF.stimulus)
  }
  analysisDF.animals$animal <- as.factor(animal)
  analysisDF <- rbind(analysisDF,analysisDF.animals)
}

ggplot(analysisDF,aes(background.mean,dff.mean)) + 
   geom_jitter(aes(color=animal)) # + facet_wrap(~stimulus)


ggplot(subset(analysisDF,animal=="AAMA"),aes(background.mean,snr.mean)) + 
   geom_jitter(aes(color=z_plane))# + facet_wrap(~stimulus)

analysisDF.subset <- subset(analysisDF,background.mean>25 & dff.mean>0.05)

# Go back to matList to find all the statistics for an ROI
matList[[1]][[1]][[1]]["ROI_1",]