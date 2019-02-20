tectumROIs <- read.csv(file.path(LSMCodeRConfig$srcdir,"stuff","tectumROI.csv"))

anatomyFiles <- dir("F:/Imaging/GCaMP7_tests/20181204-g7/",patt="^[A-Z]{3}-")
physiologyFilesSP <- dir("F:/Imaging/GCaMP7_tests/20181204-g7/",patt="[A-Z]{4}-[[:graph:]]*SP",full=T,rec=TRUE)
matFile <- "/Volumes/home/hdf5/20181204-gcamp7F-7d-SabineBars-1planeSP.mat"
file.h5 <- H5File$new(file.path(matFile), mode = "r")
imageDataSlice<-file.h5[["imagedata"]]
roiEdgeLength <- 26
# registerDoParallel(cores = 2)

# lapply(physiologyFilesSP)
matList <- list()
# set up file handling
for (myFile in physiologyFilesSP[2:6]) {
# foreach(myFile=physiologyFilesSP, .packages = "hdf5r") %dopar% {
  # get the transition frames for the mat file
  source(file.path(LSMCodeRConfig$srcdir,"parseImageForGreenFlashAndLSMandStimLogs.R"))
  currentStimulusParameters <- makeTrial(myFile)
  # set up ROIs
  matFileCode <- substring(basename(myFile),1,4)
  roiList <- getROIs(roiEdgeLength = roiEdgeLength,x=tectumROIs[tectumROIs$matfile==matFileCode,"x"],
                     y=tectumROIs[tectumROIs$matfile==matFileCode,"y"],
                     w=tectumROIs[tectumROIs$matfile==matFileCode,"w"],
                     h=tectumROIs[tectumROIs$matfile==matFileCode,"h"]) # ~10um ROIs


  statisticsList <- list()
  for (stimulation in names(currentStimulusParameters)){
    # foreach(stimulation=names(currentStimulusParameters)) %dopar% {
    print(currentStimulusParameters[[stimulation]]$start)
    rawDataForMatFileByROIs <- lapply(
      roiList,
      getROIsRawDataFromHDF5.lapply,
      hdf5Image.mat = imageDataSlice,
      frame.start = currentStimulusParameters[[stimulation]]$start+750,
      frame.end = currentStimulusParameters[[stimulation]]$end,
      offset = 399
    )
    statisticsList[[stimulation]]<- getUsefulStatisticsByROI(rawDataForMatFileByROIs,roiList,
                                                             analysisWindow=c(150:750),backgroundWindow=c(0:100))
    attr(statisticsList[[stimulation]],"ROI_Location") <- c(frame.start=currentStimulusParameters[[stimulation]]$start+750,
                                                            frame.end = currentStimulusParameters[[stimulation]]$end)
    attr(statisticsList[[stimulation]],"StimulusBlock") <- c(stimulus=currentStimulusParameters[[stimulation]]$stimulus,
                                                             block=currentStimulusParameters[[stimulation]]$block)
    
  }
  matList[[substr(stimulation,1,4)]] = statisticsList
  # clean up
  file.h5$close()
  imageDataSlice$close()
}

ggplot(statisticsList[[stimulation]],aes(xpos,ypos,fill=snr.mean)) + 
  geom_raster(interpolate = F) + 
  scale_fill_gradientn(colors = jet(20)) +
  coord_fixed() + scale_y_reverse()



# here I need to combine a few things

# for each ROI I would be able to tell when/in which frame the max occurs
sapply(roiSNR,function(x) {which.max(apply(x[,,90:150], 3, mean))+90})

ggplot(matList[[3]][[1]],aes(xpos,ypos,fill=snr.mean)) + 
  geom_raster(interpolate = F) + 
  scale_fill_gradientn(colors = jet(20)) +
  coord_fixed() + scale_y_reverse()
