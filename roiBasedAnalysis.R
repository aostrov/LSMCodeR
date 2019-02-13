

anatomyFiles <- dir("F:/Imaging/GCaMP7_tests/20181204-g7/",patt="^[A-Z]{3}-")
physiologyFilesSP <- dir("F:/Imaging/GCaMP7_tests/20181204-g7/",patt="[A-Z]{4}-[[:graph:]]*SP")
matFile <- "/Volumes/home/hdf5/20181204-gcamp7F-7d-SabineBars-1planeSP.mat"
file.h5 <- H5File$new(file.path(matFile), mode = "r")
imageDataSlice<-file.h5[["imagedata"]]


roiList <- getROIs(roiEdgeLength = 26,x=64,y=192,w=600,h=172) # ~10um ROIs
xxx <- lapply(
  roiList, 
  getROIsRawDataFromHDF5.lapply,
  hdf5Image.mat = imageDataSlice, 
  frame.start = 14, 
  frame.end = 1814, 
  offset = 399
  )

stats.xxx <- getUsefulStatisticsByROI(xxx,roiList)

# here I need to combine a few things

# for each ROI I would be able to tell when/in which frame the max occurs
sapply(roiSNR,function(x) {which.max(apply(x[,,90:150], 3, mean))+90})

ggplot(stats.xxx,aes(xpos,ypos,fill=snr)) + 
  geom_raster(interpolate = F) + 
  scale_fill_gradientn(colors = jet(20)) +
  coord_fixed() + scale_y_reverse()


matFile=myFile
trials <- list()
count=1
for (stimulus in 1:nrow(protocolList[['sabineProtocolSimple']]$presentationMatrix)) {
  for (block in 1:ncol(protocolList[['sabineProtocolSimple']]$presentationMatrix)) {
    # print(protocolList[['sabineProtocolSimple']]$presentationMatrix[stimulus,block])
    tmpdf <- subset(
      protocolList[['sabineProtocolSimple']]$stimulationSections,
      section==protocolList[['sabineProtocolSimple']]$presentationMatrix[stimulus,block]
    )
    description <- subset(
      protocolList[['sabineProtocolSimple']]$stimulationSections,
      section==protocolList[['sabineProtocolSimple']]$presentationMatrix[stimulus,block] & 
        ( description!="background" & description!="settle" )
    )$description
    
    backgroundLengthInMilliseconds <- tmpdf[tmpdf$description=="background","time"] # in ms
    start <- lsm.transition.frames[count] # in frames
    backgroundSlices <- c(start:(start + backgroundLengthInMilliseconds * 0.1))
    fromStimPresentationToEndOfStimulus <- tmpdf[tmpdf$description==description,"time"] + 
      tmpdf[tmpdf$description=="settle","time"] # in ms
    stimulusPeriod <- sum(tmpdf$time)
    end <- start + (stimulusPeriod * 0.1) + 300
    
    trials[[paste(basename(matFile),block,stimulus,sep=".")]] <- list(
      "matFile"=basename(matFile),
      "block"=block,
      "stimulus"=stimulus,
      "start"=start,
      "end"=end,
      "stimulusPeriod"=stimulusPeriod,
      "backgroundSlices"=backgroundSlices,
      "stimulusDescription"=description,
      "creationDate"=date()
    )
    
    count=count+1
  }
  
}
