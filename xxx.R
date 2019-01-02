mip<-intensityProjection(colMeansAverage)
write.nrrd(mip,file="C:/Users/Aaron/Desktop/mip2slices.nrrd")

################
## Physiology ##
################

myFirstSlice<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-SabineBars-1plane-2SP/20181204-gcamp7F-7d-SabineBars-1plane-2SP.mat"
file.h5 <- H5File$new(file.path(myNextSlice), mode = "r")
imageDataSlice<-file.h5[['imagedata']]
# file.h5$close()

rois<-list(
  soma_1=list(x=462, y=287, width=22, height=20),
  soma_2=list(x=467, y=307, width=22, height=20),
  soma_3=list(x=500, y=289, width=14, height=10),
  neuropil_1=list(x=528, y=263, width=27, height=21),
  dendrite_1=list(x=571, y=286, width=19, height=11),
  background=list(x=330,y=385,width=22,height=18)
)

# this returns a cropped and downsampled in time array
downSampledInTime<-imageDataSlice[seq(from=1,to=33000,by=10)
                                  ,,,c(rois[[2]]$y:(rois[[2]]$y+rois[[2]]$height))
                                  ,c(rois[[2]]$x:(rois[[2]]$x+rois[[2]]$width))]
write.nrrd(aperm(downSampledInTime,c(3,2,1)),file="C:/Users/Aaron/Desktop/downSampledTime1per10_ROI-dendrite_1.nrrd")

# writing out and working with ROIs
soma_2<-imageDataSlice[,,,
                       c(rois[[2]]$y:(rois[[2]]$y+rois[[2]]$height)),
                       c(rois[[2]]$x:(rois[[2]]$x+rois[[2]]$width))]
background<-imageDataSlice[seq(from=10000,to=13000,by=10)
                           ,,,c(rois[['background']]$y:(rois[['background']]$y+rois[['background']]$height))
                           ,c(rois[['background']]$x:(rois[['background']]$x+rois[['background']]$width))]
avgBackground<-apply(background,c(3,2),mean)
avgBackgroundCheating<-c(rep(403,21*23))
dim(avgBackgroundCheating)<-c(21,23)
soma_2_backgroundSubtract<-apply(soma_2[,,],1,function(x) x-avgBackgroundCheating)
dim(soma_2_backgroundSubtract)<-c(21,23,33001)
soma_2_dff<-apply(soma_2_backgroundSubtract[,,],3,function(x) x/avgBackgroundCheating)


# get full time range
firstFrame=   10.012 # in ms
lastFrame=    630473.726 # in ms
totalFrames=  63000
timePerFrame= (lastFrame-firstFrame)/totalFrames # ms per frame

getFrameFromSeconds<-function(numSeconds,timePerFrame){
  numSeconds*1000*(1/timePerFrame)
}

nrrdFiles<-"C:/Users/Aaron/Desktop/nrrdFiles/"
resizedImage<-apply(imageDataSlice[1621:(1620+1600),,,,],1,function(x) resizeImage(x,350,256))
dim(resizedImage)<-c(350,256,length(1621:(1620+1600)))
testdff<-makeDFF(resizedImage,xyzDimOrder=c(1,2,3))
write.nrrd(testdff,file.path(nrrdFiles,"testdff2.nrrd"))

##################################################################


###############
## Log Files ##
###############

logFile<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-SabineBars-1plane-2SP/logs/lsmlog_acq.xml"
logFileMetaData<-readLogFileMetaData(logFile)
logFileParsed<-readLogFileData(logFile)

# read in the table
# I will need to redo this importing to 
# try to clean up the log directly in R
crudely <- read.table("C:/Users/Aaron/Documents/R/LSMCodeR/stuff/crudelyParsedLog.txt",sep=" ")
# the stimulus shader can be read as a factor
# and the first derivative of this can tell
# us when the transitions happen
crudelyDiffed <- diff(as.integer(crudely$V1))
crudeTransitionsRedToBar <- crudely[grep(-2,crudelyDiffed),] # note that this is missing the first transition which goes from green to bar

startOfStimulations<-12
testVect <- seq(from=1,to=5, by=0.5)
filePath<-""

# I need to think of a way to tie all of the various scripts together, 
# maybe getting an input folder and being able to pass that input to 
# various sourced scripts so they can get the logs and other useful 
# things. This will require breaking some other things, but that can 
# be dealt with later.

# TODO: 
# - Get an inout dir or vector of dirs. 
# - Read the HDF5 file to find the green flashes.
# - Read the logs to get the needed offset
# - Read the logs to get the frames that start each
# - Process the stimulus presentations

file.h5 <- H5File$new(myFile, mode = "r+")
imageDataSlice<-file.h5[["imagedata"]]

#################
## Green Flash ##
#################
# - Read the HDF5 file to find the green flashes.
flashmean <- mean(imageDataSlice[1:100,,,10:20,10:20])
flashsd <- sd(imageDataSlice[1:100,,,10:20,10:20])

slicesWithBigNumbers.begin <- c()
for (i in 1:100){
  if (mean(imageDataSlice[i,,,10:20,10:20])>(flashmean+2*flashsd)){
    slicesWithBigNumbers.begin <- c(slicesWithBigNumbers.begin,i)
  }
}

startOfStimulations <- max(slicesWithBigNumbers.begin) # last slice of starting green flash
# get the salient frames for the second green flash 

slicesWithBigNumbers.end <- c()
for (i in (33001-2000):33001){
  if (mean(imageDataSlice[i,,,10:20,10:20])>(flashmean+2*flashsd)){
    slicesWithBigNumbers.end <- c(slicesWithBigNumbers.end,i)
  }
}
endOfStimulations <- min(slicesWithBigNumbers.end) # first slice of ending flash

endOfFirstGreenFlash <- startOfStimulations + 1

# - Read the logs to get the needed offset
# source(file.path(LSMCodeRConfig$srcdir,'stimLogParser.R')) # first the stim log
# source(file.path(LSMCodeRConfig$srcdir,'LSMLogParser.R')) # next the LSM log
# # there now exists an array `lsm.transition.frames` with the starting frame of each new stimuls set
# # Previously in physiologyScript.R I would create a list `presentationList2` that could in principle be used to make a best estimate/guess of these times. Now I have them explicitly, so I could set up presentation list much more easily now.

presentationList2<-list()
count=1
for (block in 0:4){
  for (stimulus in 0:3){
    print(paste("stimulus bar:",stimulus,"for stimulus block: ",block))
    x <- lsm.transition.frames[count]
    y <- x + (stimulusPeriod + 200)
    presentationList2[[count]]<-list("block"=block+1,
                                    "stimulus"=stimulus+1,
                                    "start"=x,
                                    "end"=y,
                                    "stimulusPeriod"=stimulusPeriod,
                                    "analysisWindow"=analysisWindow,
                                    "outFile"=file.path(outDir,
                                                        paste("stimulus_bar-",
                                                              stimulus+1,
                                                              "-for_stimulus_block-",
                                                              block+1,
                                                              ".nrrd",
                                                              sep="")
                                    )
                                    ,"backgroundSlices"=backgroundSlices,
                                    "resize"=c(350,256),
                                    "timeResampled"=downSampleInTime)
    print(x)
    print(y)
    count=count+1
  }
}

# outputType <- c('raw','dff')
outputType <- 'dff'
count2 = 1
for (block in 0:4){
  for (stimulus in 0:3){
    if (file.exists(presentationList2[[count2]]$outFile)) {
      print(paste(presentationList2[[count2]]$outFile,
                  "already exists. Skipping."))
      next()
    }
    print(paste("stimulus bar:",stimulus,"for stimulus block: ",block))
    x <- presentationList2[[count2]]$start
    y <- presentationList2[[count2]]$end
    print(x)
    print(y)
    if (!dryRun) {    
      rangeOfImages<-seq(from=x,to=y,by=presentationList2[[count2]]$timeResampled)
      downSampledImage<-apply(
        imageDataSlice[rangeOfImages,,,,],
        1,
        function(x) resizeImage(x,350,256))
      dim(downSampledImage)<-c(350,256,length(rangeOfImages))
      write.nrrd(
        ifelse(
          outputType == 'dff',
            makeDFF(downSampledImage,
                  xyzDimOrder = c(1,2,3),
                  backgroundSlices=presentationList2[[count2]]$backgroundSlices),
            downSampledImage),
        presentationList2[[count2]]$outFile)
    }  
    count2 <- count2 + 1
  }
}
