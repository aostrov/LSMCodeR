############
##Anatomy##
############

# open the connection to the .mat file
myFile<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-anatomyPost/20181204-gcamp7F-7d-anatomyPost.mat"
file.h5 <- H5File$new(myFile, mode = "r+")
# save a new variable that contains the image data
# this step can be avoided, and might improve performance on very large datasets
imageData<-file.h5[["imagedata"]]
# file.h5$close()
# data is in the form:
# [stacks,(channels),slices,rows, columns]

bigVector<-vector(length=imageData$dims[3]*imageData$dims[4]*imageData$dims[5])
bigVector<-as.integer(bigVector)
dim(bigVector)<-c(700,512,100)
for (i in seq(imageData$dims[1])) {
  sliceForElementInLoop<-aperm(imageData[i,,,,],c(3,2,1))
  bigVector<- bigVector + sliceForElementInLoop
}
averageImage<-bigVector/imageData$dims[1]
write.nrrd(averageImage,file="C:/Users/Aaron/Desktop/averageAnatomy.nrrd")
############
# There is some weirdness with having to discard some of the slices
# because of the way the objective moves during acquisition
# In this case, I can either reorder, or remove slices 83 or 84
# until the end; either moving them to the beginning or removing
# them all together
reorderedImage<-averageImage[,,c(84:100,1:83)]
testFileReordered<-"C:/Users/Aaron/Desktop/averageAnatomyReordered.nrrd"
write.nrrd(reorderedImage,testFileReordered)
# Even more concise:
colMeansAverage<-colMeans(imageData[,,,,],dims=1)
write.nrrd(aperm(round(colMeansAverage),c(3,2,1)),
           file.path("C:/Users/Aaron/Desktop/nrrdFiles/anatomyAverageReordered.nrrd"),
                     dtype = 'short')
# Better still
write.nrrd(frameAverageForAnatomyStacks(imageData[,,,,],reorderImageSlices = c(84:100,1:83),roundOutput = T),
           file.path("C:/Users/Aaron/Desktop/nrrdFiles/anatomyAverageReordered.nrrd"),
           dtype = 'short')



mip<-intensityProjection(colMeansAverage)
write.nrrd(mip,file="C:/Users/Aaron/Desktop/mip2slices.nrrd")

################
## Physiology ##
################

myFirstSlice<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-SabineBars-1plane-2SP/20181204-gcamp7F-7d-SabineBars-1plane-2SP.mat"
myFirstSliceLog<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-SabineBars-1plane-2SP/logs/lsmlog_acq.xml"
Plane41SP<-"F:\\Imaging\\GCaMP7_tests\\20181204-g7\\20181204-gcamp7F-7d-SabineBars-1plane-Plane41SP"
plane41SPFile<-"20181204-gcamp7F-7d-SabineBars-1plane-Plane41SP.mat"

myNextSlice<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-SabineBars-1planeSP/20181204-gcamp7F-7d-SabineBars-1planeSP.mat"

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

flashmean<-mean(imageDataSlice[1:100,,,10:20,10:20])
flashsd<-sd(imageDataSlice[1:100,,,10:20,10:20])

slicesWithBigNumbers<-c()
for (i in 1:100){
  if (mean(imageDataSlice[i,,,10:20,10:20])>(flashmean+2*flashsd)){
    slicesWithBigNumbers<-c(slicesWithBigNumbers,i)
  }
}
startOfStimulations<-max(slicesWithBigNumbers)
# get the salient frames for the second green flash 
slicesWithBigNumbers<-c()
for (i in (33001-2000):33001){
  if (mean(imageDataSlice[i,,,10:20,10:20])>(flashmean+2*flashsd)){
    slicesWithBigNumbers<-c(slicesWithBigNumbers,i)
  }
}
endOfStimulations<-min(slicesWithBigNumbers)

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

endOfFirstGreenFlash<-startOfStimulations + 1
stimulusPeriod<-1600
analysisWindow<-1800
numberOfStimuliInBlock<-4
dryRun=F
downSampleInTime<-10
backgroundSlices=c(75:85)

for (block in 0:4){
  for (stimulus in 0:3){
    print(paste("stimulus bar:",stimulus,"for stimulus block: ",block))
    x<-((endOfFirstGreenFlash+1)+(stimulusPeriod*stimulus))+
      (stimulusPeriod*numberOfStimuliInBlock*block)
    y<-((endOfFirstGreenFlash+1)+(stimulusPeriod*stimulus))+
      (stimulusPeriod*numberOfStimuliInBlock*block)+analysisWindow
    print(x)
    print(y)
    file.exists(file.path(nrrdFiles,
                          paste("stimulus_bar-",stimulus+1,"-for_stimulus_block-",block+1,".nrrd",sep="")))
    if (!dryRun) {    
      rangeOfImages<-seq(from=x,to=y,by=downSampleInTime)
      downSampledImage<-apply(
        imageDataSlice[rangeOfImages,,,,],
        1,
        function(x) resizeImage(x,350,256))
      dim(downSampledImage)<-c(350,256,length(rangeOfImages))
      write.nrrd(makeDFF(downSampledImage,xyzDimOrder = c(1,2,3),backgroundSlices=backgroundSlices),
                 file.path(nrrdFiles,
                           paste("stimulus_bar-",stimulus+1,"-for_stimulus_block-",block+1,".nrrd",sep="")))
    }  
  }
}
nrrdFiles <- file.path("C:\\Users\\Aaron\\Desktop\\nrrdFiles\\20181204-gcamp7F-7d-SabineBars-1planeSP")
presentationList<-list()
count=1
for (block in 0:4){
  for (stimulus in 0:3){
    print(paste("stimulus bar:",stimulus,"for stimulus block: ",block))
    x<-((endOfFirstGreenFlash+1)+(stimulusPeriod*stimulus))+
      (stimulusPeriod*numberOfStimuliInBlock*block)
    y<-((endOfFirstGreenFlash+1)+(stimulusPeriod*stimulus))+
      (stimulusPeriod*numberOfStimuliInBlock*block)+analysisWindow
    presentationList[[count]]<-list("block"=block+1,
                                "stimulus"=stimulus+1,
                                "start"=x,
                                "end"=y,
                                "stimulusPeriod"=stimulusPeriod,
                                "analysisWindow"=analysisWindow,
                                "outFile"=file.path(nrrdFiles,
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

# I want to be able to look through the list and get a sublist that can
# be used to perform functions on them, such as averaging.
# I should also decide if I want to store any raw data in this list, 
# although that might be too heavy a solution for memory usage. Better
# might be either lots of reading in and out from the file system nrrds, or 
# to continue using the HDF5 connection for subsetting things.
stim1<-c(1,5,9,13,17)
stim2<-c(2,6,10,14,18)
stim3<-c(3,7,11,15,19)
stim4<-c(4,8,12,16,20)
for(i in 1:4){
  print(seq(from=i,to=20,by=4))
  average<-array(data=0,dim = c(350,256,181))
  # dim(average)<-c(350,256,181)
  for (animal in seq(from=i,to=20,by=4)){
    start<-presentationList[[animal]]$start
    end<-presentationList[[animal]]$end
    rangeOfImages<-seq(from=start,to=end,by=presentationList[[animal]]$timeResampled)
    print("downsampling...\n")
    print(dim(average))
    downSampledImage<-apply(
      imageDataSlice[rangeOfImages,,,,],
      1,
      function(x) resizeImage(x,350,256))
    dim(downSampledImage)<-c(350,256,length(rangeOfImages))
    print("making the DFF...\n")
    print(dim(downSampledImage))
    average<-average+makeDFF(downSampledImage,xyzDimOrder = c(1,2,3),backgroundSlices=c(75:85))
    
  }
  average<-average/length(seq(from=i,to=20,by=4))
  write.nrrd(average,file.path(nrrdFiles,paste("Average_stim",i,".nrrd",sep="")))
  
}

###############
## Log Files ##
###############

logFile<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-SabineBars-1plane-2SP/logs/lsmlog_acq.xml"
logFileMetaData<-readLogFileMetaData(logFile)
logFileParsed<-readLogFile(logFile)
