############
##Anatomy##
############

# open the connection to the .mat file
myFile<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-anatomyPost/20181204-gcamp7F-7d-anatomyPost.mat"
file.h5 <- H5File$new(myFile, mode = "r+")
# save a new variable that contains the image data
# this step can be avoided, and might improve performance on very large datasets
imageData<-file.h5[["imagedata"]]
# data is in the form:
# [stacks,?,slices,rows, columns]

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
write.nrrd(aperm(colMeansAverage,c(3,2,1)),file="C:/Users/Aaron/Desktop/colMeansAverge.nrrd")

###############################
## Average of Anatomy Stacks ##
###############################

frameAverageForAnatomyStacks<-function(imageDataArray,dimensionOrder=c(3,2,1),reorderImageSlices=NULL) {
  averageSlices<-colMeans(imageDataArray,dims=1)
  averageSlices<-aperm(averageSlices,dimensionOrder)
  if (!is.null(reorderImageSlices)){
    averageSlices<-averageSlices[,,reorderImageSlices]
  }
  return(averageSlices)
}

doesReorder<-function(testImageArray,reorderImageSlices=c(84:100,1:83)){
  myAnatomyStack<-frameAverageForAnatomyStacks(testImageArray) - 
    frameAverageForAnatomyStacks(testImageArray,reorderImageSlices = reorderImageSlices)
  if (any(myAnatomyStack)) { 
    return(1) 
  }
  return(0)
}

##############################
## max intensity projection ##
##############################

# create an intensity projection of a 3D image 
# supports maximum, minimum, and mean projections,
# as well as standard deviation. Sd projections are comparatively slow
#
# Assumes that the images have been read in as a multidimensional vector
# The default is to assume that x is the 3rd column, and y is the 2nd column,
# but this default can be overwritten.
# Returns a 2D vector of equal length to x*y, where the z dimension has been collapsed.
#
# If a multi-frame image is passed as an array form [frames,?,slices,rows,columns]
# apply() is used to recursively call intensityProjetion.
# This returns a 3D vector that has dimensions x*y*f with the z dimension collapsed.
#
# At the moment this function does not support direct processing of H5D objects,
# these objects must first be coerced to a 4D array before being passed to the function. 
intensityProjection<-function(imageDataArray,rowColumn=c(3,2),projectionType=c("max","min","mean","sd"),fourD=FALSE) {
  if (fourD) {
    print("recursive stage reached")
    dims<-dim(imageDataArray)
    ip<-apply(imageDataArray,1,intensityProjection,rowColumn=rowColumn,projectionType=projectionType)
    dim(ip)<-c(dims[4],dims[3],dims[1])
    return(ip)    
  } else {
    print("deep state reached")
    print(dim(imageDataArray))
    projectionType=match.arg(projectionType)
    apply(imageDataArray, rowColumn, match.fun(projectionType))
  }
}
mip<-intensityProjection(colMeansAverage)
write.nrrd(mip,file="C:/Users/Aaron/Desktop/mip2slices.nrrd")

################
## Physiology ##
################

myFirstSlice<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-SabineBars-1plane-2SP/20181204-gcamp7F-7d-SabineBars-1plane-2SP.mat"
myFirstSliceLog<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-SabineBars-1plane-2SP/logs/lsmlog_acq.xml"
file.h5 <- H5File$new(myFirstSlice, mode = "r+")
imageDataSlice<-file.h5[['imagedata']]

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
startOfStimulations<-min(slicesWithBigNumbers)
# get the salient frames for the second green flash 
slicesWithBigNumbers<-c()
for (i in (33001-2000):33001){
  if (mean(imageDataSlice[i,,,200:400,])>(flashmean+2*flashsd)){
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

endOfFirstGreenFlash<-20
stimulusPeriod<-1600
analysisWindow<-1800
numberOfStimuliInBlock<-4
dryRun=F
downSampleInTime<-10

for (block in 0:4){
  for (stimulus in 0:3){
    print(paste("stimulus bar:",stimulus,"for stimulus block: ",block))
    x<-((endOfFirstGreenFlash+1)+(stimulusPeriod*stimulus))+
      (stimulusPeriod*numberOfStimuliInBlock*block)
    y<-((endOfFirstGreenFlash+1)+(stimulusPeriod*stimulus))+
      (stimulusPeriod*numberOfStimuliInBlock*block)+analysisWindow
    print(x)
    print(y)
    if (!dryRun) {    
      rangeOfImages<-seq(from=x,to=y,by=downSampleInTime)
      downSampledImage<-apply(
        imageDataSlice[rangeOfImages,,,,],
        1,
        function(x) resizeImage(x,350,256))
      dim(downSampledImage)<-c(350,256,length(rangeOfImages))
      write.nrrd(makeDFF(downSampledImage,xyzDimOrder = c(1,2,3),backgroundSlices=c(75:85)),
                 file.path(nrrdFiles,
                           paste("stimulus_bar-",stimulus+1,"-for_stimulus_block-",block+1,".nrrd",sep="")))
    }  
  }
}

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
                                ,"backgroundSlices"=c(75:85),
                                "resize"=c(350,256),
                                "timeResampled"=10)
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
stim2<-c(2,6,10,14,18)

average<-array(data=0,dim = c(350,256,181))
# dim(average)<-c(350,256,181)
for (animal in stim2){
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
average<-average/4
write.nrrd(average,file.path(nrrdFiles,"Average_stim2.nrrd"))
