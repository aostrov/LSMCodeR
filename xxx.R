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

average<-array(data=0,dim = c(350,256,1801))
# Now I want just one average that isn't downsampled in time
for (stimulus in seq(from=4,to=20,by=4)){
  start <- (presentationList2[[stimulus]]$start)
  end <- (presentationList2[[stimulus]]$end)
  rangeOfImages<-seq(from=start,to=end)
  print("downsampling...")
  downSampledImage<-apply(
    imageDataSlice[rangeOfImages,,,,],
    1,
    function(x) resizeImage(x,350,256))
  dim(downSampledImage)<-c(350,256,length(rangeOfImages))
  print("making the average ...")
  print(dim(downSampledImage))
  if (outputType == 'dff'){
    average <- average + makeDFF(downSampledImage,
                                 xyzDimOrder = c(1,2,3),
                                 backgroundSlices=c(1:100))
    
  } else {
    average <- average + downSampledImage
  }
}
average<-average/length(seq(from=4,to=20,by=4))
write.nrrd(average,file.path(outDir,outDirSubDir,paste("Average_dff_fullTime_stim","4",".nrrd",sep="")))

outDir<-"F:/Imaging/GCaMP7_tests/outputNRRDs/"
physioDirs <- dir("F:/Imaging/GCaMP7_tests/20181204-g7",patt='SP',full=T)
physiologyFiles <- dir("F:/Imaging/GCaMP7_tests/20181204-g7/",patt="[A-Z]{4}-[[:graph:]]*.mat",full=TRUE,rec=T)

physiologyFiles <- dir("F:/Imaging/GCaMP7_tests/20181204-g7/",patt="[A-Z]{4}-[[:graph:]]*SP",full=T,rec=FALSE)
special <- "examples"
outputTypes=c("snr","raw","dff")
outputType=("raw")
for (physioDir in physiologyFiles){
  lsmLogFile <- dir(file.path(physioDir,"logs"),full=T,rec=F,patt="lsmlog_")
  stimLogFile <- dir(file.path(physioDir,"logs"),full=T,rec=F,patt="stimlog_")
  myFile <- dir(physioDir,full=T,rec=F,patt=".mat")
  outDirSubDir <- paste(basename(physioDir),"_",outputType,special,"/",sep="")
  # source(file.path(LSMCodeRConfig$srcdir,"parseImageForGreenFlashAndLSMandStimLogs.R"))
  # source(file.path(LSMCodeRConfig$srcdir,"roiBasedAnalysis.R"))
  source(file.path(LSMCodeRConfig$srcdir,"physiologyScript.R"))
}


matFile <- "F:/Imaging/GCaMP7_tests/20190109-jGCaMP7fEF05-testFish1-SP/20190109-jGCaMP7fEF05-testFish1-SP.mat"
file.h5 <- H5File$new(file.path(myFile), mode = "r")
imageDataSlice<-file.h5[['imagedata']]

imageDataSlice$dims
# frames channels slices rows  columns
# 33001     1        1    512   700

rangeOfImages <- c(15000:20000)
downSampledImage<-apply(
  imageDataSlice[rangeOfImages,,,,],
  1,
  function(x) resizeImage(x,350,256))
dim(downSampledImage)<-c(350,256,length(rangeOfImages))

# testSliceBackground <- apply(imageStack[,,backgroundSlices],c(xyzDimOrder[1],xyzDimOrder[2]),mean)



xxx <- makeDFF2(downSampledImage,
                xyzDimOrder = c(1,2,3),
                backgroundSlices = 100:200)

file.h5 <- H5File$new(myFile, mode = "r+")
# save a new variable that contains the image data
# this step can be avoided, and might improve performance on very large datasets
imageDataSlice<-file.h5[["imagedata"]]





x = 348
y = 884
width = 37
height = 33

imageStack <- aperm(imageDataSlice[13:1813,,,(y:(y+height)),(x:(x+width))],c(3,2,1))
outsideFishROIBackground <- floor(mean(apply(imageDataSlice[200:2200,,,1:10,1:10],3,mean)))
myImage.cropped.baselineCorrected <- (imageStack - outsideFishROIBackground)
# dim(myImage.cropped.baselineCorrected)<-c(xdim,ydim,zdim)
roi.average <- apply(myImage.cropped.baselineCorrected,3,mean)


write.nrrd(xxx[(x+(0*w)/2):(x+(1*w)/2),(y+(0*h)/2):(y+(1*h)/2),],
           file="C:/Users/Aaron/Desktop/firstQuarter.nrrd")
write.nrrd(xxx[(x+(1*w)/2):(x+(2*w)/2),(y+(0*h)/2):(y+(1*h)/2),],
           file="C:/Users/Aaron/Desktop/secondQuarter.nrrd")

write.nrrd(xxx[(x+(0*w)/2):(x+(1*w)/2),(y+(1*h)/2):(y+(2*h)/2),],
           file="C:/Users/Aaron/Desktop/thirdQuarter.nrrd")

write.nrrd(xxx[(x+(1*w)/2):(x+(2*w)/2),(y+(1*h)/2):(y+(2*h)/2),],
           file="C:/Users/Aaron/Desktop/fourthQuarter.nrrd")

zz=10


stimDir <- dir(file.path(outDir,"20190117-fish_03-0-SabineBars_gen_A-SP_functionTest"),full=T)
fullDF <- data.frame()
for (stimFile in stimDir) {
  nrrd <- read.nrrd(stimFile)
  tempDF <- getSNRforSubROIs(nrrd,numSubunits = 40,x=53,y=44, w=280,h=280,dryRun = F)
  tempDF$file <- basename(stimFile)
  fullDF <- rbind(fullDF, tempDF)
}

ggplot(fullDF,aes(X,Y,fill=dSNR)) + 
  geom_raster(interpolate = T) + 
  facet_wrap(~file) + 
  scale_fill_gradientn(colors = jet(20))


by(fullDF[,"dSNR"],fullDF[,"file"],which.max)

# Maybe for each ROI I would look at the data and for each time point get
# mean, max, and median values



#######################
## multiplane images ##
#######################

file.h5 <- H5File$new("/Volumes/TranscendJD/Work/AAFA-gen-A-laser3-SabineSimple/AAFA-gen-A-laser3-SabineSimple.mat", 
                      mode = "r+")
imageDataSlice<-file.h5[["imagedata"]]
test.matrix <- imageDataSlice[2,,12,,]
test.matrix.melt <- melt(data = test.matrix)
ggplot(test.matrix.melt,aes(Var2,Var1,fill=value)) + 
  geom_raster(interpolate = F) + 
  scale_fill_gradientn(colors = jet(2000)) +
  coord_fixed() + scale_y_reverse()

write.nrrd(aperm(imageDataSlice[1200:1600,,13,,],c(3,2,1)),file = "C:/Users/Aaron/Desktop/multiplaneNrrd.nrrd",dtype = "short")

flashMean <- mean(imageDataSlice[1:100,,,1:200,1:200])
flashSD <- sd(imageDataSlice[1:100,,,1:200,1:200])
flashSlices <- c()
for (i in 1:20){
  flashSlices <- c(flashSlices,mean(imageDataSlice[1:10,,i,1:200,1:200]))
}
max.flash <- which.max(flashSlices)

slicesWithBigNumbers.begin <- c()
for (i in 1:20){
  if (mean(imageDataSlice[1:10,,i,1:200,1:200])>(flashMean+2*flashSD)){
    slicesWithBigNumbers.begin <- c(slicesWithBigNumbers.begin,i)
  }
}
startOfStimulations <- max(slicesWithBigNumbers.begin)
endOfFirstGreenFlash <- startOfStimulations + 1


count=0
slice.identity <- data.frame(slice=integer(0),z_plane=integer(0),time=integer(0))
for (time in 1:1650){
  for (z_plane in 1:20){
    count=count+1
    slice.identity <- rbind(slice.identity,c(count,z_plane,time))
  }
}
colnames(slice.identity) <- c("slice","z_plane","time")
slice.transitions <- slice.identity[lsm.transition.frames,]

makeTrial <- function(matFile,stimProtocol="sabineProtocolSimple",analysisWindow=300) {
  trials <- list()
  count=1
  for (stimulus in 1:nrow(protocolList[[stimProtocol]]$presentationMatrix)) {
    for (block in 1:ncol(protocolList[[stimProtocol]]$presentationMatrix)) {
      for (plane in seq(imageDataSlice.dims[['z']])) {
        tmpdf <- subset(
          protocolList[[stimProtocol]]$stimulationSections,
          section==protocolList[[stimProtocol]]$presentationMatrix[stimulus,block]
        )
        description <- subset(
          protocolList[[stimProtocol]]$stimulationSections,
          section==protocolList[[stimProtocol]]$presentationMatrix[stimulus,block] & 
            ( description!="background" & description!="settle" )
        )$description
        
        backgroundLengthInMilliseconds <- tmpdf[tmpdf$description=="background","time"] # in ms
        
        if ( imageDataSlice.dims[['z']] > 1 ) {
          start <- slice.transitions[count,"time"] # in frames
        } else {
          start <- lsm.transition.frames[count]
        }
        
        backgroundSlices <- c( 
          start : 
            (start + (backgroundLengthInMilliseconds * 0.1) / imageDataSlice.dims[['z']] ) 
          )
        stimulusPeriod <- sum(tmpdf$time) / imageDataSlice.dims[['z']] # in ms
        end <- start + 
          (stimulusPeriod * 0.1) + 
          ( analysisWindow / imageDataSlice.dims[['z']] )
        
     }

      trials[[paste(basename(matFile),stimulus,block,sep=".")]] <- list(
        "matFile"=basename(matFile),
        "stimulusProtocol"=stimProtocol,
        "block"=block,
        "stimulus"=stimulus,
        "plane"=plane,
        "start"=start,
        "end"=end,
        "stimulusPeriodPerSlice"=stimulusPeriod,
        "backgroundSlices"=backgroundSlices,
        "stimulusDescription"=description,
        "creationDate"=date()
      )

      count=count+1
    }
    
  }
  attr(trials,"imageDimensions") <- imageDataSlice.dims
  return(trials)
}


offsets <- c()
for (offset in 1:length(flashSlices)) {
  if (offset < max.flash) {
    offsets <- c(offsets,0)
  } else {
    offsets <- c(offsets,1)
  }
}


count=1
for (i in 1:20){
  write.nrrd(aperm(imageDataSlice[100,,i,,],c(2,1)),file=paste("~/Desktop/nrrdOrder/roiSelection-",i,".nrrd",sep=""),dtype = "short")
  count=count+1
}

reorderImageSlices <- c(2,1,20:5)
file.h5$close_all()
imageDataSlice$close()

someVar <- c()
for (i in ((1:20)-1) ){
  someVar <- c(someVar,length(seq(from=(13+i),to=(13+1800),by=20)))
}
first <- seq(from=(13+1),to=(13+1800),by=20)

