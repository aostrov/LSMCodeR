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
for (physioDir in physioDirs[5]){
  lsmLogFile <- dir(file.path(physioDir,"logs"),full=T,rec=F,patt="lsmlog_")
  stimLogFile <- dir(file.path(physioDir,"logs"),full=T,rec=F,patt="stimlog_")
  myFile <- dir(file.path(physioDir),full=T,rec=F,patt=".mat")
  outDirSubDir <- paste(basename(physioDir),"_functionTest/",sep="")
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

 myImage<-imageDataSlice[rangeOfImages,,,,]
 myImage <- aperm(myImage,c(2,3,1))
 myImage <- aperm(myImage,c(2,1,3))
 write.nrrd(myImage,"C:/Users/Aaron/Desktop/raw.nrrd")
 write.nrrd((myImage-407),"C:/Users/Aaron/Desktop/baselineSubtracted.nrrd")
 baselineSubtracted <- myImage-407
 background <- apply(baselineSubtracted[,,c(75:85)],c(1,2),mean)
 write.nrrd(background,"C:/Users/Aaron/Desktop/background.nrrd")
 background.raw <- apply(myImage[,,c(75:85)],c(1,2),mean)
 write.nrrd(background.raw,"C:/Users/Aaron/Desktop/background_raw.nrrd")
 df.raw <- apply(myImage[,,],3,function(x) x-round(background.raw))
 write.nrrd(df.raw,"C:/Users/Aaron/Desktop/df_raw.nrrd")
 df.offset <- apply(myImage[,,],3,function(x) x-round(background))
 write.nrrd(df.offset,"C:/Users/Aaron/Desktop/df_offset.nrrd")
 dim(df.raw) <- c(700,1024,181)
 dim(df.offset) <- c(700,1024,181)
 write.nrrd(df.raw,"C:/Users/Aaron/Desktop/df_raw.nrrd")
 write.nrrd(df.offset,"C:/Users/Aaron/Desktop/df_offset.nrrd")
 dff.raw <- apply(df.raw, 3, function(x) x/background.raw)
 dim(dff.raw) <- c(700,1024,181)
 dff.offset <- apply(df.offset, 3, function(x) x/background)
 write.nrrd(dff.offset,"C:/Users/Aaron/Desktop/dff_offset.nrrd")
 dim(dff.offset) <- c(700,1024,181)
 write.nrrd(dff.offset,"C:/Users/Aaron/Desktop/dff_offset.nrrd")
 write.nrrd(dff.raw,"C:/Users/Aaron/Desktop/dff_raw.nrrd")
 makeDFFwithBaselineSubtraction(myImage,xyzDimOrder = c(1,2,3),backgroundSlices = c(75:85))
 myImage.dff <- makeDFFwithBaselineSubtraction(myImage,xyzDimOrder = c(1,2,3),backgroundSlices = c(75:85))
 write.nrrd(myImage.dff,"C:/Users/Aaron/Desktop/myImage.dff.nrrd")
 myImage.dff <- makeDFFwithBaselineSubtraction(myImage,xyzDimOrder = c(1,2,3),backgroundSlices = c(75:85))
 write.nrrd(myImage.dff,"C:/Users/Aaron/Desktop/myImage.dff.nrrd")
 
 
makeDFF<-function(imageStack,backgroundSlices=c(75:85),xyzDimOrder=c(1,2,3)){
 zdim<-dim(imageStack)[xyzDimOrder[3]]
 xdim<-dim(imageStack)[xyzDimOrder[1]]
 ydim<-dim(imageStack)[xyzDimOrder[2]]
 testSliceBackground<-apply(imageStack[,,backgroundSlices],c(xyzDimOrder[1],xyzDimOrder[2]),mean)
 testSliceBackgroundSubtract<-apply(imageStack[,,],3,function(x) x-testSliceBackground)
 dim(testSliceBackgroundSubtract)<-c(xdim,ydim,zdim)
 testSliceDFF<-apply(testSliceBackgroundSubtract, 3, function(x) x/testSliceBackground)
 dim(testSliceDFF)<-c(xdim,ydim,zdim)
 return(testSliceDFF)
}


processSingleStimulus <- function(myList,stimulus=1,block=1,
                                  outputType=c("snr","raw","dff"),writeNRRD=FALSE) {
  # I'll probably want some kind of checking here for existing
  # files, etc
  
  outputType = match.arg(outputType)
  i <- getIndexOfStimulusBlockPair(stimulus,block)
  # print(i)
  
  rangeOfImages <- seq(from=myList[[i]]$start,
                     to=myList[[i]]$end,
                     by=myList[[i]]$timeResampled)
  
  downSampledImage <- apply(
    imageDataSlice[rangeOfImages,,,,],
    1,
    function(x) resizeImage(x,myList[[i]]$resize[1],
                            myList[[i]]$resize[2]))
  
  dim(downSampledImage)<-c(myList[[i]]$resize[1],
                           myList[[i]]$resize[2],
                           length(rangeOfImages))
  
  offsetCorrected <- returnOffsetedImage(downSampledImage,offsetValue=pixelOffset)
  
  if (outputType == 'dff'){
    offsetCorrected <- makeDFFwithBaselineSubtraction(
      offsetCorrected,
      xyzDimOrder = c(1,2,3),
      backgroundSlices=myList[[i]]$backgroundSlices
    )
    
  } else if (outputType=="snr") {
    offsetCorrected <- makeSNRByPixel(offsetCorrected,
                                      backgroundSlices=myList[[i]]$backgroundSlices
    )
  }
  
  if (writeNRRD){
    outfile <- file.path(myList[[i]]$outDir,paste(myList[[i]]$fileBaseName,"nrrd",sep="."))
    # print(paste("Writing to disk: ",outfile,sep=""))
    write.nrrd(offsetCorrected,file=outfile,dtype="short")
  }
  
  cat("+",fill=10)
  invisible(offsetCorrected)
  
  
}

write.nrrd(xxx[(x+(0*w)/2):(x+(1*w)/2),(y+(0*h)/2):(y+(1*h)/2),],
           file="C:/Users/Aaron/Desktop/firstQuarter.nrrd")
write.nrrd(xxx[(x+(1*w)/2):(x+(2*w)/2),(y+(0*h)/2):(y+(1*h)/2),],
           file="C:/Users/Aaron/Desktop/secondQuarter.nrrd")

write.nrrd(xxx[(x+(0*w)/2):(x+(1*w)/2),(y+(1*h)/2):(y+(2*h)/2),],
           file="C:/Users/Aaron/Desktop/thirdQuarter.nrrd")

write.nrrd(xxx[(x+(1*w)/2):(x+(2*w)/2),(y+(1*h)/2):(y+(2*h)/2),],
           file="C:/Users/Aaron/Desktop/fourthQuarter.nrrd")

zz=10
count=0
for (yy in (seq(zz)-1)){
  for (xx in (seq(zz)-1)){
    print(paste("xxx[(x+(",xx,"*w)/",zz,"):(x+(",xx+1,"*w)/",zz,"),(y+(",yy,"*h)/",zz,"):(y+(",yy+1,"*h)/",zz,"),]",sep=""
               ))
    
    # write.nrrd(xxx[(x+(xx*w)/zz):(x+(xx+1*w)/zz),(y+(yy*h)/zz):(y+(yy+1*h)/zz),],
    #            file=paste("C:/Users/Aaron/Desktop/temp/",xx+yy+count,".nrrd",sep=""),
    #            dtype="short")
    SNR_f0 <- mean(xxx[(x+(xx*w)/zz):(x+(xx+1*w)/zz),(y+(yy*h)/zz):(y+(yy+1*h)/zz),75:85])
    print(SNR_f0)
  }
  count=count+10
}

# so instead of writing each thing out, I can try to get a SNR_max-SNR_F0/SNR_f0
SNR_f0 <- mean(xxx[(x+(xx*w)/zz):(x+(xx+1*w)/zz),(y+(yy*h)/zz):(y+(yy+1*h)/zz),75:85])
# and I would want to look frame by frame, averaging the frames and then finding the frame
# with the max signal
SNR_max.1 <- apply(xxx[(x+(xx*w)/zz):(x+(xx+1*w)/zz),(y+(yy*h)/zz):(y+(yy+1*h)/zz),],c(1,2),mean)
