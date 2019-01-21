###################
## Magic Numbers ##
###################
# outDir
# outDirSubDir<-"20181204-gcamp7F-7d-SabineBars-1planeSP"
# outDir<-"F:/Imaging/GCaMP7_tests/outputNRRDs/"
# myFile<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-SabineBars-1plane-Plane41SP/20181204-gcamp7F-7d-SabineBars-1plane-Plane41SP.mat"

if (!dir.exists(file.path(outDir,outDirSubDir))){
  dir.create(file.path(outDir,outDirSubDir),rec=T)
} else {
  print("Dir exists")
}

framesSkipped <- 30000
stimulusPeriod<-1600
analysisWindow<-1800
numberOfStimuliInBlock<-4
dryRun=F
downSampleInTime<-10
backgroundSlices=c(75:85)
resizeFactor <- 2
pixelOffset <- 398



#################
# file handling #
#################

# open the connection to the .mat file
file.h5 <- H5File$new(myFile, mode = "r+")
# save a new variable that contains the image data
# this step can be avoided, and might improve performance on very large datasets
imageDataSlice<-file.h5[["imagedata"]]
# file.h5$close()
# data is in the form:
# [stacks,(channels),slices,rows, columns]
imageDataSlice.dims <- imageDataSlice$dims
# frames channels slices rows  columns
# 33001     1        1    512   700

#################
## Green Flash ##
#################
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
for (i in (33000-2000):33000){
  if (mean(imageDataSlice[i,,,10:20,10:20])>(flashmean+2*flashsd)){
    slicesWithBigNumbers.end <- c(slicesWithBigNumbers.end,i)
  }
}
endOfStimulations <- min(slicesWithBigNumbers.end) # first slice of ending flash

endOfFirstGreenFlash <- startOfStimulations + 1


source(file.path(LSMCodeRConfig$srcdir,"stimLogParser.R"))
source(file.path(LSMCodeRConfig$srcdir,"LSMLogParser.R"))

##################
## Doing Things ##
##################

# Make a list of the different stimuli presentations
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
                                     "numberOfSlices"=(((y-x)/downSampleInTime)+1),
                                     "outDir"=file.path(outDir,outDirSubDir),
                                     "fileBaseName"=paste("stimulus_bar-",
                                                 stimulus+1,
                                                 "-for_stimulus_block-",
                                                 block+1,
                                                 sep=""),
                                     "backgroundSlices"=backgroundSlices,
                                     "resize"=c(imageDataSlice.dims[5]/resizeFactor,imageDataSlice.dims[4]/resizeFactor),
                                     "timeResampled"=downSampleInTime,
                                     "creationDate"=date())
    print(x)
    print(y)
    count=count+1
  }
}

# Write individual files for each stimulus presentation
# outputType <- c('raw','dff','snr')
outputType <- 'snr'
blockCount <- 0
count2 = 1
# make a fragile list to help speed up making of averages
averageList <- list()
for (block in 0:4){
  for (stimulus in 0:3){
    if (blockCount==0) {
      averageList[[stimulus+1]] <- array(data=0,
                                         dim = c(presentationList2[[count2]]$resize[1],
                                                 presentationList2[[count2]]$resize[2],
                                                 presentationList2[[count2]]$numberOfSlices)
                                         )
    }
    if (file.exists(file.path(presentationList2[[count2]]$outDir,paste(presentationList2[[count2]]$fileBaseName,"_dff.nrrd",sep="")))) {
      print(paste(presentationList2[[count2]]$fileBaseName,"_dff.nrrd",
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
        function(x) resizeImage(x,presentationList2[[count2]]$resize[1],presentationList2[[count2]]$resize[2]))
      dim(downSampledImage)<-c(presentationList2[[count2]]$resize[1],presentationList2[[count2]]$resize[2],length(rangeOfImages))
      if (outputType == 'dff'){
        dffImage <- makeDFFwithBaselineSubtraction(
          downSampledImage,
          xyzDimOrder = c(1,2,3),
          backgroundSlices=presentationList2[[count2]]$backgroundSlices
          )
        averageList[[(stimulus+1)]] <- averageList[[(stimulus+1)]] + dffImage
        write.nrrd(
          dffImage,
          file.path(presentationList2[[count2]]$outDir,
                    paste(presentationList2[[count2]]$fileBaseName,"_dff.nrrd",sep="")
          )
        )
        
      } else if (outputType=="snr") {
        offsetCorrected <- returnOffsetedImage(downSampledImage,offsetValue=pixelOffset)
        write.nrrd(makeSNRByPixel(offsetCorrected,
                       backgroundSlices=presentationList2[[count2]]$backgroundSlices),
                   file=file.path(presentationList2[[count2]]$outDir,
                                  paste(presentationList2[[count2]]$fileBaseName,"_snr.nrrd",sep="")),
                                  dtype="short"
                   )
       } else {
        write.nrrd(
          offsetCorrected,
          file.path(presentationList2[[count2]]$outDir,
                    paste(presentationList2[[count2]]$fileBaseName,".nrrd",sep="")
          )
        )
        
        }
      }  
    count2 <- count2 + 1
  }
  blockCount=blockCount+1
}
mapply(function(x, i) {
  stimAvg <- x/5
  write.nrrd(stimAvg,file.path(outDir,outDirSubDir,paste("Average_dff_stim",i,".nrrd",sep="")))
  }, averageList, c(1:4))


file.h5$close()
imageDataSlice$close()
