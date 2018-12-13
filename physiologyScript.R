###################
## Magic Numbers ##
###################
# outDir
outDirSubDir<-"20181204-gcamp7F-7d-SabineBars-1plane-Plane41SP"
outDir<-"C:/Users/Aaron/Desktop/nrrdFiles/"
myFile<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-SabineBars-1plane-Plane41SP/20181204-gcamp7F-7d-SabineBars-1plane-Plane41SP.mat"

if (!dir.exists(file.path(outDir,outDirSubDir))){
  dir.create(file.path(outDir,outDirSubDir))
} else {
  print("Dir exists")
}

stimulusPeriod<-1600
analysisWindow<-1800
numberOfStimuliInBlock<-4
dryRun=F
downSampleInTime<-10
backgroundSlices=c(75:85)



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

#################
## Green Flash ##
#################

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

endOfFirstGreenFlash<-startOfStimulations + 1

##################
## Doing Things ##
##################

# Make a list of the different stimuli presentations
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

# Write individual files for each stimulus presentation
for (block in 0:4){
  for (stimulus in 0:3){
    if (file.exists(file.path(outDir,outDirSubDir,
                              paste("stimulus_bar-",
                                    stimulus+1,
                                    "-for_stimulus_block-",
                                    block+1,
                                    ".nrrd",
                                    sep="")))) {
      print(paste("stimulus_bar-",
                  stimulus+1,
                  "-for_stimulus_block-",
                  block+1,
                  ".nrrd already exists. Skipping.",
                  sep=""))
      next()
    }
    print(paste("stimulus bar:",stimulus,"for stimulus block: ",block))
    x<-((endOfFirstGreenFlash+1)+(stimulusPeriod*stimulus))+
      (stimulusPeriod*numberOfStimuliInBlock*block)
    y<-((endOfFirstGreenFlash+1)+(stimulusPeriod*stimulus))+
      (stimulusPeriod*numberOfStimuliInBlock*block)+analysisWindow
    print(x)
    print(y)
    file.exists(file.path(outDir,outDirSubDir,
                          paste("stimulus_bar-",stimulus+1,"-for_stimulus_block-",block+1,".nrrd",sep="")))
    if (!dryRun) {    
      rangeOfImages<-seq(from=x,to=y,by=downSampleInTime)
      downSampledImage<-apply(
        imageDataSlice[rangeOfImages,,,,],
        1,
        function(x) resizeImage(x,350,256))
      dim(downSampledImage)<-c(350,256,length(rangeOfImages))
      write.nrrd(makeDFF(downSampledImage,xyzDimOrder = c(1,2,3),backgroundSlices=backgroundSlices),
                 file.path(outDir,outDirSubDir,
                           paste("stimulus_bar-",stimulus+1,"-for_stimulus_block-",block+1,".nrrd",sep="")))
    }  
  }
}

# Get averages
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
  for (stimulus in seq(from=i,to=20,by=4)){
    start<-presentationList[[stimulus]]$start
    end<-presentationList[[stimulus]]$end
    rangeOfImages<-seq(from=start,to=end,by=presentationList[[stimulus]]$timeResampled)
    print("downsampling...\n")
    print(dim(average))
    downSampledImage<-apply(
      imageDataSlice[rangeOfImages,,,,],
      1,
      function(x) resizeImage(x,350,256))
    dim(downSampledImage)<-c(350,256,length(rangeOfImages))
    print("making the DFF...\n")
    print(dim(downSampledImage))
    average<-average+makeDFF(downSampledImage,xyzDimOrder = c(1,2,3),backgroundSlices=backgroundSlices)
    
  }
  average<-average/length(seq(from=i,to=20,by=4))
  write.nrrd(average,file.path(outDir,outDirSubDir,paste("Average_stim",i,".nrrd",sep="")))
  
}

file.h5$close()