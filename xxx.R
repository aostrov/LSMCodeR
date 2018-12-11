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
