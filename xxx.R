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

