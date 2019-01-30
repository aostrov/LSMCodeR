init.singleTrial <- function(matFile, outDir, 
                             block, stimulus, start, end,
                             stimulusPeriod, analysisWindow, 
                             startOfStimulusInSeconds = 9, timeForF0Integration = 2, 
                             resizeFactor=2, timeResampled=10, 
                             imageDims) {
  
  
  getBackgroundSlices <- function(startOfStimulusInSeconds = startOfStimulusInSeconds, 
                                  timeForF0Integration = timeForF0Integration, 
                                  timeResampled = timeResampled){
    timePerFrame <- 0.01 * timeResampled # seconds /frame
    # seconds / seconds * frame == frames
    timeForF0IntegrationInFrames <- timeForF0Integration / timePerFrame 
    startOfStimulusInFrames <- startOfStimulusInSeconds / timePerFrame
    backgroundSlices <- ceiling(
      c((startOfStimulusInFrames - timeForF0IntegrationInFrames):startOfStimulusInFrames)
    )
    return(backgroundSlices)
  }
  backgroundSlices <- getBackgroundSlices(startOfStimulusInSeconds = startOfStimulusInSeconds, 
                                          timeForF0Integration = timeForF0Integration, 
                                          timeResampled = timeResampled)
  
  singleTrial <- list("matFile"=basename(matFile),
                      "block"=block+1,
                      "stimulus"=stimulus+1,
                      "start"=start,
                      "end"=end,
                      "stimulusPeriod"=stimulusPeriod,
                      "analysisWindow"=analysisWindow,
                      "numberOfSlices"=(((end-start)/downSampleInTime)+1),
                      "outDir"=file.path(outDir,outDirSubDir),
                      "fileBaseName"=paste("stimulus_bar-",
                                           stimulus+1,
                                           "-for_stimulus_block-",
                                           block+1,
                                           sep=""),
                      "backgroundSlices"=backgroundSlices,
                      "resize"=c(imageDims[5]/resizeFactor,
                                 imageDims[4]/resizeFactor),
                      "timeResampled"=timeResampled,
                      "creationDate"=date())
  return(singleTrial)
}


init.Animal <- function(pathToSpreadsheet) {
  Animal <- read.csv(pathToSpreadsheet)
  # animal_ID, genotype_injection, date_injected, date_run
  return(Animal)
}

init.MatFile <- function(animalID, pathToMatFile, protocolUsed) {
  x=unique(protocolList[[protocolUsed]]$stimulationSections[
    protocolList[[protocolUsed]]$stimulationSections$speed == 0 &
      !is.na(protocolList[[protocolUsed]]$stimulationSections$speed),
    ]$time)
  if (length(x)>1) stop("Your protocol is too complex")
  MatFile <- list(animalID = animalID,
                  matID = basename(pathToMatFile),
                  stimulusProtocol = names(protocolList[protocolUsed]),
                  transitionFrames= lsm.transition.frames
  )
  MatFile
}

init.Trial <- function(matID,) {
  
}