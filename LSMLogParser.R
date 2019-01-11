# LSMLogParser.R

# lsmLogFile <- file.path(LSMCodeRConfig$logDir,
                     # "20181204-gcamp7F-7d-SabineBars-1plane-2SP",
                     # "lsmlog_acq.xml")
# lsmLogFile <- file.path(logdir,"lsmlog_acq.xml")
# read in some data
logFileParsed <- readLogFileData(lsmLogFile)
# decide if I need to discard some first set of frames
lsmLogFileShort <- logFileParsed[(nrow(logFileParsed) - framesSkipped):nrow(logFileParsed),]

# logFileMetaData<-readLogFileMetaData(logFile)
# logFileParsed<-readLogFileData(logFile)

endOfStimulations.frame <- (startOfStimulations + 1) # 12 # dummy variable that I'll need to grab from a script later
lsm2stim.offset.ms <- lsmLogFileShort[endOfStimulations.frame,'time'] - 
  (stimdf[stimdf$shader=='green',][2,'seconds']*1000) # get an offset based on the end of the green flashes

transitions.stimdf.bar <- stimdf[newRows,]
transitions.stimdf.bar$correctedMilliseconds <- transitions.stimdf.bar$seconds*1000 + lsm2stim.offset.ms

# get the frames in the lsmlog at which a transition to the bar shader occurs
lsm.transition.frames <- c()
for (i in transitions.stimdf.bar$correctedMilliseconds) {
  lsm.transition.frames <- c(lsm.transition.frames,findInterval(i,lsmLogFileShort$time))
}
# not sure if I want the vector of frames as is, or if i want to offset it by +1 to account for getting the next frame...

