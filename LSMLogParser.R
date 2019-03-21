# LSMLogParser.R

# lsmLogFile <- file.path(LSMCodeRConfig$logDir,
                     # "20181204-gcamp7F-7d-SabineBars-1plane-2SP",
                     # "lsmlog_acq.xml")
# lsmLogFile <- file.path(logdir,"lsmlog_acq.xml")
lsmLogFile <- dir(file.path(dirname(myFile),"logs"),full=T,patt="lsmlog_")
if (length(lsmLogFile)!=1) stop(paste("Check your lsmLogFile for file: ", myfile, sep=""))

# read in some data
logFileParsed <- readLogFileData(lsmLogFile)
# decide if I need to discard some first set of frames
lsmLogFileShort <- logFileParsed[(nrow(logFileParsed) - protocolList$sabineProtocolSimple$framesSkipped):nrow(logFileParsed),]

# logFileMetaData<-readLogFileMetaData(logFile)
# logFileParsed<-readLogFileData(logFile)

endOfStimulations.frame <- (startOfStimulations + 1) # 12 # dummy variable that I'll need to grab from a script later
if (nrow(stimdf[stimdf$shader=='green',])==4) {
  # proceed as normal
  someDodgyFix <- stimdf[stimdf$shader=='green',][2,'seconds']
} else {
  if (diff(round(stimdf[stimdf$shader=='green',"seconds"]))[1] > diff(round(stimdf[stimdf$shader=='green',"seconds"]))[2]) {
    someDodgyFix <- stimdf[stimdf$shader=='green',"seconds"][1]
  } else {
    someDodgyFix <- stimdf[stimdf$shader=='green',"seconds"][2]
  }
}


lsm2stim.offset.ms <- lsmLogFileShort[endOfStimulations.frame,'time'] - 
  (someDodgyFix*1000) # get an offset based on the end of the green flashes

transitions.stimdf.bar <- stimdf[newRows,]
transitions.stimdf.bar$correctedMilliseconds <- transitions.stimdf.bar$seconds*1000 + lsm2stim.offset.ms

# get the frames in the lsmlog at which a transition to the bar shader occurs
lsm.transition.frames <- c()
for (i in transitions.stimdf.bar$correctedMilliseconds) {
  lsm.transition.frames <- c(lsm.transition.frames,findInterval(i,lsmLogFileShort$time))
}
# not sure if I want the vector of frames as is, or if i want to offset it by +1 to account for getting the next frame...

count=0
slice.identity <- data.frame(slice=integer(0),z_plane=integer(0),time=integer(0))
for (time in 1:1650){ # FIXME
  for (z_plane in 1:20){ # FIXME
    count=count+1
    slice.identity <- rbind(slice.identity,c(count,z_plane,time))
  }
}
colnames(slice.identity) <- c("slice","z_plane","time")
slice.transitions <- slice.identity[lsm.transition.frames,]
