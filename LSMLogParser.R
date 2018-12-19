# LSMLogParser.R

lsmLogFile <- file.path(LSMCodeRConfig$logDir,
                     "20181204-gcamp7F-7d-SabineBars-1plane-2SP",
                     "lsmlog_acq.xml")

# read in some data
logFileParsed<-readLogFileData(lsmLogFile)
# decide if I need to discard some first set of frames

# logFileMetaData<-readLogFileMetaData(logFile)
# logFileParsed<-readLogFileData(logFile)

startOfStimulations<-12 # dummy variable that I'll need to grab from a script later


