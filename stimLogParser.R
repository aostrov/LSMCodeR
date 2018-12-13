logdir<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-SabineBars-1plane-2SP/logs/"
logFile <- "stimlog_2018_12_04_11_51_15.txt"

stimLog <- readLines(file.path(logdir,logFile), warn=FALSE)
stimLogCleaning1 <- stimLog[
  grep("@log",stimLog):length(stimLog)][grep("[[:graph:]]* [[:graph:]]* [[:graph:]]*|bar$|green$|red$",
                                                                       stimLog[grep("@log",stimLog):length(stimLog)])]
stimdf <- data.frame(matrix(stimLogCleaning1,ncol=2,byrow = T))
stimdf$X2<-as.numeric(sub(" [[:graph:]]* [[:graph:]]*$","",stimdf$X2))
names(stimdf) <- c("shader","seconds")
stimdiff <- diff(as.integer(stimdf$shader))
stimTransitions <- stimdf[grep(-2,stimdiff),]

 # TODO: stimTransitions is missing the first instance of the transition
stimdf[grep(-1,stimdiff),][2,] # it should look something like this, though this returns one frame before what I want 

# the next step here is to find out how much offset there is between the stimlog file (this one) and the LSM log file
# once I know this offset, I can find the transitions in the stimlog file, and then find the frame in the LSM log
# that most closely resembles this time (rounding up!). Then I can define the frames exactly, without having
# to fuck about with guessing frames from seconds elapsed according to the Protocol.