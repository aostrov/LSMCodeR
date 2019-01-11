stimLog <- readLines(stimLogFile, warn=FALSE)
stimLogCleaning1 <- stimLog[
  grep("@log",stimLog):length(stimLog)][grep("[[:graph:]]* [[:graph:]]* [[:graph:]]*|bar$|green$|red$",
                                                                       stimLog[grep("@log",stimLog):length(stimLog)])]

stimdf <- data.frame(matrix(stimLogCleaning1,ncol=2,byrow = T))
stimdf$X2<-as.numeric(sub(" [[:graph:]]* [[:graph:]]*$","",stimdf$X2))
names(stimdf) <- c("shader","seconds")
stimdiff <- diff(as.integer(stimdf$shader))
stimTransitions <- stimdf[grep(-2,stimdiff),]

stimTransitionsFull<-rbind(stimdf[grep(1,stimdiff),][2,],stimTransitions) # this records all of the Stim frames which occur just before the transition to a new shader. I will always want to take the time of the next frame
newRows <- c()
for (i in 1:nrow(stimTransitionsFull)){
  newRows <- c(newRows,as.integer(row.names(stimTransitionsFull[i,]))+1)
}
if (!all(stimdf[newRows,"shader"]=="bar")){
  print("Not all of the start times are Bar shaders :(")
}
