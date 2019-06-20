trial2Animal <- function(trial){
  return(unique(substr(as.character(trial),1,3)))
}

sampleDFRownamesByAnimal <- function(numberOfSamples,animal,df){
  rows <- rownames(df[grep(paste("^",animal,sep = ""),df$animal),])
  if (numberOfSamples > length(rows)) {
    warning(paste(numberOfSamples," is greater than the number of rows (",length(rows),"). Setting the number of rows to 1/3 of the total rows",sep=""))
    numberOfSamples <- ceiling(length(rows)*0.33)
    
  }
  return(sample(rows,numberOfSamples,replace = F))
}

ggplot(data=tempDF.animals[is.finite(tempDF.animals$dff.max),],aes(background.mean,dff.max)) + 
  geom_point() + 
  geom_point(data=tempDF.animals[is.finite(tempDF.animals$dff.max) & tempDF.animals$background.mean > 11,],aes(color="red")) +
  xlim(c(0,80))


animals2=unique(substr(unique(analysisDF$animal),1,3))
tempDF.animals95 <- data.frame()
for (animal in animals2){
  tempDF.animals <- analysisDF[grepl(paste("^",animal,sep=""),analysisDF$animal),]
  dff95=quantile(tempDF.animals$dff.max[is.finite(tempDF.animals$dff.max) & 
                                          tempDF.animals$dff.max!=0],
                 probs = seq(0,1,0.01))["95%"]
  back95=quantile(tempDF.animals$background.mean,probs = seq(0,1,0.01))["90%"]
  xxx=tempDF.animals[tempDF.animals$background.mean >= back95,]
  yyy=tempDF.animals[tempDF.animals$dff.max >= dff95 & is.finite(tempDF.animals$dff.max),]
  tempDF.animals <- tempDF.animals[intersect(rownames(xxx),rownames(yyy)),]
  if (nrow(tempDF.animals)>0){
    tempDF.animals95 <- rbind(tempDF.animals95,tempDF.animals)
  } else {
    print(paste("Skipping", animal,"as it has no example points to use."))
  }
  
}




animals <- trial2Animal(tempDF.animals95$animal)
examples <- c()
for (animal in animals) {
  print(animal)
  rowname.sample <- sampleDFRownamesByAnimal(2,animal,tempDF.animals95)
  print(rowname.sample)
  for (row in 1:length(rowname.sample)){
    getFramesFromStimParamListFromAnalysisDF(rowname.sample[row],
                                             file=file.path("C:/Users/Aaron/Desktop/exampleTraces/",
                                                            paste(animal,
                                                                  "_",
                                                                  rowname.sample[row],
                                                                  "-90percentBackground.nrrd",sep="")))
    
  }
  
  examples <- c(examples,rowname.sample)
}
names(examples) <- animals



ggplot(data=tempDF.animals95,aes(background.mean,dff.max)) + geom_point(aes(color=geno)) + xlim(c(0,100))
