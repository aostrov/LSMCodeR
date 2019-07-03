# require('signal')
# 
# aaaa <- readRDS(file.path(LSMCodeRConfig$srcdir,"toDiscardEventually/olderThings/AAAA-20190211-SabineSimple-laser_3-SP.mat.1.1.RDS"))
# z="z_6"

animalAverageDF <- data.frame()
animalAverageDF.raw <- data.frame()
notAbove_dFF_threshold <- c()
notAbove_sdThreshold <- c()
sdThreshold = 7
dffThreshold = 0.1

fileDir <- "/Volumes/Prospero/traceFiles/"

# animals <- dir(file.path(LSMCodeRConfig$srcdir,"toDiscardEventually"),patt="1.RDS$",full=T)
animals <- dir(fileDir,full=T)
animals.clean <- animals[!substr(basename(animals),1,3)%in%fishGenos$Fish[fishGenos$Use=="No"]]
for (x in animals.clean) {
  animal <- readRDS(x)
  genotype <- as.character(unique(fishGenos[fishGenos$Fish==substr(basename(x),1,3),]$Full_geno))
  
  avg <- data.frame()
  avg.raw <- data.frame()
  stats <- data.frame(z=character(),roi=character(),back.mean=numeric(),back.sd=numeric())
  stats.raw <- data.frame(z=character(), roi=character(), raw.dff=numeric() )

  for (z in names(animal)){
    # print(z)
    for (roi in names(animal[[z]])) {
      filtered <- sgolayfilt(animal[[z]][[roi]]$signal,p=2)
      back <- mean(animal[[z]][[roi]]$background)
      back.delta <- ( (animal[[z]][[roi]]$background - back) / back )
      back.mean <- mean(back.delta)
      back.sd <- sd(back.delta)
      delta <- filtered-back
      deltaff <- delta/back
      rawdff <- ( (animal[[z]][[roi]]$signal - back) / back )
      if (length(names(animal))==1) {
        deltaff <- deltaff[seq(1,length(deltaff),by = 20)]
        rawdff <- rawdff[seq(1,length(rawdff),by = 20)]
      }
      if( any(deltaff[5:15]>dffThreshold) ) {
        # print(roi)
        stimBlock=substr(basename(x),
                         regexpr(pattern = ".mat.",basename(x))+5,
                         regexpr(".RDS",basename(x))-1)
        avg <- rbind(avg,deltaff)
        avg.raw <- rbind(avg.raw,rawdff)
        stats <- rbind(stats,data.frame(z=z,roi=roi,trialName=substr(basename(x),1,4),
                                        stimBlock=stimBlock,
                                        back.mean=back.mean,
                                        back.sd=back.sd))
        stats.raw <- rbind(stats.raw, data.frame(z=z,roi=roi,trialName=substr(basename(x),1,4),
                                                 stimBlock=stimBlock,
                                                 raw.dff=rawdff) )
        
        # count2=count2+1
      }
    }
    # count=count+1
    ncol.filtered <- length(filtered)
  }
  if (nrow(avg)==0 | nrow(stats)==0) {
    print(paste("skipping",basename(x),"because it has 0 rows with a df/f above",dffThreshold))
    notAbove_dFF_threshold <- c(notAbove_dFF_threshold,x)
  } else {
    df=cbind(avg,stats)
    df.raw=cbind(avg.raw,stats.raw)
    cutoff <- c()
    colnames(df) <- c(paste("X",1:length(df[,grepl("X",colnames(df))]),sep="."),
                      "z","roi",
                      "animalTrial","stimBlock",
                      "back.mean","back.sd")
    
    colnames(df.raw) <- c(paste("X",1:length(df.raw[,grepl("X",colnames(df.raw))]),sep="."),
                      "z","roi",
                      "animalTrial","stimBlock",
                      "raw.dff")
    
    length.data=length(paste("X",1:length(df[,grepl("X",colnames(df))])))

    df3=df[apply(X = df,MARGIN = 1,FUN = function(x) {
      max(as.numeric(x[1:31])) > ( as.numeric(x[length(x)-1]) + (sdThreshold * as.numeric(x[length(x)])) )
    }),]
    df3.raw=df.raw[apply(X = df,MARGIN = 1,FUN = function(x) {
      max(as.numeric(x[1:31])) > ( as.numeric(x[length(x)-1]) + (sdThreshold * as.numeric(x[length(x)])) )
    }),]
    
    ############################################################################
    ## Write out a new file that has the traces being slightly better aligned ##
    ############################################################################
    testData <- subset(df3.raw,X.7>0.3)
    if (nrow(testData)==0) {
      print(paste("No rows for X.7>0.3 for",x))
    } else {
      write.nrrd(x = as.matrix(testData[nrow(testData)/2, grepl(x = colnames(testData),pattern = "^X")]) , file = "~/Desktop/tmp/referenceTrace_61.nrrd")
      # unlink(fixed)
      
      registeredDF <- data.frame()
      for (i in 1:nrow(testData)) { #1:nrow(df3)
        write.nrrd(x = as.matrix(testData[i, grepl(x = colnames(testData),
                                                   pattern = "^X")
                                          ]
                                 ), 
                          file = "~/Desktop/tmp/moving.nrrd"
                   )
        print(paste("working on row",i))
        system("sh ~/Desktop/2dReg-test.sh",ignore.stdout = TRUE)
        registered=c(read.nrrd(file="~/Desktop/tmp/reformatted.nrrd"))
        registeredDF <- rbind(registeredDF,registered)
      }
      saveRDS(registeredDF,
              file = file.path(fileDir,
                               paste(basename(x),"-aligned.RDS",sep="")))
      
    }
     
    pdf(file = file.path(fileDir,
                         paste(basename(x),".pdf",sep="")))
      plot(x=seq(0,6,0.2),
           y = registeredDF[1,],
           type = "l",ylim=c(-0.05,max(registeredDF[,6:18])))
      
      for ( j in seq(2,nrow(registeredDF),10) ) {
        lines(x=seq(0,6,0.2),y=registeredDF[j,])
      }
      lines(x=seq(0,6,0.2),y=colMeans(registeredDF),col="green",lwd=2)
      lines(x=seq(0,6,0.2),y=colMeans(testData[,grepl("X",colnames(testData))]),col="magenta",lwd=2)
    dev.off()
    
    ############################################################################
    ############################################################################
    
  #   if (nrow(df3)>0) {
  #     # get the average for the animal
  #     trialAverage <- as.numeric(colMeans(df3[,grepl("X",colnames(df3))]))
  #     trialAverageDF <- data.frame(trailAverage=trialAverage,animalTrial=unique(df3$animalTrial),stimBlock=unique(df3$stimBlock))
  #     trialAverageDF$time <- (as.numeric(row.names(trialAverageDF))/5 - 0.2) # let's assume this is true
  #     trialAverageDF$genotype <- genotype
  #     animalAverageDF <- rbind(animalAverageDF,trialAverageDF)
  #     
  #     # plot all ROIs and all Z planes
  #     maxDFF <- max(apply(X = df3[,grepl("X",colnames(df3))],MARGIN = 2, max))
  #     maxDFF <- maxDFF + (maxDFF * 0.1)
  #     # pdf(file = file.path(LSMCodeRConfig$srcdir,"toDiscardEventually",paste(basename(x),".pdf",sep="")))
  #     #   if (length(names(animal))==1) {
  #     #     plot(x=seq(0,6,0.2),y=animal[[1]][[1]]$signal[seq(1,length(animal[[1]][[1]]$signal),by = 20)],
  #     #          ylim = c(0,maxDFF),type = "n", axes = T,frame.plot = T,xlab = "Response window (seconds)",ylab = "dF/F")
  #     #   } else {
  #     #     plot(x=seq(0,6,0.2),y=animal[[1]][[1]]$signal,
  #     #          ylim = c(0,maxDFF),type = "n", axes = T,frame.plot = T,xlab = "Response window (seconds)",ylab = "dF/F")
  #     #     
  #     #   }
  #     #   
  #     #   zByColor <- data.frame(z=unique(df$z),color=viridis(length(unique(df$z))))
  #     #   for (i in seq(1,nrow(df3))) {
  #     #     lines(x=seq(0,6,0.2),y=as.numeric(df3[i,grepl("X",colnames(df))]),col=zByColor[df3$z[i],"color"])
  #     #   }
  #     #   lines(x=seq(0,6,0.2),y=as.numeric(colMeans(df3[,grepl("X",colnames(df3))])),col="magenta",lwd=4)
  #     # dev.off()
  #     
  #     # raw
  #     trialAverage.raw <- as.numeric(colMeans(df3.raw[,grepl("X",colnames(df3.raw))]))
  #     trialAverageDF.raw <- data.frame(trailAverage=trialAverage.raw,
  #                                      animalTrial=unique(df3.raw$animalTrial),
  #                                      stimBlock=unique(df3.raw$stimBlock))
  #     trialAverageDF.raw$time <- (as.numeric(row.names(trialAverageDF.raw))/5 - 0.2) # let's assume this is true
  #     trialAverageDF.raw$genotype <- genotype
  #     animalAverageDF.raw <- rbind(animalAverageDF.raw,trialAverageDF.raw)
  #     
  #     # plot all ROIs and all Z planes
  #     maxDFF <- max(apply(X = df3.raw[,grepl("X",colnames(df3.raw))],MARGIN = 2, max))
  #     maxDFF <- maxDFF + (maxDFF * 0.1)
  #     pdf(file = file.path(LSMCodeRConfig$srcdir,"toDiscardEventually",paste("raw-",basename(x),".pdf",sep="")))
  #       if (length(names(animal))==1) {
  #         plot(x=seq(0,6,0.2),y=animal[[1]][[1]]$signal[seq(1,length(animal[[1]][[1]]$signal),by = 20)],
  #              ylim = c(0,maxDFF),type = "n", axes = T,frame.plot = T,xlab = "Response window (seconds)",ylab = "dF/F")
  #       } else {
  #         plot(x=seq(0,6,0.2),y=animal[[1]][[1]]$signal,
  #              ylim = c(0,maxDFF),type = "n", axes = T,frame.plot = T,xlab = "Response window (seconds)",ylab = "dF/F")
  #         
  #       }
  #       
  #       zByColor <- data.frame(z=unique(df.raw$z),color=viridis(length(unique(df.raw$z))))
  #       for (i in seq(1,nrow(df3.raw))) {
  #         lines(x=seq(0,6,0.2),y=as.numeric(df3.raw[i,grepl("X",colnames(df))]),col=zByColor[df3.raw$z[i],"color"])
  #       }
  #       lines(x=seq(0,6,0.2),y=as.numeric(colMeans(df3.raw[,grepl("X",colnames(df3.raw))])),col="magenta",lwd=4)
  #     dev.off()
  #     
  #   } else {
  #     print(paste("Skipping",basename(x),"because it does not have a response above",sdThreshold,"sd"))
  #     notAbove_sdThreshold <- c(notAbove_sdThreshold,x)
  #   }
  #   
  #   
  #   
  }
}

attr(notAbove_sdThreshold,"sdThreshold") <- sdThreshold
attr(notAbove_dFF_threshold,"dFF_Threshold") <- dffThreshold
saveRDS(animalAverageDF,file=file.path(LSMCodeRConfig$srcdir,"objects","animalAverageTraces.RDS"))

newDF.big <- data.frame()
for (j in unique(substring(animalAverageDF.raw$animalTrial,1,3))){
  byAnimal=animalAverageDF.raw[grepl(paste("^",j,sep=""),animalAverageDF.raw$animalTrial),]
  newDF <- data.frame()
  for (i in seq(0,6,0.2)){
    newDF <- rbind(newDF,data.frame(average=mean(byAnimal[as.character(byAnimal$time)==as.character(i),"trailAverage"]),time=i))
    
  }
  newDF$genotype <- unique(byAnimal$genotype)
  newDF$animal <- j
  newDF.big <- rbind(newDF.big,newDF)
}
saveRDS(newDF.big,file=file.path(LSMCodeRConfig$srcdir,"objects","animalAverageTracesByAnimalNotROI.RDS"))
