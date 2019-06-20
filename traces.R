# require('signal')
# 
# aaaa <- readRDS(file.path(LSMCodeRConfig$srcdir,"toDiscardEventually/olderThings/AAAA-20190211-SabineSimple-laser_3-SP.mat.1.1.RDS"))
# z="z_6"

animalAverageDF <- data.frame()
animals <- dir(file.path(LSMCodeRConfig$srcdir,"toDiscardEventually"),patt="1.RDS$",full=T)
for (x in animals) {
  animal <- readRDS(x)
  genotype <- as.character(unique(fishGenos[fishGenos$Fish==substr(basename(x),1,3),]$Full_geno))
  
  avg <- data.frame()
  stats <- data.frame(z=character(),roi=character(),back.mean=numeric(),back.sd=numeric())
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
      if (length(names(animal))==1) {
        deltaff <- deltaff[seq(1,length(deltaff),by = 20)]
      }
      if( any(deltaff[5:15]>0.1) ) {
        # print(roi)
        stimBlock=substr(basename(x),
                         regexpr(pattern = ".mat.",basename(x))+5,
                         regexpr(".RDS",basename(x))-1)
        avg <- rbind(avg,deltaff)
        stats <- rbind(stats,data.frame(z=z,roi=roi,trialName=substr(basename(x),1,4),
                                        stimBlock=stimBlock,
                                        back.mean=back.mean,
                                        back.sd=back.sd))
        
        # count2=count2+1
      }
    }
    # count=count+1
    ncol.filtered <- length(filtered)
  }
  if (nrow(avg)==0 | nrow(stats)==0) {
    print(paste("skipping",basename(x),"because it has 0 rows"))
  } else {
    df=cbind(avg,stats)
    cutoff <- c()
    colnames(df) <- c(paste("X",1:length(df[,grepl("X",colnames(df))]),sep="."),
                      "z","roi",
                      "animalTrial","stimBlock",
                      "back.mean","back.sd")
    length.data=length(paste("X",1:length(df[,grepl("X",colnames(df))])))
    # df2=df[apply(X = df[,grepl("X",colnames(df))],MARGIN = 1,FUN = function(x) max(x)>0.5),]
    # 
    # df3=df[apply(X = df,MARGIN = 1,FUN = function(x) {
    #   max(x[1:length.data]) > (x[34] + (3 * x[35]))
    # }),]
    
    df3=df[apply(X = df,MARGIN = 1,FUN = function(x) {
      max(as.numeric(x[1:31])) > ( as.numeric(x[length(x)-1]) + (7 * as.numeric(x[length(x)])) )
    }),]
    
    if (nrow(df3)>0) {
      # get the average for the animal
      trialAverage <- as.numeric(colMeans(df3[,grepl("X",colnames(df3))]))
      trialAverageDF <- data.frame(trailAverage=trialAverage,animalTrial=unique(df3$animalTrial),stimBlock=unique(df3$stimBlock))
      trialAverageDF$time <- as.numeric(row.names(trialAverageDF))/5 # let's assume this is true
      trialAverageDF$genotype <- genotype
      animalAverageDF <- rbind(animalAverageDF,trialAverageDF)
      
      # plot all ROIs and all Z planes
      maxDFF <- max(apply(X = df3[,grepl("X",colnames(df3))],MARGIN = 2, max))
      maxDFF <- maxDFF + (maxDFF * 0.1)
      pdf(file = file.path(LSMCodeRConfig$srcdir,"toDiscardEventually",paste(basename(x),".pdf",sep="")))
      if (length(names(animal))==1) {
        plot(x=seq(0.2,6.2,0.2),y=animal[[1]][[1]]$signal[seq(1,length(animal[[1]][[1]]$signal),by = 20)],
             ylim = c(0,maxDFF),type = "n", axes = T,frame.plot = T,xlab = "Response window (seconds)",ylab = "dF/F")
      } else {
        plot(x=seq(0.2,6.2,0.2),y=animal[[1]][[1]]$signal,
             ylim = c(0,maxDFF),type = "n", axes = T,frame.plot = T,xlab = "Response window (seconds)",ylab = "dF/F")
        
      }
      
      zByColor <- data.frame(z=unique(df$z),color=viridis(length(unique(df$z))))
      for (i in seq(1,nrow(df3))) {
        lines(x=seq(0.2,6.2,0.2),y=as.numeric(df3[i,grepl("X",colnames(df))]),col=zByColor[df3$z[i],"color"])
      }
      lines(x=seq(0.2,6.2,0.2),y=as.numeric(colMeans(df3[,grepl("X",colnames(df3))])),col="magenta",lwd=4)
      dev.off()
      
    } else {
      print(paste("Skipping",basename(x),"because it does not have a response above 7sd"))
    }
    
    
    
  }
}
saveRDS(animalAverageDF,file=file.path(LSMCodeRConfig$srcdir,"objects","animalAverageTraces.RDS"))
