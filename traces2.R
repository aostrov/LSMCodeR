#################
### Functions ###
#################

# this function commits the horrible crime of 
# relying on a dataframe in the global environment
# and is therefore single purpose and not very useful
makeAverageTrace <- function(plotTraces=F){
  df.sample.all <- data.frame()
  for (stim in unique(df.raw$stimBlock)) {
    # print(paste("sampling stim:",stim))
    df.sampled <- sample.dataframe(
      df.raw[df.raw$stimBlock==stim,],
                                   ifelse(sum(df.raw$stimBlock==stim)<200,
                                          sum(df.raw$stimBlock==stim)*0.75,
                                          200)
      )[,1:31]
    df.sample.all <- rbind(df.sample.all,df.sampled)
  }
  
  df.final <- data.frame()
  
  # Apply some thresholding
  for (i in 1:nrow(df.sample.all)) {
    if (sum(df.sample.all[i,5:20]>0.25)>3){
      df.final <- rbind(df.final,df.sample.all[i,])
      cat("+")
    } else {
      cat("-")
    }
  }
  
  if (nrow(df.final)==0) return(data.frame())
  
  avg <- as.numeric(colMeans(df.final))
  df.shifted <- data.frame()
  for (row.shifted in 1:nrow(df.final)) {
    shift <- which.max(sgolayfilt(avg,p=2)) - which.max(sgolayfilt(as.numeric(df.final[row.shifted,],p=2)))
    df.shifted <- rbind(df.shifted,shift1Darray( as.numeric(df.final[row.shifted,]) , ifelse(abs(shift)>5,0,shift) ))
  }
  
  avg2 <- as.numeric(colMeans(df.shifted))
  
  if (plotTraces){
    plot(as.numeric(df.shifted[1,]),type = "n",ylim = c(0,(max(df.shifted)+0.1)))
    
    for (i in 1:nrow(df.shifted)) {
      lines(as.numeric(df.shifted[i,]),col=ifelse(sum(df.shifted[i,5:20]>0.25)>3,"green","red"))
      # if (sum(df.sampled[i,5:20]>0.25)>3){
      #   df.final <- rbind(df.final,df.sampled[i,])
      # }
    }
    
  }
  
  xxx=data.frame(averageValue=avg2)
  xxx$time <- seq(0,(0.2*nrow(xxx)-0.2),0.2)
  return(xxx)
}

# shiftBy > 0 shifts the plot to the right,
# shiftBy < 0 shifts the plot to the left...
# 0s are used to keep the plot the same length as the 
# input, not sure if an average would be better
shift1Darray <- function(array2Shift,shiftBy){
  if (shiftBy==0) return(array2Shift)
  if (shiftBy > 0) {
    newArray <- c(rep(0,abs(shiftBy)),
                  array2Shift[-((length(array2Shift)-abs(shiftBy)+1):length(array2Shift))])
    return(newArray)
  }
  if (shiftBy < 0) {
    newArray <- c(array2Shift[-(1:(0+abs(shiftBy)))],
                  rep(0,abs(shiftBy)))
    return(newArray)
  }
}

#################
#################

# A bit of setup
dffThreshold = 0.1
notAbove_dFF_threshold <- c()
notAbove_sdThreshold <- c()
sdThreshold = 7
# stats.raw <- data.frame(z=character(), roi=character(), raw.dff=numeric() )
# avg.raw <- data.frame()
fileDir <- "/Volumes/Prospero/traceFiles/"


animals <- dir(fileDir,full=T,patt="1.RDS$")
animals.clean <- animals[
  !substr(basename(animals),1,3)%in%fishGenos$Fish[fishGenos$Use=="No"] #&
    #!substr(basename(animals),1,3)%in%fishGenos$Fish[fishGenos$Full_geno=="iGABASnFr"]
  ]
animalBaseName <- unique(substr(basename(animals.clean),1,3))

df.averages <- data.frame()

for (uniqueAnimal in animalBaseName[]){#27:length(animalBaseName)
  stats.raw <- data.frame(z=character(), roi=character(), raw.dff=numeric() )
  avg.raw <- data.frame()
  
  for (x in animals.clean[grep(paste("^",uniqueAnimal,sep=""),basename(animals.clean))]) {
    print(x)
    animal <- readRDS(x)
    genotype <- as.character(unique(fishGenos[fishGenos$Fish==substr(basename(x),1,3),]$Full_geno))
    
    for (z in names(animal)){
      print(".")
      for (roi in names(animal[[z]])) {
        back <- mean(animal[[z]][[roi]]$background)
        rawdff <- ( (animal[[z]][[roi]]$signal - back) / back )
        if (length(names(animal))==1) {
          rawdff <- rawdff[seq(1,length(rawdff),by = 20)]
        }
        if( any(rawdff[5:15]>dffThreshold) ) {
          stimBlock=substr(basename(x),
                           regexpr(pattern = ".mat.",basename(x))+5,
                           regexpr(".RDS",basename(x))-1)
          avg.raw <- rbind(avg.raw,rawdff)
          
          stats.raw <- rbind(stats.raw, data.frame(z=z,roi=roi,trialName=substr(basename(x),1,4),
                                                   stimBlock=stimBlock,
                                                   raw.dff=rawdff) )
          
        }
      }
    }
    if (nrow(avg.raw)==0 | nrow(stats.raw)==0) {
      print(paste("skipping",basename(x),"because it has 0 rows with a df/f above",dffThreshold))
      notAbove_dFF_threshold <- c(notAbove_dFF_threshold,x)
    } else {
      df.raw=cbind(avg.raw,stats.raw)
      colnames(df.raw) <- c(paste("X",1:length(df.raw[,grepl("X",colnames(df.raw))]),sep="."),
                            "z","roi",
                            "animalTrial","stimBlock",
                            "raw.dff")
    }
    
  }
  
  
  at <- makeAverageTrace(plotTraces = F)
  if (nrow(at)==0){
    notAbove_dFF_threshold <- c(notAbove_dFF_threshold,uniqueAnimal)
  } else {
    at$animal <- uniqueAnimal
    at$genotype <- genotype
    df.averages <- rbind(df.averages,at)
    
  }
}


