## Functions ##

get2SD<- function(sd.level=2){
  dff.mean <- mean(tempDF.animals.subset.bestStim$dff.max[
    is.finite(tempDF.animals.subset.bestStim$dff.max)
    ])
  dff.sd <- sd(tempDF.animals.subset.bestStim$dff.max[
    is.finite(tempDF.animals.subset.bestStim$dff.max)
    ])
  if (is.na(dff.sd)) dff.sd <- 0
  return(dff.mean + sd.level*dff.sd)
}

## Setup ##

completedMats <- dir(file.path(LSMCodeRConfig$srcdir,"objects"),patt=".mat",full=T)

completedMats.list <- list()
for (mats in completedMats) {
  # print(substr(basename(mats),1,4))
  # print(names(readRDS(mats)))
  completedMats.list[[substr(basename(mats),1,4)]] <- readRDS(mats)
}

fishFullGeno <- c()
for (i in unique(substr(names(completedMats.list),1,3))){
  fishFullGeno <- c(fishFullGeno,as.character(unique(fishGenos[fishGenos$Fish==i,"Full_geno"])))
}
names(fishFullGeno) <- unique(substr(names(completedMats.list),1,3))



analysisDF <- c()
for (trial in 1:length(completedMats.list)){
  analysisDF.animals <- c()
  animal <- names(completedMats.list[trial])
  cat("",sep="\n")
  print(paste("animal:",animal))
  for (stimRepetition in 1: length(completedMats.list[[trial]][[1]])){
    analysisDF.stimulus <- c()
    stimulus <- substr(names(completedMats.list[[trial]][[1]][stimRepetition]),
                       nchar(names(completedMats.list[[trial]][[1]][stimRepetition]))-2,
                       nchar(names(completedMats.list[[trial]][[1]][stimRepetition])))
    # print(paste("stimulus:",stimulus))
    for (z in 1:length(completedMats.list[[trial]][[1]][[stimRepetition]])) {
      analysisDF.zplane <- c()
      z_plane <- names(completedMats.list[[trial]][[1]][[stimRepetition]][z])
      # print(paste("z_plane:",z_plane))
      cat(".")
      analysisDF.zplane <- 
        completedMats.list[[trial]][[1]][[stimRepetition]][[z]][,c("snr.mean",
                                                                   "background.mean",
                                                                   "dff.mean",
                                                                   "snr.max",
                                                                   "dff.max")]
      analysisDF.zplane$z_plane <- as.factor(z_plane)
      analysisDF.zplane$roi <- row.names(analysisDF.zplane)
      analysisDF.stimulus <- rbind(analysisDF.stimulus,analysisDF.zplane)
    }
    analysisDF.stimulus$stimulus <- as.factor(stimulus)
    analysisDF.animals <- rbind(analysisDF.animals,analysisDF.stimulus)
  }
  analysisDF.animals$animal <- as.factor(animal)
  # analysisDF.animals$date <- fishGenos[trial,"Date"]
  analysisDF.animals$laser <- fishGenos[trial,"Laser"]
  analysisDF.animals$use <- fishGenos[trial,"Use"]
  analysisDF <- rbind(analysisDF,analysisDF.animals)
  
}


for (fish in names(fishFullGeno)) {
  analysisDF$geno[substr(analysisDF$animal,1,3)==fish] <- fishFullGeno[fish]
  analysisDF$date[grepl(paste("^",fish,sep=""),analysisDF$animal)] <- unique(fishGenos[fishGenos$Fish==fish,"Date"])
}

analysisDF$geno <- as.factor(analysisDF$geno)
saveRDS(analysisDF,file=file.path(LSMCodeRConfig$srcdir,"objects",paste("analysisDF",".RDS",sep="")))

analysisDF.subset <- subset(analysisDF,background.mean > 20)

# analysisDF.block <- analysisDF[substr(analysisDF$stimulus,1,1)==1,]
# analysisDF.stimulus <- analysisDF[substr(analysisDF$stimulus,3,3)==3,]
ggplot(subset(analysisDF,
              background.mean>50 &
                background.mean<400 & 
                dff.max > 0.75 &
                snr.mean > 3),
       aes(background.mean,dff.max)) + 
  geom_jitter(aes(color=animal,size=snr.mean),alpha=0.3) + 
  facet_wrap(~geno, ncol=3) + 
  ylim(c(-1,8)) + xlim(c(-10,320)) + scale_color_viridis_d()

ggplot(subset(analysisDF,
              background.mean>50 &
                background.mean<400 & 
                dff.max > 0.75 &
                snr.mean > 3),
       aes(geno,snr.max)) + 
  geom_boxplot(notch = F)

ggplot(analysisDF) + geom_histogram(aes(geno),stat="count")



ggplot(subset(animalSummaryDF2, animalSummaryDF2$animal!="AAM" ),
       aes(geno,dff.75quantile)) + 
  geom_boxplot(notch = F) +
  xlab("Genotype") + ylab("DF/F (75th quantile)") + 
  theme(legend.position = "none") + geom_jitter(aes(color=animalKey))

ggplot(countingDF) + geom_histogram(aes(geno),stat="count") + scale_y_continuous(breaks = seq(0, 12))

# This is a bit more straightforward, and it provides what I think is a bit more of a clean
# analysis by dropping out all ROIs that have no change in fluorescence after stimulus
# presentation.
#
#
summaryDF.noZeros <- data.frame(dff.mean=double(),
                                background.mean=double(),
                                snr.mean=double(),
                                geno=character(),
                                date=character(),
                                animalKey=character()
)

animals2=unique(substr(unique(analysisDF$animal),1,3))

for (animal in animals2){
  tempDF.animals <- analysisDF[grepl(paste("^",animal,sep=""),analysisDF$animal),]
  tempDF.animals.subset <- subset(tempDF.animals,
                                  background.mean>20 &
                                    dff.max!=0)
  
  
  
  if (nrow(tempDF.animals.subset)>0) {
    print(paste(animal,unique(tempDF.animals.subset$date)))
    tempSummary=data.frame(dff.mean=mean(tempDF.animals.subset$dff.max[is.finite(tempDF.animals.subset$dff.max)]),
                           background.mean=mean(tempDF.animals.subset$background.mean[is.finite(tempDF.animals.subset$background.mean)]),
                           snr.mean=mean(tempDF.animals.subset$snr.max[is.finite(tempDF.animals.subset$snr.max)]),
                           geno=unique(tempDF.animals.subset$geno),
                           date=unique(tempDF.animals.subset$date),
                           animalKey=animal)
    
    summaryDF.noZeros <- rbind(summaryDF.noZeros,tempSummary)
  }
  
}

ggplot(subset(summaryDF.noZeros, summaryDF.noZeros$animalKey!="AAM" & geno!="iGABASnFr"),
       aes(geno,dff.mean)) + 
  geom_boxplot(notch = F) + 
  geom_point(aes(color=animalKey)) + 
  facet_wrap(~date) +
  theme(legend.position = "none")

# Unfortunately, there doesn't seem to be much of a difference
pairwise.wilcox.test(x = summaryDF.noZeros$dff.mean, 
                     g = summaryDF.noZeros$geno, 
                     p.adjust.method = "bonferroni")

# But we can go p-hunting and get a single significant difference
# between 7f EF05 and 6s when looking at SNR
pairwise.wilcox.test(x = summaryDF.noZeros$snr.mean, 
                     g = summaryDF.noZeros$geno, 
                     p.adjust.method = "bonferroni")


# subset the dataset so that only the top responders are plotted
zzz=getArbitraryTopROIsPerZ(aaia,threshold = 0.5)
zz=topRespondersDF(zzz)
ggplot(subset(analysisDF.subset,animal=="AAIA" & z_plane%in%unique(zz$z) & roi%in%zz$roi),
       aes(background.mean,dff.max)) + 
  geom_jitter(aes(color=z_plane)) + facet_wrap(~stimulus, ncol=4)

