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
for (k in 1:length(completedMats.list)){
  analysisDF.animals <- c()
  animal <- names(completedMats.list[k])
  cat("",sep="\n")
  print(paste("animal:",animal))
  for (j in 1: length(completedMats.list[[k]][[1]])){
    analysisDF.stimulus <- c()
    stimulus <- substr(names(completedMats.list[[k]][[1]][j]),
                       nchar(names(completedMats.list[[k]][[1]][j]))-2,
                       nchar(names(completedMats.list[[k]][[1]][j])))
    # print(paste("stimulus:",stimulus))
    for (i in 1:length(completedMats.list[[k]][[1]][[j]])) {
      analysisDF.zplane <- c()
      z_plane <- names(completedMats.list[[k]][[1]][[j]][i])
      # print(paste("z_plane:",z_plane))
      cat(".")
      analysisDF.zplane <- 
        completedMats.list[[k]][[1]][[j]][[i]][,c("snr.mean","background.mean","dff.mean","snr.max","dff.max")]
      analysisDF.zplane$z_plane <- as.factor(z_plane)
      analysisDF.zplane$roi <- row.names(analysisDF.zplane)
      analysisDF.stimulus <- rbind(analysisDF.stimulus,analysisDF.zplane)
    }
    analysisDF.stimulus$stimulus <- as.factor(stimulus)
    analysisDF.animals <- rbind(analysisDF.animals,analysisDF.stimulus)
  }
  analysisDF.animals$animal <- as.factor(animal)
  analysisDF.animals$date <- fishGenos[k,"Date"]
  analysisDF.animals$laser <- fishGenos[k,"Laser"]
  analysisDF.animals$use <- fishGenos[k,"Use"]
  analysisDF <- rbind(analysisDF,analysisDF.animals)
  
}


for (fish in 1:length(fishFullGeno)) {
  analysisDF$geno[substr(analysisDF$animal,1,3)==names(fishFullGeno[fish])] <- fishFullGeno[fish]
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

# Set things up to get average data by animal
animalSummaryDF2 <- data.frame(dff.mean=double(),
                               dff.75quantile=double(),
                               dff.2sd=double(),
                               snr.mean=double(),
                               background.mean=double(),
                               #animal=character(),
                               geno=character(),
                               animalKey=character(),
                               bestStim=character())

countingDF <- data.frame(animal=character(),
                         geno=character())

animals2=unique(substr(unique(analysisDF$animal),1,3))

for (animal in animals2){
  tempDF.animals <- analysisDF[grepl(paste("^",animal,sep=""),analysisDF$animal),]
  tempDF.animals.subset <- subset(tempDF.animals,
                                  background.mean>50 &
                                    background.mean<400)
  # stimuls.block
  # get block 1
  block1 <- tempDF.animals.subset[grepl("[1-4].1",tempDF.animals.subset$stimulus),]
  block1.sub <- subset(block1,
                       background.mean>10 &
                         background.mean<400)
  stims <- unique(block1.sub$stimulus)
  response=0
  bestResponse=character()
  for (stim in stims){
    x=block1.sub[block1.sub$stimulus==stim,"dff.max"]
    if (mean(x)>response) {
      response = mean(x)
      bestResponse = stim
    }
  }
  
  tempDF.animals.subset.bestStim <- subset(tempDF.animals.subset, stimulus==bestResponse)
  
  if (nrow(tempDF.animals.subset.bestStim)>0) {
    tempSummary=data.frame(dff.mean=mean(tempDF.animals.subset.bestStim$dff.max[is.finite(tempDF.animals.subset.bestStim$dff.max)]),
                           dff.75quantile=quantile(tempDF.animals.subset$dff.max[is.finite(tempDF.animals.subset.bestStim$dff.max)])[["75%"]],
                           dff.2sd=get2SD(4),
                           snr.mean=mean(tempDF.animals.subset.bestStim$snr.max[is.finite(tempDF.animals.subset.bestStim$snr.max)]),
                           background.mean=mean(tempDF.animals.subset.bestStim$background.mean[is.finite(tempDF.animals.subset.bestStim$background.mean)]),
                           #animal=unique(tempDF.animals.subset.bestStim$animal),
                           geno=unique(tempDF.animals.subset.bestStim$geno),
                           animalKey=animal,
                           bestStim=bestResponse)
    
    animalSummaryDF2 <- rbind(animalSummaryDF2,tempSummary)
    countingDF <- rbind(countingDF,data.frame(animal=animal,geno=unique(animalSummaryDF2[animalSummaryDF2$animalKey==animal,"geno"])))
  }
  
}

get2SD<- function(sd.level=2){
  dff.mean <- mean(tempDF.animals.subset.bestStim$dff.max[is.finite(tempDF.animals.subset.bestStim$dff.max)])
  dff.sd <- sd(tempDF.animals.subset.bestStim$dff.max[is.finite(tempDF.animals.subset.bestStim$dff.max)])
  if (is.na(dff.sd)) dff.sd <- 0
  return(dff.mean + sd.level*dff.sd)
}

ggplot(subset(animalSummaryDF2, animalSummaryDF2$animal!="AAM" ),
       aes(geno,dff.75quantile)) + 
  geom_boxplot(notch = F) +
  xlab("Genotype") + ylab("DF/F (75th quantile)") + theme(legend.position = "none") + geom_jitter(aes(color=animalKey))

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
                                animalKey=character()
)

animals2=unique(substr(unique(analysisDF$animal),1,3))

for (animal in animals2){
  tempDF.animals <- analysisDF[grepl(paste("^",animal,sep=""),analysisDF$animal),]
  tempDF.animals.subset <- subset(tempDF.animals,
                                  background.mean>20 &
                                    dff.max!=0)
  if (nrow(tempDF.animals.subset)>0) {
    tempSummary=data.frame(dff.mean=mean(tempDF.animals.subset$dff.max[is.finite(tempDF.animals.subset$dff.max)]),
                           background.mean=mean(tempDF.animals.subset$background.mean[is.finite(tempDF.animals.subset$background.mean)]),
                           snr.mean=mean(tempDF.animals.subset$snr.max[is.finite(tempDF.animals.subset$snr.max)]),
                           geno=unique(tempDF.animals.subset$geno),
                           animalKey=animal)
    
    summaryDF.noZeros <- rbind(summaryDF.noZeros,tempSummary)
  }
  
}
ggplot(subset(summaryDF.noZeros, summaryDF.noZeros$animalKey!="AAM" ),
       aes(geno,snr.mean)) + 
  geom_boxplot(notch = T) + 
  geom_point(aes(color=animalKey)) + 
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
ggplot(subset(analysisDF.subset,animal=="AAIA" & z_plane%in%unique(zz$z) & roi%in%zz$roi),aes(background.mean,dff.max)) + 
  geom_jitter(aes(color=z_plane)) + facet_wrap(~stimulus, ncol=4)

