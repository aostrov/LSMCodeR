completedMats <- dir(file.path(LSMCodeRConfig$srcdir,"objects"),patt=".mat",full=T)

completedMats.list <- list()
for (mats in completedMats) {
  # print(substr(basename(mats),1,4))
  # print(names(readRDS(mats)))
  completedMats.list[[substr(basename(mats),1,4)]] <- readRDS(mats)
}

fishGenos <- read.csv(file.path(LSMCodeRConfig$srcdir,"models","DSLM-fish.csv"))
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
      analysisDF.stimulus <- rbind(analysisDF.stimulus,analysisDF.zplane)
    }
    analysisDF.stimulus$stimulus <- as.factor(stimulus)
    analysisDF.animals <- rbind(analysisDF.animals,analysisDF.stimulus)
  }
  analysisDF.animals$animal <- as.factor(animal)
  analysisDF <- rbind(analysisDF,analysisDF.animals)
}


for (fish in 1:length(fishFullGeno)) {
  analysisDF$geno[substr(analysisDF$animal,1,3)==names(fishFullGeno[fish])] <- fishFullGeno[fish]
}
analysisDF$geno <- as.factor(analysisDF$geno)

# analysisDF.block <- analysisDF[substr(analysisDF$stimulus,1,1)==1,]
# analysisDF.stimulus <- analysisDF[substr(analysisDF$stimulus,3,3)==3,]
ggplot(subset(analysisDF,
              background.mean>50 &
                background.mean<400 & 
                dff.max > 0.5 &
                snr.mean > 3),
       aes(background.mean,dff.max)) + 
  geom_jitter(aes(color=animal,size=snr.mean),alpha=0.3) + 
  facet_wrap(~geno, ncol=3) + 
  ylim(c(-1,8)) + xlim(c(-10,320))

ggplot(subset(analysisDF,
              background.mean>20),
       aes(geno,snr.max)) + 
  geom_boxplot(notch = T)