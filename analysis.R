objDir <- "C:/Users/Aaron/Documents/R/LSMCodeR/objects/"
matList <- readRDS(file.path(objDir,"matList1.RDS"))
stimParList <- readRDS(file.path(objDir,'stimParList.RDS'))
# myFile <- readRDS("AAAB-20190211-SabineSimple-laser_3-SP.mat.1.1.RDS")
fishGenos <- read.csv(file.path(LSMCodeRConfig$srcdir,"models","DSLM-fish.csv"))
fishFullGeno <- c()
for (i in unique(substr(names(matList),1,3))){
  fishFullGeno <- c(fishFullGeno,as.character(unique(fishGenos[fishGenos$Fish==i,"Full_geno"])))
}
names(fishFullGeno) <- unique(substr(names(matList),1,3))



analysisDF <- c()
for (k in 1:length(matList)){
  analysisDF.animals <- c()
  animal <- names(matList[k])
  print(paste("animal:",animal))
  for (j in 1: length(matList[[k]])){
    analysisDF.stimulus <- c()
    stimulus <- substr(names(matList[[k]][j]),nchar(names(matList[[k]][j]))-2,nchar(names(matList[[k]][j])))
    print(paste("stimulus:",stimulus))
    for (i in 1:length(matList[[k]][[j]])) {
      analysisDF.zplane <- c()
      z_plane <- names(matList[[k]][[j]][i])
      analysisDF.zplane <- matList[[k]][[j]][[i]][,c("snr.mean","background.mean","dff.mean")]
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
              background.mean>20 &
                background.mean<400 &
                snr.mean > 3),
       aes(background.mean,dff.mean)) + 
  geom_jitter(aes(color=animal,size=snr.mean),alpha=0.3) + 
  facet_wrap(~geno, ncol=4)

ggplot(subset(analysisDF,
              background.mean>20),
       aes(geno,snr.mean)) + 
  geom_boxplot(notch = T)