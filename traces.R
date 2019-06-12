require('signal')

abra <- readRDS(file.path(LSMCodeRConfig$srcdir,"toDiscardEventually/ABRA-gen_A-laser-1-SabineSimple.mat.1.1.RDS"))
z="z_6"


avg <- data.frame()
stats <- data.frame(z=character(),roi=character())
for (z in names(abra)){
  # print(z)
  for (roi in names(abra[[z]])) {
    filtered <- sgolayfilt(abra[[z]][[roi]]$signal,p=2)
    back <- mean(abra[[z]][[roi]]$background)
    delta <- filtered-back
    deltaff <- delta/back
    if( any(deltaff>0.1) ) {
      # print(roi)
      avg <- rbind(avg,deltaff)
      stats <- rbind(stats,data.frame(z=z,roi=roi))
      
      # count2=count2+1
    }
  }
  # count=count+1
  ncol.filtered <- length(filtered)
}
df=cbind(avg,stats)
colnames(df) <- c(paste("X",1:length(df[,grepl("X",colnames(df))]),sep="."),"z","roi")
df2=df[apply(X = df[,grepl("X",colnames(df))],MARGIN = 1,FUN = function(x) max(x)>0.5),]

# plot all ROIs and all Z planes
maxDFF <- max(apply(X = df2[,grepl("X",colnames(df2))],MARGIN = 2, max))
maxDFF <- maxDFF + (maxDFF * 0.1)
plot(abra[[z]][["ROI_50"]]$signal,
     ylim = c(0,maxDFF),type = "n", axes = T,frame.plot = T,xlab = "Response window",ylab = "dF/F")
zByColor <- data.frame(z=unique(df$z),color=viridis(length(unique(df$z))))
for (i in seq(1,nrow(df2))) {
  lines(as.numeric(df2[i,grepl("X",colnames(df))]),col=zByColor[df2$z[i],"color"])
}
lines(as.numeric(colMeans(df2[,grepl("X",colnames(df2))])),col="magenta",lwd=4)


