require('signal')

abra <- readRDS(file.path(LSMCodeRConfig$srcdir,"toDiscardEventually/ABRA-gen_A-laser-1-SabineSimple.mat.1.1.RDS"))
z="z_6"

selectedROI=79
plot(abra[[z]][["ROI_50"]]$signal,ylim = c(0,1.5),type = "n", axes = T,frame.plot = T,xlab = "Response window",ylab = "dF/F")
count=1
count2=0
avg <- c()
for (z in names(abra)){
  print(z)
  for (roi in names(abra[[z]])) {
    filtered <- sgolayfilt(abra[[z]][[roi]]$signal,p=2)
    back <- mean(abra[[z]][[roi]]$background)
    delta <- filtered-back
    deltaff <- delta/back
    if( any(deltaff>0.5) & ( roi%in%paste("ROI_",c((selectedROI-1):(selectedROI+1)),sep="") ) ) {
      lines(deltaff,col=rainbow(length(names(abra)))[count])
      print(roi)
      avg <- c(avg,deltaff)
      count2=count2+1
    } else {
      # lines(deltaff,col="black")
    }
    
  }
  count=count+1
}
lines(avg/count2,col="black",lwd=2)

lines(x = abra$z_14$ROI_50$signal)
lines(sgolayfilt(abra$z_14$ROI_50$signal,p=2),col="red")
