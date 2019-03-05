mip<-intensityProjection(colMeansAverage)
write.nrrd(mip,file="C:/Users/Aaron/Desktop/mip2slices.nrrd")

################
## Physiology ##
################


outDir<-"F:/Imaging/GCaMP7_tests/outputNRRDs/"
physioDirs <- dir("F:/Imaging/GCaMP7_tests/20181204-g7",patt='SP',full=T)
physiologyFiles <- dir("F:/Imaging/GCaMP7_tests/20181204-g7/",patt="[A-Z]{4}-[[:graph:]]*.mat",full=TRUE,rec=T)

physiologyFiles <- dir("F:/Imaging/GCaMP7_tests/20181204-g7/",patt="[A-Z]{4}-[[:graph:]]*SP",full=T,rec=FALSE)
special <- "examples"
outputTypes=c("snr","raw","dff")
outputType=("raw")
for (physioDir in physiologyFiles){
  lsmLogFile <- dir(file.path(physioDir,"logs"),full=T,rec=F,patt="lsmlog_")
  stimLogFile <- dir(file.path(physioDir,"logs"),full=T,rec=F,patt="stimlog_")
  myFile <- dir(physioDir,full=T,rec=F,patt=".mat")
  outDirSubDir <- paste(basename(physioDir),"_",outputType,special,"/",sep="")
  # source(file.path(LSMCodeRConfig$srcdir,"parseImageForGreenFlashAndLSMandStimLogs.R"))
  # source(file.path(LSMCodeRConfig$srcdir,"roiBasedAnalysis.R"))
  source(file.path(LSMCodeRConfig$srcdir,"physiologyScript.R"))
}

#######################
## multiplane images ##
#######################

file.h5 <- H5File$new("/Volumes/TranscendJD/Work/AAFA-gen-A-laser3-SabineSimple/AAFA-gen-A-laser3-SabineSimple.mat", 
                      mode = "r+")
imageDataSlice<-file.h5[["imagedata"]]
test.matrix <- imageDataSlice[2,,12,,]
test.matrix.melt <- melt(data = test.matrix)
ggplot(test.matrix.melt,aes(Var2,Var1,fill=value)) + 
  geom_raster(interpolate = F) + 
  scale_fill_gradientn(colors = jet(2000)) +
  coord_fixed() + scale_y_reverse()

write.nrrd(aperm(imageDataSlice[1200:1600,,13,,],c(3,2,1)),file = "C:/Users/Aaron/Desktop/multiplaneNrrd.nrrd",dtype = "short")

flashMean <- mean(imageDataSlice[1:100,,,1:200,1:200])
flashSD <- sd(imageDataSlice[1:100,,,1:200,1:200])
flashSlices <- c()
for (i in 1:20){
  flashSlices <- c(flashSlices,mean(imageDataSlice[1:10,,i,1:200,1:200]))
}
max.flash <- which.max(flashSlices)

slicesWithBigNumbers.begin <- c()
for (i in 1:20){
  if (mean(imageDataSlice[1:10,,i,1:200,1:200])>(flashMean+2*flashSD)){
    slicesWithBigNumbers.begin <- c(slicesWithBigNumbers.begin,i)
  }
}
startOfStimulations <- max(slicesWithBigNumbers.begin)
endOfFirstGreenFlash <- startOfStimulations + 1


count=0
slice.identity <- data.frame(slice=integer(0),z_plane=integer(0),time=integer(0))
for (time in 1:1650){
  for (z_plane in 1:20){
    count=count+1
    slice.identity <- rbind(slice.identity,c(count,z_plane,time))
  }
}
colnames(slice.identity) <- c("slice","z_plane","time")
slice.transitions <- slice.identity[lsm.transition.frames,]

}

reorderImageSlices <- c(2,1,20:5)
file.h5$close_all()
imageDataSlice$close()

someVar <- c()
for (i in ((1:20)-1) ){
  someVar <- c(someVar,length(seq(from=(13+i),to=(13+1800),by=20)))
}
first <- seq(from=(13+1),to=(13+1800),by=20)

