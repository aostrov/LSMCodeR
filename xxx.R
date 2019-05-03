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

for (file in multiplaneFiles){
  writeNrrdForROISelection(file,"C:/Users/Aaron/Desktop/nrrdOrder/")
}

reorderImageSlices <- c(2,1,20:5)
file.h5$close_all()
imageDataSlice$close()

someVar <- c()
for (i in ((1:20)-1) ){
  someVar <- c(someVar,length(seq(from=(13+i),to=(13+1800),by=20)))
}
first <- seq(from=(13+1),to=(13+1800),by=20)

# write test images to make ROIs
dirByDate <- function(date="2019-03-05",directory=imageDir){
  subsettedDirByDate <- dir(directory,full=T)[
    grepl(date,file.info(dir(directory,full=T))$mtime)
    ]
  return(subsettedDirByDate)
}
recent <- dirByDate()
recentPhysio <- recent[!grepl("Anatomy",recent)]

for (image in recentPhysio) {
  print(image)
  writeNrrdForROISelection(file.path(image,paste(basename(image),"mat",sep=".")),"C:/Users/Aaron/Desktop/nrrdOrder/")
}

# save the relevant data from each ROI.
# full data still exists as .mat files
# but it might be nice to have the background
# and signal areas saved externally
test <- lapply(myFile, function(x) {
  lapply(x,
         function(y) {
           background <- apply(y[,,750:850],3,mean)
           signal <- apply(y[,,900:1200],3,mean)
           return(list(background=background,signal=signal))
         })
})
##############
## Analysis ##
##############
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
      analysisDF.zplane <- matList[[k]][[j]][[i]][,c("dff.mean","background.mean")]
      analysisDF.zplane$z_plane <- as.factor(z_plane)
      analysisDF.stimulus <- rbind(analysisDF.stimulus,analysisDF.zplane)
    }
    analysisDF.stimulus$stimulus <- as.factor(stimulus)
    analysisDF.animals <- rbind(analysisDF.animals,analysisDF.stimulus)
  }
  analysisDF.animals$animal <- as.factor(animal)
  analysisDF <- rbind(analysisDF,analysisDF.animals)
}


ggplot(analysisDF.subset,aes(background.mean,dff.mean)) + 
  geom_jitter(aes(color=animal))

# plot a raster of the image
ggplot(matList[[1]][[5]][[1]],aes(xpos,ypos,fill=dff.mean)) + 
  geom_raster(interpolate = F) + 
  scale_fill_gradientn(colors = jet(20)) +
  coord_fixed() + scale_y_reverse()

# background > 25 and dff > 0.05
analysisDF.subset <- subset(analysisDF,background.mean>25 & dff.mean>0.05)

for (plane in 1:20) {
  yyy=lapply(
    currentStimulusParameters, # list being worked on
    processSingleStimulus.lapply, # function doing the work
    outputType="dff",writeNRRD=T,downSampleImage=T,resizeFactor=2,image=imageDataSlice,z_plane=plane # other variables
    )
}

# time series stuff
aaia <- readRDS(file.path(LSMCodeRConfig$srcdir,"toDiscardEventually/AAIA-gen_A_laser-3_SabineSimple.mat.1.1.RDS"))

make_norm_dist <- function(x, mean, sd){
  norm = c()
  for (i in seq(x)) {
    norm = c(norm, 1.0/(sd * sqrt(2*pi)) * exp(-(i - mean)^2/(2*sd^2)) )
  }
  return(norm)
}
test=make_norm_dist(1000,500,100)
x=sort(sample(length(test),200))
y=test[x]
fit4 <- lm(y~poly(x,4,raw=TRUE))

hm=test[which.max(test)]/2

xmax=which.max(test[200:800])
testdf=data.frame(y=test[200:800],x=seq(length(test[200:800])))
x1 <- testdf$x[testdf$x < xmax][which.min(abs(testdf$y[testdf$x < xmax]-max(testdf$y)/2))]
x2 <- testdf$x[testdf$x > xmax][which.min(abs(testdf$y[testdf$x > xmax]-max(testdf$y)/2))]

points(c(x1, x2), c(testdf$y[testdf$x==x1], testdf$y[testdf$x==x2]), col="red")

rois <- aaia[["z_5"]]
getTopROIperZ <- function(rawSubsettedDataByROI) {
  # rois <- aaia[["z_5"]]
  rois.sort <- sort(sapply(rawSubsettedDataByROI, function(x) max(x$signal)))
  rois.sort.max <- which.max(rois.sort)
  roi.sort.max.background <- mean(rawSubsettedDataByROI[[names(rois.sort.max)]]$background)
  roi.sort.max.dff <- (rawSubsettedDataByROI[[names(rois.sort.max)]]$signal -
                         roi.sort.max.background) / 
                          roi.sort.max.background
  # roi.sort.max.dff.xmax <- which.max(roi.sort.max.dff)
  return(roi.sort.max.dff)
}



rois.sort <- sort(sapply(rois, function(x) max(x$signal)))
rois.sort.max <- which.max(rois.sort)
roi.sort.max.background <- mean(rois[[names(rois.sort.max)]]$background)
roi.sort.max.dff <- (rois[[names(rois.sort.max)]]$signal-roi.sort.max.background)/roi.sort.max.background
roi.sort.max.dff.xmax <- which.max(roi.sort.max.dff)
fit4=lm(roi.sort.max.dff[1:roi.sort.max.dff.max] ~ poly(1:length(roi.sort.max.dff[1:roi.sort.max.dff.max]),4,raw=TRUE))

fourth_order <- function(newdist, model) {
  coefs <- coef(model)
  #y = d + cx + bx^2 + ax^3
  res <- coefs[1] + (coefs[2] * newdist) + (coefs[3] * newdist^2) + (coefs[4] * newdist^3) + (coefs[4] * newdist^4)
  return(res)
}

# try to figure out where the start of the rise occurs
zzz=lapply(aaia,getArbitraryTopROIsPerZ)
zzz.max <- lapply(zzz,function(x) lapply(x,which.max))
zzz.hist=hist(sapply(zzz,function(x) which.max(diff(x))),breaks = length(zzz))
timeRiseStarts <- zzz.hist$breaks[which.max(zzz.hist$counts)]

# maybe I want to try to compare each of the plots to a model???
stimParListRDS <- file.path(LSMCodeRConfig$srcdir,"objects","stimParListForParallel.RDS")
stimulusParametersList <- readRDS(file=stimParListRDS)
myFile=file.path(imageDir,"AABA-20190211-SabineSimple-laser_3-GCaMP7b-SP","AABA-20190211-SabineSimple-laser_3-GCaMP7b-SP.mat")
file.h5 <- H5File$new(myFile, mode = "r")
imageDataSlice<-file.h5[["imagedata"]]

processSingleStimulus.lapply(myList=stimulusParametersList[["AABA"]][[3]],outputType = "dff",writeNRRD = T,image=imageDataSlice,resizeFactor = 2,downSampleImage=T)
