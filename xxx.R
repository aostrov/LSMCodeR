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

######################################
##### more parallel, less clever #####
######################################


# for (myFile in sample(dir(imageDir,patt="[A-Z]{4}-[[:graph:]]",full=T,rec=TRUE)[1:2],
#                       length(dir(imageDir,patt="[A-Z]{4}-[[:graph:]]",full=T,rec=TRUE)[1:2]))) {
for (myFile in dir(imageDir,patt="[A-Z]{4}-[[:graph:]]",full=T,rec=TRUE)[1:2]) {
  rdsFile <- file.path(LSMCodeRConfig$srcdir,"objects",basename(myFile))
  lockFile <- sub(".mat",".lock",rdsFile)
  matFileCode <- substring(basename(myFile),1,4)
  if (file.exists(rdsFile) | file.exists(lockFile)) {
    print(paste(myFile,"exists or is being worked on. Skipping"))
  } else {
    cat("some text here",file=lockFile)

    print("Making ROI list")
    matFileROIListByZ <- list()
    for (tectumROIforZ in 1:nrow(tempDF <- tectumROIs[as.character(tectumROIs$matfile)==matFileCode,])) {
      print(tectumROIforZ)
      matFileROIListByZ[[paste("z",tempDF[tectumROIforZ,"z"],sep="_")]] = getROIs(roiEdgeLength = roiEdgeLength,
                                                                                  x=tempDF[tectumROIforZ,"x"],
                                                                                  y=tempDF[tectumROIforZ,"y"],
                                                                                  w=tempDF[tectumROIforZ,"w"],
                                                                                  h=tempDF[tectumROIforZ,"h"])
    }
    
    statisticsList <- list()
    
    print("Extracting raw data and useful statistics")
    
    # open the connection to the .mat file
    file.h5 <- H5File$new(myFile, mode = "r")
    # save a new variable that contains the image data
    # this step can be avoided, and might improve performance on very large datasets
    imageDataSlice<-file.h5[["imagedata"]]
    file.h5$close()
    # file.h5$close()
    # data is in the form:
    # [stacks,(channels),slices,rows, columns]
    imageDataSlice.dims <- imageDataSlice$dims
    names(imageDataSlice.dims) <- c('t','c','z','y','x')
    print(imageDataSlice.dims)
    
    for (stimulation in names(stimulusParametersList[[matFileCode]][[1]])){ # for a matfile get the raw data as Z plane and ROI
      # foreach(stimulation=names(currentStimulusParameters)) %dopar% {
      cat(".")
      print(stimulusParametersList[[matFileCode]][[1]][[stimulation]]$start)
      rawDataForMatFileByROIs <- mapply( # outer lapply breaks down roiList by z-plane
        function(roiList,z) 
        {
          lapply( # inner lapply gets the raw statistics for each ROI
            roiList,
            getROIsRawDataFromHDF5.lapply,
            z,
            hdf5Image.mat = imageDataSlice,
            frame.start = stimulusParametersList[[matFileCode]][[1]][[stimulation]]$start, # 750,
            frame.end = stimulusParametersList[[matFileCode]][[1]][[stimulation]]$end,
            offset = 399
          )
        },
        roiList=matFileROIListByZ,z=as.integer(sub("z_","",names(matFileROIListByZ))),SIMPLIFY = F
      )
      # names(rawDataForMatFileByROIs) <- names(matFileROIListByZ)
      # saveRDS(rawDataForMatFileByROIs,file=file.path(LSMCodeRConfig$srcdir,"objects",paste(stimulation,".RDS",sep="")),compress = T)
      # statisticsList[[stimulation]]<- getUsefulStatisticsByROI(rawDataForMatFileByROIs,matFileROIListByZ,
      #                                                          analysisWindow=c(150:750),backgroundWindow=c(0:100))
      cat("+")
      z_planes=attr(stimulusParametersList[[matFileCode]][[1]],"imageDimensions")[['z']]
      statisticsList[[stimulation]]<- lapply(
        rawDataForMatFileByROIs,getUsefulStatisticsByROI,matFileROIListByZ,
        analysisWindow=c((900/z_planes):(1500/z_planes)),
        backgroundWindow=c((700/z_planes):(900/z_planes))
      )
      
      attr(statisticsList[[stimulation]],"ROI_Location") <- c(frame.start=stimulusParametersList[[matFileCode]][[1]][[stimulation]]$start,
                                                              frame.end = stimulusParametersList[[matFileCode]][[1]][[stimulation]]$end)
      attr(statisticsList[[stimulation]],"StimulusBlock") <- c(stimulus=stimulusParametersList[[matFileCode]][[1]][[stimulation]]$stimulus,
                                                               block=stimulusParametersList[[matFileCode]][[1]][[stimulation]]$block)
      attr(statisticsList[[stimulation]],"matFile") <- myFile

    }
    matList[[substr(stimulation,1,4)]] <- statisticsList
    saveRDS(matList,file=matListRDS,compress = TRUE)
    unlink(lockFile)
    # clean up
    file.h5$close()
    imageDataSlice$close()
    
  } 
  
  # stimulusParametersList[[matFileCode]] <- animal
  # saveRDS(stimulusParametersList,file=stimParListRDS,compress = TRUE)
}
