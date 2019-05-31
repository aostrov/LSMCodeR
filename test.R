# stimParListRDS <- file.path(LSMCodeRConfig$srcdir,"objects","stimParListForParallel.RDS")
# stimulusParametersList <- readRDS(file=stimParListRDS)
# analysisDF <- readRDS(file=file.path(LSMCodeRConfig$srcdir,"objects",paste("analysisDF",".RDS",sep="")))

getFramesFromStimParamListFromAnalysisDF <- function(row,file,writeEntireNrrd=T,...) {
  animal <- as.character(analysisDF[row,"animal"])
  print(animal)
  stim <- analysisDF[row,"stimulus"]
  z_plane <- as.integer(sub("z_","",analysisDF[row,"z_plane"]))
  roi <- as.integer(sub("ROI_","",analysisDF[row,"roi"]))
  
  matfile <- dir(imageDir,pattern = paste("^",animal,sep=""),full = T, recursive = T)
  print(matfile)
  myStim <- stimulusParametersList[[animal]][[grep(stim,names(stimulusParametersList[[animal]]))]]
  myROIs <- tectumROIs[tectumROIs$matfile==animal & tectumROIs$z==z_plane,]
  rois <- getROIs(roiEdgeLength = roiEdgeLength,x=myROIs$x,y=myROIs$y,w=myROIs$w,h=myROIs$h)
  
  y_dims <- rois[[roi]]$yRange
  x_dims <- rois[[roi]]$xRange
  file.h5 <- H5File$new(matfile,mode = "r")
  imageDataSlice<-file.h5[["imagedata"]]
  on.exit(file.h5$close())
  on.exit(imageDataSlice$close(),add = T)
  imageDataSlice.dims <- imageDataSlice$dims
  print(imageDataSlice.dims)
  names(imageDataSlice.dims) <- c('t','c','z','y','x')
  if (writeEntireNrrd) {
    if (imageDataSlice.dims[['z']]==1){
      print("downsampling")
      write.nrrd(aperm(imageDataSlice[seq(myStim$start,myStim$end,by=20),,z_plane,,],c(3,2,1)),file=file,...)
    } else {
      print("Using full image.")
      write.nrrd(aperm(imageDataSlice[myStim$start:myStim$end,,z_plane,,],c(3,2,1)),file=file,...)
    }
    writeRectForROI(file=paste(file,"pdf",sep = "."),imageDataSlice.dims[['x']],imageDataSlice.dims[['y']],x_dims,y_dims)
    
    ggplot(
      analysisDF[sampleAnalysisDFRownamesByGenotype(500,genotype = analysisDF[row,"geno"],20,0.75),],
      aes(background.mean,dff.max)
    ) +
      geom_jitter(aes(color=animal,size=snr.mean),alpha=0.1) +
      facet_wrap(~geno, ncol=3) +
      ylim(c(-1,8)) + xlim(c(-10,320)) + scale_color_viridis_d() +
      geom_point(data=analysisDF[row,], aes(background.mean,dff.max,
                                                      fill = "black",size = snr.mean)) +
      theme(legend.position = "none")
    ggsave(paste(file,"ggplot.pdf",sep = "."),width = 5,height = 5,units = "cm")
    
  } else {
    write.nrrd(aperm(imageDataSlice[myStim$start:myStim$end,,z_plane,y_dims,x_dims],c(3,2,1)),file=file,...)
  }
  # return(paste())
  
}

writeRectForROI <- function(file,fullImageWidth,fullImageHeight,roiXRange,roiYRange) {
  pdf(file=file,
      title = "",
      width = (7*fullImageWidth/fullImageHeight),
      height = 7*(fullImageHeight/fullImageWidth))
  plot(x = c(1,fullImageWidth), y = c(1,fullImageHeight), 
       xlim = c(1,fullImageWidth), ylim = c(fullImageHeight,1),
       type = "n",ann = F,axes = F, frame.plot = T)
  rect(xleft = roiXRange[1],xright = roiXRange[length(roiXRange)],
       ytop = roiYRange[1],ybottom = roiYRange[length(roiYRange)])
  dev.off()
}

# TODO: Try to get an outlier, something with a large dff but small background, for each animal.
# OR, I get the average response for an animal, and try to find a response that is pretty close
# to that.
sampleAnalysisDFRownamesByAnimal <- function(numberOfSamples,animal,backgroundFloor,dffMax){
  rows <- rownames(analysisDF[analysisDF$background.mean>backgroundFloor & 
                               analysisDF$animal==animal & 
                               analysisDF$dff.max > dffMax,])
  if (numberOfSamples > length(rows)) {
    warning(paste(numberOfSamples," is greater than the number of rows (",length(rows),"). Setting the number of rows to 1/3 of the total rows",sep=""))
    numberOfSamples <- length(rows)*0.33
    
  }
  return(sample(rows,numberOfSamples))
}

sampleAnalysisDFRownamesByGenotype <- function(numberOfSamples,genotype,backgroundFloor,dffMax){
  rows <- rownames(analysisDF[analysisDF$background.mean>backgroundFloor & 
                                analysisDF$geno==genotype & 
                                analysisDF$dff.max > dffMax,])
  if (numberOfSamples > length(rows)) {
    warning(paste(numberOfSamples," is greater than the number of rows (",length(rows),"). Setting the number of rows to 1/3 of the total rows",sep=""))
    numberOfSamples <- length(rows)*0.33
  }
  return(sample(rows,numberOfSamples))
}


animals <- as.character(unique(analysisDF[with(analysisDF,background.mean>50 & dff.max > 0.75),"animal"]))
# examples <- c()
for (animal in animals) {
  print(animal)
  rowname.sample <- sampleAnalysisDFRownamesByAnimal(1,animal,backgroundFloor = 50,dffMax = 0.75)
  getFramesFromStimParamListFromAnalysisDF(rowname.sample,
                                           file=file.path("C:/Users/Aaron/Desktop/exampleTraces/",
                                                          paste(animal,
                                                                "_",
                                                                rowname.sample,
                                                                ".nrrd",sep="")))
  
  examples <- c(examples,rowname.sample)
}
names(examples) <- animals


# A monster call. It tries to find use the average response per 'animal' to find the actual ROI so that I can
# single out that ROI and use it for display purposes
animalKeyDF <- data.frame()
for (animalKey in as.character(animalSummaryDF2$animalKey)) {
  analysisDF.sorted <- sort(analysisDF[grep(animalKey,analysisDF$animal),"dff.max"])
  analysisDF.sorted <- analysisDF.sorted[-is.infinite(analysisDF.sorted)]
  ak <- analysisDF[ # I want to pull a full row from analysisDF
    analysisDF$dff.max==analysisDF.sorted[ # I am subsetting by dff.max, and by Animal trial
      findInterval(animalSummaryDF2[animalSummaryDF2$animalKey==animalKey,"dff.2sd"], # find the interval lower bound for the average response
                   analysisDF.sorted) # I use the sorted analysisDF for the animal 
      ],
    ]
  animalKeyDF <- rbind(animalKeyDF,ak)
}

for (row in rownames(animalKeyDF)) {
  getFramesFromStimParamListFromAnalysisDF(row,
                                           file=file.path("C:/Users/Aaron/Desktop/exampleTraces/",
                                                          paste(animalKeyDF[row,"animal"],
                                                                "_",
                                                                row,
                                                                "-represenative2sd.nrrd",sep="")))
  
}

# > examples
# AAAA            AAAB            AABA            AABB            AABC            AACA            AACB 
# "ROI_8510"      "ROI_9421"     "ROI_33142"      "ROI_2803"     "ROI_17184"      "ROI_1795"     "ROI_57216" 
# AADA            AAEA            AAFA            AAGA            AAGB            AAHB            AAIA 
# "ROI_407112"   "ROI_1526142"     "ROI_34627"   "ROI_6110164"     "ROI_70208"     "ROI_69038"   "ROI_1341206" 
# AAJA            AAKA            AALA            AAMA            AANA            AAOA            AAPA 
# "ROI_921187"    "ROI_711906"    "ROI_130567"  "ROI_31101910"    "ROI_136280"    "ROI_942315"   "ROI_5916109" 
# AAQA            AARA            AARB            AASA            AATA            AAUA            AAWA 
# "ROI_10491811"    "ROI_506015"    "ROI_201900"   "ROI_5118615"   "ROI_1820220"   "ROI_2102017"    "ROI_248253" 
# AAXA            AAXB            AAYA            AAZA            ABAA            ABBA            ABCA 
# "ROI_19110147"    "ROI_179600"    "ROI_512399"   "ROI_6761225"  "ROI_11051523"   "ROI_1143424"   "ROI_1007425" 
# ABDA            ABEC            ABFA            ABHA            ABIA            ABKA            ABLA 
# "ROI_486331"   "ROI_4621328"    "ROI_306246"   "ROI_1967145"    "ROI_972275"  "ROI_36121028" "ROI_175101112" 
# ABMA            ABNA            ABOA 
# "ROI_1095245"   "ROI_1393251"   "ROI_2156138" 

