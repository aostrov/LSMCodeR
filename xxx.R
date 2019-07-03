

prepForDisplay <- function(image,z,t,verbose=FALSE) {
  file.h5 <- H5File$new(file.path(imageDir,image), mode = "r+")
  imageDataSlice<-file.h5[["imagedata"]]
  on.exit(file.h5$close_all())
  on.exit(imageDataSlice$close(),add=T)
  test.matrix <- makeDFF(aperm((imageDataSlice[c(1:120),,z,,] - 399),c(3,2,1)),
                         backgroundSlices = c(20:40),
                         xyzDimOrder = c(1,2,3),verbose = verbose)
  test.matrix.melt <- melt(data = test.matrix)
  names(test.matrix.melt) <- c("x","y","t","dff")
  temp <- test.matrix.melt[test.matrix.melt$t==t,]
  temp[is.infinite(temp$dff)| is.na(temp$dff),"dff"] <- 0
  temp[temp$dff<quantile(temp$dff,na.rm=TRUE,probs = seq(0,1,0.25))["25%"] | 
         temp$dff>quantile(temp$dff,na.rm=TRUE,probs = seq(0,1,0.25))["75%"],"dff"] <- 0
  return(temp)
}

ggplot(prepForDisplay("/AAFA-gen-A-laser3-SabineSimple/AAFA-gen-A-laser3-SabineSimple.mat",z=6,t=95),
       aes(x,y,fill=dff)) + 
  geom_raster(interpolate = F) + 
  scale_fill_gradientn(colors = viridis(6)) +
  coord_fixed() + scale_y_reverse()

##############
## Analysis ##
##############

# time series stuff
aaia <- readRDS(file.path(LSMCodeRConfig$srcdir,"toDiscardEventually/AAIA-gen_A_laser-3_SabineSimple.mat.1.1.RDS"))
abra <- readRDS(file.path(LSMCodeRConfig$srcdir,"toDiscardEventually/ABRA-gen_A-laser-1-SabineSimple.mat.1.1.RDS"))
aaaa1.1 <- readRDS(file.path(LSMCodeRConfig$srcdir,"toDiscardEventually/AAAA-20190211-SabineSimple-laser_3-SP.mat.1.1.RDS"))

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


# For an analysisDF ROI, can I:
#   1) go back and get the particular nrrd that corresponds to an interesting stimulus
#   2) figure out what part of the image corresponds to the ROI
#   3) replot the ggplot to make this single example stand out 
# 
# I should be able to do this in one or two functions

subsetArbitraryDF <- function(df,row){
  return(df[row,])
}





# cumulative sum

analysisDF <- readRDS(file.path(LSMCodeRConfig$srcdir,"objects",paste("analysisDF",".RDS",sep="")))
animals <- unique(substr(analysisDF$animal,1,3))

cumSumDF <- data.frame()
for (animal in animals) {
  animalAnalysisDF <- analysisDF[grepl(paste("^",animal,sep=""),analysisDF$animal),]
  dffs <- animalAnalysisDF[
    is.finite(animalAnalysisDF$dff.max) & animalAnalysisDF$background.mean > quantile(animalAnalysisDF$background.mean)["25%"],
    "dff.max"
    ]
  names(dffs) <- rownames(animalAnalysisDF[
    is.finite(animalAnalysisDF$dff.max) & animalAnalysisDF$background.mean > quantile(animalAnalysisDF$background.mean)["25%"],
    ])
  dffs.sort <- sort(dffs)
  dffs.cumsum <- cumsum(dffs.sort)
  dffs.cumsum.percent <- dffs.cumsum/(max(dffs.cumsum))
  dffs.df <- data.frame(sorteddFF=dffs.sort,cumSumPercent=dffs.cumsum.percent)
  dffs.df$animal <- animal
  dffs.df$genotype <- unique(as.character(animalAnalysisDF$geno))
  cumSumDF <- rbind(cumSumDF,dffs.df)
}


ggplot(data=dffs.df,aes(sorteddFF,cumSumPercent)) + 
  geom_point() + 
  geom_point(data=dffs.df[findInterval(0.5,dffs.df$cumSumPercent),],color="red",size=2)

ggplot(data=cumSumDF[(cumSumDF$animal!="AAM" |  cumSumDF$animal!="ABJ") &
                       cumSumDF$genotype!="iGABASnFr" & 
                       cumSumDF$cumSumPercent>0.45 & 
                       cumSumDF$cumSumPercent<0.55,],
       aes(sorteddFF,cumSumPercent)) + 
  geom_point(aes(color=genotype)) + 
  xlim(c(0,5)) #+ ylim(c(45,55))

# TODO: maybe get an average for each animal, of each 1 or 10% of the dffs. 
# See how much that differs at 50% from using the full range of numbers. 
# Then use that to average across animals per genotype.
deciles <- quantile(cumSumDF$cumSumPercent,probs = seq(0,1,.1))
bottom <- cumSumDF[cumSumDF$cumSumPercent>deciles[1] & cumSumDF$cumSumPercent<deciles[2],]

decile <- 0
percentDF <- data.frame()
for (animal in animals) {
  animalAnalysisDF2 <- analysisDF[grepl(paste("^",animal,sep=""),analysisDF$animal),]
  animalAnalysisDF2 <- animalAnalysisDF2[is.finite(animalAnalysisDF2$dff.max) & 
                                           animalAnalysisDF2$background.mean > quantile(animalAnalysisDF2$background.mean)["25%"] &
                                           animalAnalysisDF2$dff.max!=0,]
  percents <- quantile(animalAnalysisDF2$dff.max,probs = seq(0,1,.01),names = FALSE,na.rm = TRUE)
  percent.cumsum <- cumsum(percents)/max(cumsum(percents))
  

    
  tempDF <- data.frame(PercentdFF=percents,
                       cumsum=percent.cumsum,
                       percentile=c(0:100),
                       animal=animal,
                       genotype=unique(animalAnalysisDF2$geno))
  percentDF <- rbind(percentDF,tempDF)
}

genotypes=as.character(unique(percentDF$genotype))
averageDFF <- data.frame()
for (genotype in genotypes){
  genos=percentDF[percentDF$genotype==genotype,]
  animals=as.character(unique(genos$animal))
  averageDFF.wide <- data.frame(percentile=c(0:100))
  for (animal in animals) {
    averageDFF.wide <- cbind(averageDFF.wide,genos[genos$animal==animal,"PercentdFF"])
    means <- c()
    SDs <- c()
    for (row in 1:nrow(averageDFF.wide)) {
      means <- c(means,mean(unlist(averageDFF.wide[row,c(2:length(averageDFF.wide))])))
      SDs <- c(SDs,sd(unlist(averageDFF.wide[row,c(2:length(averageDFF.wide))])))
    }
    averageDFF.long <- data.frame(percentile=c(0:100),
                                  means=means,
                                  SDs=SDs,
                                  genotype=genotype)
  }
  averageDFF <- rbind(averageDFF,averageDFF.long)
}


# averageDFF=averageDFF/length(animals)

ggplot(data=percentDF,
       aes(PercentdFF,cumsum)) + 
  geom_point(aes(color=genotype))

## Time Series Nonsense ##

xx=readRDS("/Users/aostrov/projects/R/OrgerLab/LSMCodeR/toDiscardEventually/AAFA-gen-A-laser3-SabineSimple.mat.1.1.RDS")
ts.test <- xx[[1]][[1]]$signal
plot(ts(ts.test))

testData <- subset(df3,X.5>0.2)
maxInstantaneousDerivative <- c()
for (i in 1:nrow(testData)) {
  maxInstantaneousDerivative <- c(maxInstantaneousDerivative,
                                  which.max(diff(unlist(c(testData[i,grepl(x = colnames(df3),pattern = "^X")])))))
}
testData$maxDiff=maxInstantaneousDerivative

zByColor <- data.frame(maxDiff=unique(testData$maxDiff),color=viridis(length(unique(testData$maxDiff))))

plot(x=seq(0,6,0.2),
     y = testData[1,grepl(x = colnames(testData),pattern = "^X")],
     type = "l",ylim=c(-0.05,max(testData[,6:18])),
     col=zByColor[testData$maxDiff[1],"color"])

for (j in 2:nrow(testData)){
  lines(x=seq(0,6,0.2),y=testData[j,grepl(x = colnames(testData),pattern = "^X")],col=as.character(zByColor[testData$maxDiff[j],"color"]))
}



plot(x=seq(0,6,0.2),
     y = testData[1,grepl(x = colnames(testData),pattern = "^X")],
     type = "l",ylim=c(-0.05:max(testData[,6:18])))

for (j in 1:nrow(testData)){
  lines(x=seq(0,6,0.2),y=subset(testData,maxDiff==3)[j,grepl(x = colnames(testData),pattern = "^X")])
}

for (something in zByColor$maxDiff) {
  xxx <- subset(testData,maxDiff==something)
  plot(x=seq(0,6,0.2),
       y = (xxx[1,grepl(x = colnames(testData),pattern = "^X")] / 
              max(xxx[1,grepl(x = colnames(testData),pattern = "^X")])),
       type = "l",ylim = c(-0.1,1.1),
       xlab = "Seconds", ylab = "Normalized dFF",main = paste("maxDiff of", something))
  
  for (j in 1:nrow(xxx)){
    lines(x=seq(0,6,0.2),
          y=xxx[j,grepl(x = colnames(testData),pattern = "^X")] / 
            max(xxx[j,grepl(x = colnames(testData),pattern = "^X")]),
          col=as.character(zByColor[xxx$maxDiff[j],"color"]))
  }
  
  
}


testdata.melt <- melt(testData,id.vars = c("z","roi","maxDiff"),measure.vars = grep(x = colnames(testData),pattern = "^X"))

# aligning the peaks of different traces

# reference nrrd

rigidAlignment <- function(df,subsetByThreshold=T,threshold=0.2,outDir="~/Desktop/tmp/",refRow=61){
  if (subsetByThreshold){
    df <- subset(df,X.5>0.2)
  }
  fixed=write.nrrd(x = as.matrix(df[refRow, grepl(x = colnames(df),pattern = "^X")]) , 
                   file = file.path(outDir,"referenceTrace_61.nrrd"))
  # unlink(fixed)
  registeredDF <- data.frame()
  for (i in 1:nrow(df)) { 
    moving=write.nrrd(x = as.matrix(df[i, grepl(x = colnames(df),pattern = "^X")]) , 
                      file = file.path(outDir,"moving.nrrd"))
    print(paste("working on row",i))
    system("sh ~/Desktop/2dReg-test.sh",ignore.stdout = TRUE)
    registered=c(read.nrrd(file=file.path(outDir,"reformatted.nrrd")))
    registeredDF <- rbind(registeredDF,registered)
  }
  return(registeredDF)
  
}


testData <- subset(df3,X.5>0.2)

fixed=write.nrrd(x = as.matrix(testData[61, grepl(x = colnames(testData),pattern = "^X")]) , file = "~/Desktop/tmp/referenceTrace_61.nrrd")
# unlink(fixed)

registeredDF <- data.frame()
for (i in 1:nrow(testData)) { #1:nrow(df3)
  moving=write.nrrd(x = as.matrix(testData[i, grepl(x = colnames(testData),pattern = "^X")]) , 
                    file = "~/Desktop/tmp/moving.nrrd")
  print(paste("working on row",i))
  system("sh ~/Desktop/2dReg-test.sh",ignore.stdout = TRUE)
  registered=c(read.nrrd(file="~/Desktop/tmp/reformatted.nrrd"))
  registeredDF <- rbind(registeredDF,registered)
}

plot(x=seq(0,6,0.2),
     y = registeredDF[1,],
     type = "l",ylim=c(-0.05,max(registeredDF[,6:18])))

for (j in 2:nrow(testData)){
  lines(x=seq(0,6,0.2),y=registeredDF[j,])
}
lines(x=seq(0,6,0.2),y=colMeans(registeredDF),col="green",lwd=2)
lines(x=seq(0,6,0.2),y=colMeans(testData[,grepl("X",colnames(testData))]),col="magenta",lwd=2)



plot(x=seq(0,6,0.2),y=colMeans(registeredDF),col="green", ylab = "dF/F", xlab = "Seconds", main = "Comparison of average traces before and after rigid alignment", type = "l")

# try some groupwise averages to get a better fixed image...
testdata.sample <- sample.dataframe(testData,50)
for (nrrd in 1:nrow(testdata.sample)){
  write.nrrd(x = as.matrix(testdata.sample[nrrd, 
                                           grepl(x = colnames(testdata.sample),
                                                 pattern = "^X")
                                           ]), 
             file = paste("~/Desktop/tmp/forTemplate_",nrrd,".nrrd",sep=""))
  
}
plot(x=seq(0,6,0.2),
      y=c(read.nrrd("/Users/aostrov/Desktop/tmp/reformatted_template0.nii-1.nrrd")),
      col="green")
for (j in 1:nrow(testdata.sample)){
  lines( x=seq(0,6,0.2),y=testdata.sample[j,grepl("X",colnames(testdata.sample))] ) 
}
