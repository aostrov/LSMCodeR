# Basic functions for dealing with microscopy data

readLogFileData<-function(logFile){
  lines<-readLines(logFile,n=2)
  magic<-"<!-- LSM acquisition log file | Orger's Lab | Champalimaud Center for the Unknown -->"
  if (lines[2]!=magic){
    stop("This doesn't appear to be a log file from the LSM!")
  }
  lines<-readLines(logFile)
  framedata<-lines[(grep("<framedata>",lines)+1):(grep("</framedata>",lines)-1)]
  bigFrame<-lapply(strsplit(framedata,"\\t"), function(x) as.numeric(x[c(1,2,9)]))
  bigDataFrame<-data.frame(matrix(unlist(bigFrame), nrow=length(bigFrame), byrow=T))
  colnames(bigDataFrame)<-c("frame","time","otherTime")
  return(bigDataFrame)
}

readLogFileMetaData<-function(logFile){
  lines<-readLines(logFile,n=2)
  magic<-"<!-- LSM acquisition log file | Orger's Lab | Champalimaud Center for the Unknown -->"
  if (lines[2]!=magic){
    stop("This doesn't appear to be a log file from the LSM!")
  }
  lines<-readLines(logFile)
  framedata<-lines[1:10]
  framedata<-gsub("<[[:graph:]]*?>","",framedata)[c(3,5,6,7,8)]
  names(framedata)<-c("Date","exposure","numSheets","verticalSpan","frameCount")
  return(framedata)
}


makeDFF<-function(imageStack,backgroundSlices=c(50:101),xyzDimOrder=c(3,2,1)){
  zdim<-dim(imageStack)[xyzDimOrder[3]]
  xdim<-dim(imageStack)[xyzDimOrder[1]]
  ydim<-dim(imageStack)[xyzDimOrder[2]]
  testSliceBackground<-apply(imageStack[,,backgroundSlices],c(xyzDimOrder[1],xyzDimOrder[2]),mean)
  testSliceBackgroundSubtract<-apply(imageStack[,,],3,function(x) x-testSliceBackground)
  dim(testSliceBackgroundSubtract)<-c(xdim,ydim,zdim)
  testSliceDFF<-apply(testSliceBackgroundSubtract, 3, function(x) x/testSliceBackground)
  dim(testSliceDFF)<-c(xdim,ydim,zdim)
  return(testSliceDFF)
}

makeDFFwithBaselineSubtraction <- function(imageStack,baselineSlices=c(50:101),backgroundSlices=c(50:101),xyzDimOrder=c(3,2,1)){
  zdim<-dim(imageStack)[xyzDimOrder[3]]
  xdim<-dim(imageStack)[xyzDimOrder[1]]
  ydim<-dim(imageStack)[xyzDimOrder[2]]
  
  outsideFishROIBackground <- mean(apply(imageStack[1:10,1:10,baselineSlices],3,mean))
  testSliceBackgroundSubtractROI <- (imageStack[,,] - outsideFishROIBackground)
  dim(testSliceBackgroundSubtractROI)<-c(xdim,ydim,zdim)
  
  testSliceBackground<-apply(testSliceBackgroundSubtractROI[,,backgroundSlices],c(xyzDimOrder[1],xyzDimOrder[2]),mean)
  testSliceBackgroundSubtract<-apply(testSliceBackgroundSubtractROI[,,],3,function(x) x-testSliceBackground)
  dim(testSliceBackgroundSubtract)<-c(xdim,ydim,zdim)
  
  testSliceDFF<-apply(testSliceBackgroundSubtract, 3, function(x) x/testSliceBackground)
  dim(testSliceDFF)<-c(xdim,ydim,zdim)
  
  return(testSliceBackgroundSubtract)
}

resizeImage = function(img, new_width, new_height, func=c("spline","linear")) {
  func=match.arg(func)
  if (func=="linear") func<-"approx"
  func=match.fun(func)
  new_img = apply(img, 2, function(y) {return (func(y, n = new_height)$y)})
  new_img = t(apply(new_img, 1, function(y) {return (func(y, n = new_width)$y)}))
  
  new_img[new_img < 0] = 0
  new_img = round(new_img)

  return (t(new_img))
}

downSampleMean2x2<-function(imageData){
  dimImageData<-dim(imageData)
  print(dimImageData)
  smallArray=c()
  for (x in seq(1,dimImageData[2],by=2)){
    for (y in seq(1,dimImageData[1],by=2)){
      # print(paste(y,x))
      smallArray<-c(smallArray,mean(c(imageData[y,x],imageData[y,(x+1)],imageData[(y+1),x],imageData[(y+1),(x+1)])))
    }
  }
  dim(smallArray)<-dimImageData/2
  return(t(round(smallArray)))
}

###############################
## Average of Anatomy Stacks ##
###############################

# ... can be used to pass a digits = ... to round()
frameAverageForAnatomyStacks<-function(imageDataArray,
                                       dimensionOrder=c(3,2,1),
                                       reorderImageSlices=NULL,
                                       roundOutput=F,
                                       digits=1) {
  averageSlices<-colMeans(imageDataArray,dims=1)
  averageSlices<-aperm(averageSlices,dimensionOrder)
  if (!is.null(reorderImageSlices)){
    averageSlices<-averageSlices[,,reorderImageSlices]
  }
  if (roundOutput){
    averageSlices <- round(averageSlices,digits)
  }
  return(averageSlices)
}

doesReorder<-function(testImageArray,reorderImageSlices=c(84:100,1:83)){
  myAnatomyStack<-frameAverageForAnatomyStacks(testImageArray) - 
    frameAverageForAnatomyStacks(testImageArray,reorderImageSlices = reorderImageSlices)
  if (any(myAnatomyStack)) { 
    return(1) 
  }
  return(0)
}

##############################
## max intensity projection ##
##############################

# create an intensity projection of a 3D image 
# supports maximum, minimum, and mean projections,
# as well as standard deviation. Sd projections are comparatively slow
#
# Assumes that the images have been read in as a multidimensional vector
# The default is to assume that x is the 3rd column, and y is the 2nd column,
# but this default can be overwritten.
# Returns a 2D vector of equal length to x*y, where the z dimension has been collapsed.
#
# If a multi-frame image is passed as an array form [frames,?,slices,rows,columns]
# apply() is used to recursively call intensityProjetion.
# This returns a 3D vector that has dimensions x*y*f with the z dimension collapsed.
#
# At the moment this function does not support direct processing of H5D objects,
# these objects must first be coerced to a 4D array before being passed to the function. 
intensityProjection<-function(imageDataArray,rowColumn=c(3,2),projectionType=c("max","min","mean","sd"),fourD=FALSE) {
  if (fourD) {
    print("recursive stage reached")
    dims<-dim(imageDataArray)
    ip<-apply(imageDataArray,1,intensityProjection,rowColumn=rowColumn,projectionType=projectionType)
    dim(ip)<-c(dims[4],dims[3],dims[1])
    return(ip)    
  } else {
    print("deep state reached")
    print(dim(imageDataArray))
    projectionType=match.arg(projectionType)
    apply(imageDataArray, rowColumn, match.fun(projectionType))
  }
}

#########################
## Display of analysis ##
#########################

# When presented with a list of values over time, make a silly animated plot
# using ggplot
#
# takes a two column CSV file (preferably exported directly out of imageJ),
# with a first column X that has frames from a time series, and second
# column Y that has the grey level indicated
# also takes an outdir, where the animated gif will be saved
# Warning: requires imagemagick to be installed, and will not warn if it is
# not installed and in the path, so I have no idea how this will work on
# windows!!!
# There are also lots of other things that can be tweaked when either making
# the ggplots (size and color of the dot, output file type, x&y labels),
# or when making the gif (duration of frames, looping characteristics,
# and output size). 
# ... also provides options for read.csv().
#
# At the moment this saves each frame to a temporary directory and then deletes
# the directory at the end of the function.
# This should create Gifs of up to 9999 frames without getting things out
# of order, though I haven't tested that yet.
animatedTimeSeries <- function(timeSeriesFile,
                               outdir,
                               outFileName=NULL,
                               ylab="dF/F",
                               xlab="Frame (100ms/frame)",
                               color="red",
                               size=2,
                               fps=10,
                               outputX=480,
                               outputY=270,
                               resolution=96,
                               returnGIF=F,
                               ...) {
  csv<-read.csv(timeSeriesFile,...)
  outfile <- ifelse(is.null(outFileName),
         file.path(outdir,paste(basename(timeSeriesFile),".gif",sep="")),
         file.path(outdir,paste(outFileName,".gif",sep="")))
  img <- image_graph(outputX, outputY, res = resolution)
  out <- sapply(c(1:nrow(csv)),makeGIFWithMagick,csv,ylab=ylab,xlab=xlab,color=color,size=size)
  dev.off()
  animation <- image_animate(img, fps = fps)
  image_write(animation, outfile)
  if (returnGIF)  return(animation)
}

makeGIFWithMagick <- function(frameNumber,dataCSV,ylab="Y",xlab="X",color="yellow",size=5) {
  test <- ggplot(data=dataCSV,aes(X,Y)) +
    geom_point() + 
    ylab(ylab) +
    xlab(xlab) +
    geom_point(data=dataCSV[frameNumber,], color=color,size=size)
  print(test)
}