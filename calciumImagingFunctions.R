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

# for presentationList2 it will return the index of which list 
# element contains a particular stimulus and block pair
# This is an incredibly stupid function at the moment
getIndexOfStimulusBlockPair <- function(stimulus,block){
  index <- which(
    unlist(
      lapply(
        presentationList2, function(x) {
          x$block==block & x$stimulus==stimulus
        }
      )
    )
  )
  if (length(index) < 1) stop("You done fucked up")
  return(index)
}

# a little something to make my life easier when I want to either
# loop through all of the sitmuli and make either dFF or SNR
# or just reprocess a single one.
processSingleStimulus <- function(myList,stimulus=1,block=1,
                                  outputType=c("snr","raw","dff"),writeNRRD=FALSE,
                                  downSampleImage=FALSE, setFloor = FALSE, ...) {
  # I'll probably want some kind of checking here for existing
  # files, etc
  
  outputType = match.arg(outputType)
  i <- getIndexOfStimulusBlockPair(stimulus,block)
  # print(i)
  
  rangeOfImages <- seq(from=myList[[i]]$start,
                       to=myList[[i]]$end,
                       by=myList[[i]]$timeResampled)
  
  myImage <- image[rangeOfImages,,,,]
  
  if (downSampleImage){
    myImage <- apply(
    image[rangeOfImages,,,,],
    1,
    function(x) resizeImage(x,myList[[i]]$resize[1],
                            myList[[i]]$resize[2]))
  
  dim(myImage)<-c(myList[[i]]$resize[1],
                           myList[[i]]$resize[2],
                           length(rangeOfImages))
  }

  myImage <- returnOffsetedImage(myImage,offsetValue=pixelOffset, setFloor = setFloor, ...)

  if (outputType == 'dff'){
    myImage <- makeDFF(
      myImage,
      xyzDimOrder = c(1,2,3),
      backgroundSlices=myList[[i]]$backgroundSlices
    )
    
  } else if (outputType=="snr") {
    myImage <- makeSNRByPixel(myImage,
                              backgroundSlices=myList[[i]]$backgroundSlices
    )
  }
  
  if (writeNRRD){
    outfile <- file.path(myList[[i]]$outDir,paste(myList[[i]]$fileBaseName,"_",outputType,".nrrd",sep=""))
    # print(paste("Writing to disk: ",outfile,sep=""))
    write.nrrd(myImage,file=outfile) #,dtype="short"
  }
  
  invisible(myImage)
  
  
}

processSingleStimulus.lapply <- function(myList,
                                  outputType=c("snr","raw","dff"),writeNRRD=FALSE,
                                  downSampleImage=FALSE, setFloor = FALSE, image,
                                  downSampleInTime=1,z_plane=1,
                                  resizeFactor=1, outDir=file.path(LSMCodeRConfig$srcdir,"images"),
                                  offsetValue=399,
                                  ...) {
  # I'll probably want some kind of checking here for existing
  # files, etc
  
  cat("starting",sep="\n")
  outputType = match.arg(outputType)
  image.dims=image$dims
  names(image.dims) <- c('t','c','z','y','x')

  rangeOfImages <- seq(from=myList$start,
                       to=myList$end,
                       by=downSampleInTime)
  
  if (downSampleImage){
    myImage <- apply(
      image[rangeOfImages,,z_plane,,],
      1,
      function(x) resizeImage(x,
                              image.dims[['x']]/resizeFactor,
                              image.dims[['y']]/resizeFactor)
      )
    
    dim(myImage)<-c(image.dims[['x']]/resizeFactor,
                    image.dims[['y']]/resizeFactor,
                    length(rangeOfImages))
  } else {
    myImage <- aperm(image[rangeOfImages,,z_plane,,],c(3,2,1))
  }
  
  myImage <- returnOffsetedImage(myImage,offsetValue=offsetValue, setFloor = setFloor, ...)
  
  if (outputType == 'dff'){
    myImage <- makeDFF(
      myImage,
      xyzDimOrder = c(1,2,3),
      backgroundSlices=c(0:length(myList$backgroundSlices))
    )
    
  } else if (outputType=="snr") {
    myImage <- makeSNRByPixel(myImage,
                              backgroundSlices=c(0:length(myList$backgroundSlices))
    )
  }
  
  if (writeNRRD){
    outfile <- file.path(outDir,paste(myList$matFile,"_stimulus-",myList$stimulus,
                                      "_block-",myList$block,"_",
                                      outputType,
                                      "_zplane-",z_plane,".nrrd",
                                      sep=""))
    # print(paste("Writing to disk: ",outfile,sep=""))
    write.nrrd(myImage,file=outfile) #,dtype="short"
  }
  
  invisible(myImage)
  
  
}

# takes a SNR file and turns it into subROIs that show the maximum change SNR over background
getSNRforSubROIs <- function(imageData,numSubunits,x,y,w,h,
                             f0Window=c(75:85),responseWindow=c(90:120),
                             writeNRRD=FALSE,dryRun=FALSE){
  SNR_f0 <- c()
  SNR_max <- c()
  x.array <- c()
  y.array <- c()
  count=0
  
  for (yy in (seq(numSubunits)-1)){
    for (xx in (seq(numSubunits)-1)){
      cat(".")
      if (dryRun){
        print(paste("XX:",xx))
        print(paste("YY:",yy))
        print(paste("xxx[(x+(",xx,"*w)/",numSubunits,"):(x+(",xx+1,"*w)/"
                    ,numSubunits,"),(y+(",yy,"*h)/",
                    numSubunits,"):(y+(",yy+1,"*h)/",
                    numSubunits,"),]",sep=""
        )
        )
      } else {
        if (writeNRRD){
          write.nrrd(imageData[(x+(xx*w)/numSubunits):(x+((xx+1)*w)/numSubunits),
                               (y+(yy*h)/numSubunits):(y+((yy+1)*h)/numSubunits),],
                     file=paste("C:/Users/Aaron/Desktop/temp/numSubunits_",
                                numSubunits,"_",xx+yy+count,".nrrd",sep=""),
                     dtype="short")
        }
        
        SNR_f0 = c(SNR_f0,mean(imageData[(x+(xx*w)/numSubunits):(x+((xx+1)*w)/numSubunits),
                                         (y+(yy*h)/numSubunits):(y+((yy+1)*h)/numSubunits),
                                         f0Window]))
        
        SNR_max = c(SNR_max,max(apply(
          imageData[(x+(xx*w)/numSubunits):(x+((xx+1)*w)/numSubunits),
                    (y+(yy*h)/numSubunits):(y+((yy+1)*h)/numSubunits),responseWindow],
          c(1,2),mean)))
        x.array <- c(x.array,xx)
        y.array <- c(y.array,yy)
        
        
        
      }
    }
    count=count+numSubunits
  }
  # Return block
  if (!dryRun){
    df <- data.frame(SNR_f0,SNR_max)
    df$dSNR <- df$SNR_max - df$SNR_f0
    df$X <- x.array
    df$Y <- y.array
    return(df)
  }
  
}

# Given a rectangular area, return a list of coordinates of
# sub-rectangles that either divide up the full area based
# on the number of units requested, or based on edge length
# of the ROIs
# If the number of subunits is selected, the subunits inherit
# the aspect ratio of the outer box.
# If the edge length is specified, the ROIs are square.
# If the outer box is not perfectly divisble by the number or
# size of subunits desired, the outer box is adjusted to the 
# next smallest size to evenly accomodate everything.
# TODO: Right now the subROIs overlap by 1 pixel
getROIs <- function(numSubunits=NULL,roiEdgeLength=NULL,
                    x=1,y=1,w=300,h=400){
  if ( (!is.null(numSubunits) && !is.null(roiEdgeLength)) || (is.null(numSubunits) && is.null(roiEdgeLength)) ) {
    stop("Choose either a number of ROIs or the edge length of square ROIs")
  } 
  divisor <- ifelse(is.null(numSubunits),roiEdgeLength,numSubunits)
  w <- roundToNearestDivisibleInteger(w,divisor)
  h <- roundToNearestDivisibleInteger(h,divisor)
  if ( w%%divisor!=0 | h%%divisor!=0 ) stop("AAAHHHH!")
  count=1
  roiList <- list()
  for ( yy in seq( ( ifelse(is.null(numSubunits),h/roiEdgeLength,numSubunits) ) )-1 ) {
    for ( xx in seq( ( ifelse(is.null(numSubunits),w/roiEdgeLength,numSubunits) ) )-1 ){
      xRange <- c(
        (x+(xx* ifelse(is.null(numSubunits),roiEdgeLength,w/divisor) )):
          (x+((xx+1)* ifelse(is.null(numSubunits),roiEdgeLength,w/divisor) )))
      yRange <- c(
        (y+(yy* ifelse(is.null(numSubunits),roiEdgeLength,h/divisor) )):
          (y+((yy+1) * ifelse(is.null(numSubunits),roiEdgeLength,h/divisor) ))
      )
      roiList[[paste("ROI",count,sep="_")]] <- list(xRange=xRange,yRange=yRange,xPosition=xx,yPosition=yy)
      count=count+1
    }
  }
  attrVector <- c(x,y,w,h,divisor)
  names(attrVector) <- c("X_origin","Y_origin","width","height",ifelse(is.null(numSubunits),"roiEdgeLength","numSubunits"))
  attr(roiList, "ROI_info" ) = attrVector
  return(roiList)
}


getROIsRawDataFromHDF5.lapply <- function(roiListElement,hdf5Image.mat,frame.start,frame.end,z=1,offset=399) {
  tempY <- roiListElement$yRange
  tempX <- roiListElement$xRange
  # get hdf5 image and figure out which slices I need
  # I would get this information from presentationList2
  # or whatever ends up replacing it
  tempImage <- hdf5Image.mat[c(frame.start:frame.end),,z,tempY,tempX]
  tempImage <- aperm(tempImage,c(3,2,1)) - offset
  # print(dim(tempImage))
  imageAttributes <- c(frame.start,frame.end,offset)
  names(imageAttributes) <- c("frame.start","frame.end","pixel_offset")
  attr(tempImage, "imageAttributes") <- imageAttributes
  roiAttributes <- c(xPosition=roiListElement$xPosition,yPosition=roiListElement$yPosition)
  attr(tempImage,"roiAttributes") <- roiAttributes
  return(tempImage)
}

getUsefulStatisticsByROI <- function(rawDataByROI,roiList,analysisWindow=c(900:1500),backgroundWindow=c(750:850)) {
  usefulStatisticsByROI.DF <- data.frame(
    # raw but offset correct background
    background.mean=sapply(rawDataByROI,function(x) mean(x[,,backgroundWindow])),
    background.max=sapply(rawDataByROI,function(x) max(x[,,backgroundWindow])),
    background.sd=sapply(rawDataByROI,function(x) sd(x[,,backgroundWindow])),
    # raw but offset correct data
    raw.mean=sapply(rawDataByROI,function(x) mean(x[,,analysisWindow])),
    raw.max=sapply(rawDataByROI,function(x) max(x[,,analysisWindow])),
    raw.sd=sapply(rawDataByROI,function(x) sd(x[,,analysisWindow])),
    # offset correct SnR data
    snr.mean=sapply(rawDataByROI,function(x) {
      mean(apply(makeSNRByPixel(x,backgroundSlices=backgroundWindow)[,,analysisWindow],3,mean))
    }),
    snr.max=sapply(rawDataByROI,function(x) {
      max(apply(makeSNRByPixel(x,backgroundSlices=backgroundWindow)[,,analysisWindow],3,mean))
    }),
    snr.sd=sapply(rawDataByROI,function(x) {
      sd(apply(makeSNRByPixel(x,backgroundSlices=backgroundWindow)[,,analysisWindow],3,mean))
    }),
    # dff on offset data, with a small fudge for when NaNs appear after 0/0
    dff.mean=sapply(rawDataByROI, function(x){
      dff <- mean(apply(makeDFF(x,backgroundSlices=backgroundWindow,xyzDimOrder=c(1,2,3))[,,analysisWindow],3,mean))
      dff[is.nan(dff)]=0
      return(dff)
    }),
    dff.max=sapply(rawDataByROI, function(x){
      dff <- max(apply(makeDFF(x,backgroundSlices=backgroundWindow,xyzDimOrder=c(1,2,3))[,,analysisWindow],3,mean))
      dff[is.nan(dff)]=0
      return(dff)
    }),
    dff.sd=sapply(rawDataByROI, function(x){
      dff <- sd(apply(makeDFF(x,backgroundSlices=backgroundWindow,xyzDimOrder=c(1,2,3))[,,analysisWindow],3,mean))
      dff[is.nan(dff)]=0
      return(dff)
    }),
    # X and Y positions
    xpos=sapply(rawDataByROI,function(x) attr(x,"roiAttributes")['xPosition']),
    ypos=sapply(rawDataByROI,function(x) attr(x,"roiAttributes")['yPosition']),
    frame.start=sapply(rawDataByROI,function(x) attr(x,"imageAttributes")["frame.start"])
  )
  return(usefulStatisticsByROI.DF)
}



roundToNearestDivisibleInteger <- function(numerator,divisor,roundUp=FALSE){
  if (numerator%%divisor != 0){
    remainder <- numerator%%divisor
    if (roundUp) {
      difference <- divisor - remainder
      divisible <- numerator + difference
    } else {
      divisible <- numerator - remainder
    }
    return(divisible)
  } else {
    return(numerator)
  }
}



makeDFF<-function(imageStack,backgroundSlices=c(50:101),xyzDimOrder=c(3,2,1)){
  # establish what the actual dimensions are of the image in X, Y, and Z
  zdim<-dim(imageStack)[xyzDimOrder[3]]
  xdim<-dim(imageStack)[xyzDimOrder[1]]
  ydim<-dim(imageStack)[xyzDimOrder[2]]

  # get background
  testSliceBackground<-apply(imageStack[,,backgroundSlices],c(xyzDimOrder[1],xyzDimOrder[2]),mean)
  
  # subtract background
  testSliceBackgroundSubtract<-apply(imageStack[,,],3,function(x) x-testSliceBackground)
  dim(testSliceBackgroundSubtract)<-c(xdim,ydim,zdim)
  
  # divide by background
  testSliceDFF<-apply(testSliceBackgroundSubtract, 3, function(x) x/testSliceBackground)
  dim(testSliceDFF)<-c(xdim,ydim,zdim)
  
  return(testSliceDFF)
}

makeDFFwithBaselineSubtraction <- function(imageStack,backgroundSlices=c(50:101),xyzDimOrder=c(3,2,1),...){
  zdim<-dim(imageStack)[xyzDimOrder[3]]
  xdim<-dim(imageStack)[xyzDimOrder[1]]
  ydim<-dim(imageStack)[xyzDimOrder[2]]

  offset <- getBaseline(imageStack,...)
  print(offset)
  imageStack.offset <- (imageStack - offset)
  dim(imageStack.offset)<-c(xdim,ydim,zdim)
  
  imageStack.offset.F0<-apply(imageStack.offset[,,backgroundSlices],c(1,2),mean)
  imageStack.offset.F0.F<-apply(imageStack.offset[,,],3,function(x) x-round(imageStack.offset.F0))
  dim(imageStack.offset.F0.F)<-c(xdim,ydim,zdim)
  
  imageStack.offset.F0.F.dff<-apply(imageStack.offset.F0.F, 3, function(x) x/imageStack.offset.F0)
  dim(imageStack.offset.F0.F.dff)<-c(xdim,ydim,zdim)
  
  return(imageStack.offset.F0.F.dff)
}

makeSNRByPixel<-function(imageStack,backgroundSlices=c(75:85),xyzDimOrder=c(1,2,3)){
  zdim<-dim(imageStack)[xyzDimOrder[3]]
  xdim<-dim(imageStack)[xyzDimOrder[1]]
  ydim<-dim(imageStack)[xyzDimOrder[2]]
  
  # backgroundSlices <- c(0:length(backgroundSlices))
  noise <- apply(imageStack[,,backgroundSlices],c(xyzDimOrder[1],xyzDimOrder[2]),calcSNR.noise)
  snr <- apply(imageStack,xyzDimOrder[3],function(x) {x/noise} )
  dim(snr) <- c(xdim,ydim,zdim)
  return(snr)
}

makeSNRByPixel.lapply<-function(rois,imageStack,backgroundSlices=c(75:85),xyzDimOrder=c(1,2,3),asIfDFF=FALSE){
  zdim<-dim(imageStack)[xyzDimOrder[3]]
  xdim<-dim(imageStack)[xyzDimOrder[1]]
  ydim<-dim(imageStack)[xyzDimOrder[2]]
  
  
  noise <- apply(imageStack[rois$xRange,rois$yRange,
                            backgroundSlices]
                 ,c(xyzDimOrder[1],xyzDimOrder[2]),calcSNR.noise)
  snr <- apply(imageStack[rois$xRange,rois$yRange,],xyzDimOrder[3],function(x) {x/noise} )
  dim(snr) <- c(length(rois$xRange),length(rois$yRange),zdim)
  # if (asIfDFF) {
  #   snr <- apply(snr,)
  # }
  return(snr)
}

calcSNR.noise <- function(vector) {
  noise <- ((sd(diff(vector)))/sqrt(2))
  return(noise)
}

makeTrial <- function(matFile,stimProtocol="sabineProtocolSimple",analysisWindow=300) {
  trials <- list()
  count=1
  for (block in 1:nrow(protocolList[[stimProtocol]]$presentationMatrix)) {
    for (stimulus in 1:ncol(protocolList[[stimProtocol]]$presentationMatrix)) {
      for (plane in seq(imageDataSlice.dims[['z']])) {
        tmpdf <- subset(
          protocolList[[stimProtocol]]$stimulationSections,
          section==protocolList[[stimProtocol]]$presentationMatrix[block,stimulus]
        )
        description <- subset(
          protocolList[[stimProtocol]]$stimulationSections,
          section==protocolList[[stimProtocol]]$presentationMatrix[block,stimulus] & 
            ( description!="background" & description!="settle" )
        )$description
        
        backgroundLengthInMilliseconds <- tmpdf[tmpdf$description=="background","time"] # in ms
        
        if ( imageDataSlice.dims[['z']] > 1 ) {
          start <- slice.transitions[count,"time"] # in frames
        } else {
          start <- lsm.transition.frames[count]
        }
        
        backgroundSlices <- c( 
          start : 
            (start + (backgroundLengthInMilliseconds * 0.1) / imageDataSlice.dims[['z']] ) 
        )
        
        stimulusPeriod <- sum(tmpdf$time) / imageDataSlice.dims[['z']] # in ms
        
        end <- start + 
          (stimulusPeriod * 0.1) + 
          ( analysisWindow / imageDataSlice.dims[['z']] )
        
      }
      
      trials[[paste(basename(matFile),stimulus,block,sep=".")]] <- list(
        "matFile"=basename(matFile),
        "stimulusProtocol"=stimProtocol,
        "block"=block,
        "stimulus"=stimulus,
        "plane"=plane, #at the moment plane is basically ignored, and will always 
                       # reflect the total number of planes without any indication of how each plane might differ.
        "start"=start,
        "end"=end,
        "stimulusPeriodPerSlice"=stimulusPeriod,
        "backgroundSlices"=backgroundSlices,
        "stimulusDescription"=description,
        "creationDate"=date()
      )
      
      count=count+1
    }
    
  }
  attr(trials,"imageDimensions") <- imageDataSlice.dims
  return(trials)
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

returnOffsetedImage <- function(imageData,xyzDimOrder=c(1,2,3),offsetValue,setFloor=FALSE,floor=0) {
  # zdim<-dim(imageData)[xyzDimOrder[3]]
  # xdim<-dim(imageData)[xyzDimOrder[1]]
  # ydim<-dim(imageData)[xyzDimOrder[2]]
  
  imageData <- imageData - offsetValue

  if (setFloor) {
    imageData[imageData < 0] = floor
  }
  return(imageData)
}

getBaseline <- function(imageStack,rangeOfFrames=NULL,xyzDimOrder=c(3,2,1),x=1,width=10,y=1,height=10,dimension=1) {
  zdim<-dim(imageStack)[xyzDimOrder[3]]
  xdim<-dim(imageStack)[xyzDimOrder[1]]
  ydim<-dim(imageStack)[xyzDimOrder[2]]

  if (is.null(rangeOfFrames)){
    rangeOfFrames <- c( round((zdim[1] * 0.25)) , round((zdim[1] * 0.75)) )
  } 
  if (!all(rangeOfFrames%in%c(1:zdim[1]))) {
    rangeOfFrames <- c( round((zdim[1] * 0.25)) , round((zdim[1] * 0.75)) )
    warning("\n'rangeOfFrames' is outside of the bounds of the array.\nSetting rangeOfFrames to the quartiles")
  }
  
  return(floor(mean(imageStack[rangeOfFrames,(y:(y+height)),(x:(x+width))])))
}

writeNrrdForROISelection <- function(matFile,outPath,testTime=100,numZSlices=20,bitDepth="short"){
  file.h5 <- H5File$new(matFile, 
                        mode = "r")
  imageDataSlice<-file.h5[["imagedata"]]
  imageDataSlice.dims <- imageDataSlice$dims
  names(imageDataSlice.dims) <- c('t','c','z','y','x')
  animal <- substr(basename(matFile),1,4)
  count=1
  for (i in 1:imageDataSlice.dims[['z']]){
    write.nrrd(aperm(imageDataSlice[testTime,,i,,],c(2,1)),
               file=file.path(outPath,paste(animal,"_roiSelection-",i,".nrrd",sep="")),
               dtype = bitDepth)
    count=count+1
  }
  imageDataSlice$close()
  file.h5$close_all()
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
# If a multi-frame image is passed as an array form [frames,channels,slices,rows,columns]
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