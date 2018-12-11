# Basic functions for dealing with microscopy data

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
