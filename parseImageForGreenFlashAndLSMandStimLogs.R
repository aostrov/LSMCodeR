#################
# file handling #
#################

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
# frames channels slices rows  columns
# 33001     1        1    512   700

#################
## Green Flash ##
#################

if (imageDataSlice.dims[['z']] > 1 ) {
  print("Z is bigger than 1, first check in green flash parser")
  flashMean <- mean(imageDataSlice[1:100,,,1:200,1:200])
  flashSD <- sd(imageDataSlice[1:100,,,1:200,1:200])
  
  slicesWithBigNumbers.begin <- c()
  for (i in 1:imageDataSlice.dims[['z']]){
    if (mean(imageDataSlice[1:10,,i,1:200,1:200])>(flashMean+2*flashSD)){
      slicesWithBigNumbers.begin <- c(slicesWithBigNumbers.begin,i)
    }
  }

} else {
  flashmean <- mean(imageDataSlice[1:100,,,10:20,10:20])
  flashsd <- sd(imageDataSlice[1:100,,,10:20,10:20])
  
  slicesWithBigNumbers.begin <- c()
  for (i in 1:100){
    if (mean(imageDataSlice[i,,,,])>(flashmean+2*flashsd)){
      slicesWithBigNumbers.begin <- c(slicesWithBigNumbers.begin,i)
    }
  }
  
  
}

if (is.null(slicesWithBigNumbers.begin)) {
  print(paste(matFileCode, "does not appear to have a green flash :( Setting slicesWithBigNumbers.begin to 0."))
  slicesWithBigNumbers.begin <- 1
}
startOfStimulations <- max(slicesWithBigNumbers.begin) # last slice of starting green flash


endOfFirstGreenFlash <- startOfStimulations + 1
print(paste("end of first green flash:",endOfFirstGreenFlash))


source(file.path(LSMCodeRConfig$srcdir,"stimLogParser.R"))
source(file.path(LSMCodeRConfig$srcdir,"LSMLogParser.R"))
