############
##Anatomy##
############

myFile<-"F:/Imaging/GCaMP7_tests/20181204-g7/20181204-gcamp7F-7d-anatomyPost/20181204-gcamp7F-7d-anatomyPost.mat"
outDir<-"C:/Users/Aaron/Desktop/nrrdFiles/"
outDirSubDir<-"20181204-gcamp7F-7d-anatomyPost"
reorderImageSlices=NULL # or c(84:100,1:83)
roundOutput = TRUE
dtype = 'short'

if (!dir.exists(file.path(outDir,outDirSubDir))){
  dir.create(file.path(outDir,outDirSubDir))
} else {
  print("Dir exists")
}



# open the connection to the .mat file
file.h5 <- H5File$new(myFile, mode = "r")
# save a new variable that contains the image data
# this step can be avoided, and might improve performance on very large datasets
imageData<-file.h5[["imagedata"]]
# file.h5$close()
# data is in the form:
# [stacks,(channels),slices,rows, columns]


write.nrrd(frameAverageForAnatomyStacks(imageData[,,,,],
                                        reorderImageSlices = reorderImageSlices,
                                        roundOutput = roundOutput),
           file.path(outDir,outDirSubDir,"anatomyAverageReordered.nrrd"),
           dtype = dtype)



file.h5$close()