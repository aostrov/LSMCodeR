####################
## Anatomy take 2 ##
####################

anatomyFiles <- dir("/Volumes/home/anatomy",rec=T,patt="[A|a]natomy.mat",full=T)
anatomyFiles <- sample(anatomyFiles,length(anatomyFiles))

for (anatomyFile in anatomyFiles) {
  checkForExistingFile(anatomyFile,writeAnatomyStacksFromHD5,"~/Desktop/tmpNRRD/")
}

# Another single use function, this takes an Anatomy.mat file and makes
# and average image of the slices, re-orders them as needed, and also
# makes a max-intensity z-projection. No frills, no ability to easily
# change paths. Both files are gzip encoded nrrds. Lots of problems with
# this code, but it works for now.
writeAnatomyStacksFromHD5 <- function(matFileWithPath,outdir,outFileBaseName) {
  file.h5 <- H5File$new(matFileWithPath,mode = "r")
  imageData=file.h5[["imagedata"]]
  imageData.dims <- imageData$dims
  names(imageData.dims) <- c("f","c","z","y","x")
  
  reorderImageSlices <- switch (as.character(imageData.dims[['z']]),
                                "50"=c(32:1,50:35),
                                "100"=c(84:100,1:83),
                                NULL
  )
  
  header=list(type="uint16",
              encoding="gzip",
              endian="little",
              dimension=3,
              sizes=c(imageData.dims[['x']],
                      imageData.dims[['y']],
                      imageData.dims[['z']]), # voxel dims
              `space dimensions`=matrix(c(pixelSize,0,0,
                                          0,pixelSize,0,
                                          0,0,(150/imageData.dims[['z']])),
                                        nrow=3,byrow = T), # voxel sizes
              `space units`=rep("microns",3)
  )
  
  anatomyStack <- frameAverageForAnatomyStacks(imageData[,,,,],
                                               reorderImageSlices = reorderImageSlices,
                                               roundOutput = T)
  attr(anatomyStack, "header") <- header
  
  write.nrrd(anatomyStack,
             file.path(outdir,outFileBaseName),
             dtype = "short")
  write.nrrd(intensityProjection(anatomyStack,c(1,2),"max"),
             file.path(outdir,paste("MAX_",outFileBaseName,sep="")),
             dtype = "short")
  file.h5$close()
  imageDataSlice$close()
  
}
  


# A very single use function at the moment, it checks if a file exists,
# or if a lock file based on that file exists. This would in principle
# allow for the parallelization of the creation of the anatomy stacks
# This function automates the creation and deletion of lock files, using
# cat and on.exit(unlink()) respectively
checkForExistingFile <- function(inputFileWithFullPath,cmd,outputDirectory,...) {
  input <- basename(inputFileWithFullPath)
  fishName=substr((basename(inputFileWithFullPath)),1,3)
  geno=sub(" ","_",fishGenos[fishGenos$Fish==fishName,"Full_geno"])
  fishNameFull <- paste(fishName,"_",geno,"-Anatomy.nrrd",sep="")
  # lockfile needs to be generalized
  lockFile <- file.path(outputDirectory,paste(fishNameFull,".lock",sep=""))
  if (file.exists(file.path(outputDirectory,fishNameFull))) return(print("File exists. Skipping."))
  if (file.exists(lockFile)) return("File is being worked on. Skipping.")

  cat("some text here",file=lockFile)
  on.exit(unlink(lockFile))

  cmd(inputFileWithFullPath,outputDirectory,fishNameFull)
}
