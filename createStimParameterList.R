if (!exists("LSMCodeRConfig")) stop("Please start up the project...")

# TODO: Read in all of the physiology files names from disk and compare with
# stimParListForParallel on disk and/or in memory. If there are any missing
# go ahead and load those.
stimParListRDS <- file.path(LSMCodeRConfig$srcdir,"objects","stimParListForParallel.RDS")

tectumROIs <- read.csv(file.path(LSMCodeRConfig$srcdir,"stuff","tectumROI.csv"))
myFiles <- dir(imageDir,patt="[A-Z]{4}-[[:graph:]]",full=T,rec=TRUE)


for (myFile in myFiles) {
  # foreach(myFile=physiologyFilesSP, .packages = "hdf5r") %dopar% {
  matFileCode <- substring(basename(myFile),1,4)
  print(paste("Starting to parse",matFileCode))
  if (T){
    # get the transition frames for the mat file
    # parseImageForGreenFlashAndLSMandStimLogs.R uses **myFile** from the global environment
    source(file.path(LSMCodeRConfig$srcdir,"parseImageForGreenFlashAndLSMandStimLogs.R"))
    print("Passed parseImageForGreenFlash...")
    print("Setting up the trial information")
    # animal[[basename(myFile)]] <- makeTrial(myFile)
    stimulusParametersList[[matFileCode]] <- makeTrial(myFile)
    saveRDS(stimulusParametersList,file=file.path(LSMCodeRConfig$srcdir,"objects","stimParListForParallel.RDS"),compress = TRUE)
    file.h5$close()
    imageDataSlice$close()
    
  }
}
