if (!exists("LSMCodeRConfig")) stop("Please start up the project...")

# TODO: Read in all of the physiology files names from disk and compare with
# stimParListForParallel on disk and/or in memory. If there are any missing
# go ahead and load those.
stimParListRDS <- file.path(LSMCodeRConfig$srcdir,"objects","stimParListForParallel.RDS")

myFiles <- dir(imageDir,patt="[A-Z]{4}-[[:graph:]]",full=T,rec=TRUE)


for (myFile in myFiles) {
  # foreach(myFile=physiologyFilesSP, .packages = "hdf5r") %dopar% {
  matFileCode <- substring(basename(myFile),1,4)
  print(paste("Starting to parse",matFileCode))
  if (!matFileCode%in%names(stimulusParametersList)){
    # get the transition frames for the mat file
    # parseImageForGreenFlashAndLSMandStimLogs.R uses **myFile** from the global environment
    source(file.path(LSMCodeRConfig$scriptsUtility,"parseImageForGreenFlashAndLSMandStimLogs.R"))
    print("Passed parseImageForGreenFlash...")
    print("Setting up the trial information")
    # animal[[basename(myFile)]] <- makeTrial(myFile)
    stimulusParametersList[[matFileCode]] <- makeTrial(myFile)
    file.h5$close()
    imageDataSlice$close()
    
  } else {
    print(paste(matFileCode,"already exists in stimulusParametersList. Skipping."))
  }
}

# saving the list inside the loop doesn't make a whole lot of sense...
saveRDS(stimulusParametersList,file=file.path(LSMCodeRConfig$srcdir,"objects","stimParListForParallel.RDS"),compress = TRUE)
