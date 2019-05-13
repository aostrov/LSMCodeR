dirByDate <- function(date="2019-03-05",directory=imageDir){
  subsettedDirByDate <- dir(directory,full=T)[
    grepl(date,file.info(dir(directory,full=T))$mtime)
    ]
  return(subsettedDirByDate)
}

recent <- dirByDate()
recentPhysio <- recent[!grepl("Anatomy",recent)]

for (image in recentPhysio) {
  print(image)
  writeNrrdForROISelection(file.path(image,paste(basename(image),"mat",sep=".")),"C:/Users/Aaron/Desktop/nrrdOrder/")
}
