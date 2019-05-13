dirByDate <- function(date=format(Sys.time(), "%Y-%m-%d"),directory=imageDir){
  subsettedDirByDate <- dir(directory,full=T)[
    grepl(date,file.info(dir(directory,full=T))$mtime)
    ]
  return(subsettedDirByDate)
}

recent <- dirByDate()
recentPhysio <- recent[!grepl("[A|a]natomy",recent)]

for (image in recentPhysio) {
  print(image)
  writeNrrdForROISelection(file.path(image,paste(basename(image),"mat",sep=".")),"C:/Users/Aaron/Desktop/nrrdOrder/")
}
