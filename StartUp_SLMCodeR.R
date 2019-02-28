# Basic set up
# packages=c('hdf5r','nat','magick','ggplot2')
# sapply(packages,require)
# install.packages(c("xml2","hdf5r"))
require("xml2")
require('hdf5r')
require('nat')
require("magick")
require('ggplot2')
require('squash')
require('digest')

LSMCodeRConfig<-list()
LSMCodeRConfig$srcdir<-normalizePath(dirname(attr(body(function() {}),'srcfile')$filename))
LSMCodeRConfig$maindir<-dirname(LSMCodeRConfig$srcdir)
LSMCodeRConfig$logDir <- file.path(LSMCodeRConfig$srcdir,"logs")
LSMCodeRConfig$protocolDir <- file.path(LSMCodeRConfig$srcdir,"protocolCSVs")

source(file.path(LSMCodeRConfig$srcdir,"calciumImagingFunctions.R"))

pixelOffset=399
protocolList <- list(
  sabineProtocolSimple = list(
    presentationMatrix = matrix(rep(c(3,4,5,6),5),nrow = 5,byrow=TRUE),
    stimulationSections = read.csv(file.path(LSMCodeRConfig$protocolDir,"sabineProtocolSimple")),
    framesSkipped = 30000
  )
)

if (dir.exists("F:\\Imaging\\GCaMP7_tests\\20181204-g7")) {
  print("Setting imageDir to F:/Imaging/GCaMP7_tests/20181204-g7")
  imageDir <- "F:/Imaging/GCaMP7_tests/20181204-g7"
} else if (dir.exists("/Volumes/TranscendJD/Work/")) {
  print("Setting imageDir to /Volumes/TranscendJD/Work/")
  imageDir <- "/Volumes/TranscendJD/Work/"
} else {
  print("Please set the variable 'imageDir' to something sensible.")
}