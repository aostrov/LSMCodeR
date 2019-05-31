# Basic set up
# require("xml2")
# require('hdf5r')
# require('nat')
# require("magick")
# require('ggplot2')
# require('squash')
# require('digest')
# require(doParallel)
# require(doSNOW)
# require("viridis")
# require('signal') for sgolay filter

packages <- c("viridis"=require("viridis"),
              "hdf5r"=require('hdf5r'),
              "nat"=require('nat'),
              "magick"=require('magick'),
              "ggplot2"=require('ggplot2'),
              "reshape2"=require('reshape2'),
              "signal"=require('signal'))

# this should be a function that can call itself recursively
for (package in 1:length(packages)) {
  if (!packages[package]) {
    prompt <- readline(prompt=paste("Load required package",names(packages[package]) ,"? [y/n]: "))
    switch (prompt,
      "y" = install.packages(names(packages[package])),
      "Y" = print("WHOOPIE"),
      "n" = print("Things will fail"),
      "N" = print("Capital things will fail"),
      print("choose y or n...")
    )
  }
}

# there must be a better way of doing this these days!
LSMCodeRConfig<-list()
LSMCodeRConfig$srcdir<-normalizePath(dirname(attr(body(function() {}),'srcfile')$filename))
LSMCodeRConfig$maindir<-dirname(LSMCodeRConfig$srcdir)
LSMCodeRConfig$logDir <- file.path(LSMCodeRConfig$srcdir,"logs")
LSMCodeRConfig$protocolDir <- file.path(LSMCodeRConfig$srcdir,"protocolCSVs")

source(file.path(LSMCodeRConfig$srcdir,"calciumImagingFunctions.R"))

protocolList <- list(
  sabineProtocolSimple = list(
    presentationMatrix = matrix(rep(c(3,4,5,6),5),nrow = 5,byrow=TRUE),
    stimulationSections = read.csv(file.path(LSMCodeRConfig$protocolDir,"sabineProtocolSimple")),
    framesSkipped = 30000,
    background.start = 700,
    background.end = 900,
    signal.start = 900,
    signal.end = 1500
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



##################
## magic numbers##
##################
offset <- 399
pixelSize <- 0.8
roiEdgeLength <- 26

######################
## On disk metadata ##
######################
tectumROIs <- read.csv(file.path(LSMCodeRConfig$srcdir,"stuff","tectumROI.csv"))
fishGenos <- read.csv(file.path(LSMCodeRConfig$srcdir,"models","DSLM-fish.csv"))
