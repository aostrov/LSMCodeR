# Basic set up
# packages=c('hdf5r','nat','magick','ggplot2')
# sapply(packages,require)
# install.packages(c("xml2","hdf5r"))
require("xml2")
require('hdf5r')
require('nat')
require("magick")
require('ggplot2')

LSMCodeRConfig<-list()
LSMCodeRConfig$srcdir<-normalizePath(dirname(attr(body(function() {}),'srcfile')$filename))
LSMCodeRConfig$maindir<-dirname(LSMCodeRConfig$srcdir)
LSMCodeRConfig$logDir <- file.path(LSMCodeRConfig$srcdir,"logs")

source(file.path(LSMCodeRConfig$srcdir,"calciumImagingFunctions.R"))
