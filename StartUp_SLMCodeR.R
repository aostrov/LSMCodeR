# Basic set up

# install.packages(c("xml2","hdf5r"))
library("xml2")
library('hdf5r')
library('nat')

LSMCodeRConfig<-list()
LSMCodeRConfig$srcdir<-normalizePath(dirname(attr(body(function() {}),'srcfile')$filename))
LSMCodeRConfig$maindir<-dirname(LSMCodeRConfig$srcdir)
LSMCodeRConfig$logDir <- file.path(LSMCodeRConfig$srcdir,"logs")

source(file.path(LSMCodeRConfig$srcdir,"calciumImagingFunctions.R"))
