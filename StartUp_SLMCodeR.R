# Basic set up

# install.packages(c("xml2","hdf5r"))
library("xml2")
library('hdf5r')
library('nat')

LSMCodeRConfig<-list()
LSMCodeRConfig$srcdir<-normalizePath(dirname(attr(body(function() {}),'srcfile')$filename))
LSMCodeRConfig$maindir<-dirname(LSMCodeRConfig$srcdir)

source(file.path(LSMCodeRConfig$maindir,"calciumImagingFunctions.R"))


source("")