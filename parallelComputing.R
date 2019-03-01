library("hdf5r")
library("parallel")
install.packages("doParallel")
library(doParallel)

myMatrix <- matrix(sample(c(1:50),1000,replace = TRUE))
dim(myMatrix) <- c(10,10,10)
testH5.file <- H5File$new("~/Desktop/testH5.h5", mode = "a")
testH5.file[["fakeData"]] <- myMatrix

# Parallel package
no_cores <- detectCores() - 1
# no_cores <- 10
# Initiate cluster
cl <- makeCluster(no_cores)

clusterExport(cl,c("testH5.file","myMatrix"))
clusterEvalQ(cl,library("hdf5r"))
parSapply(cl,c(1:10),function(x) mean(myMatrix[,,x])) # works
parSapply(cl,c(1:10),function(x) mean(testH5.file[["fakeData"]][,,x])) # fails
stopCluster(cl)

testH5.file$close_all()

# doParallel package

foreach(i=1:10) %do% mean(myMatrix[,,i])
foreach((i=1:10)) %dopar% mean(myMatrix[,,i])
foreach(i=1:10) %dopar% mean(testH5.file[["fakeData"]][,,i]) # Works!!!!

numCores=3
if (Sys.info()[["sysname"]] == "Darwin") {
  registerDoParallel(cores = numCores)
  cat("Registering parallel execution for MacOS")
} else if (Sys.info()[["sysname"]] == "Windows") {
  registerDoSNOW(makeCluster(numCores,type="SOCK")) 
  getDoParWorkers() 
  cat("Registering parallel execution for Windows")
} else {
  cat("I can't determine the OS!!!")
}

# oddly enough, this seems to work on my mac. I need to see if I have to
# set it up a bit differently on a windows machine :(
physiologyFilesSP <- dir(imageDir,patt="[A-Z]{4}-[[:graph:]]*SP",full=T,rec=TRUE)
foreach(myFile=physiologyFilesSP, .packages = c("hdf5r",'nat')) %dopar% {
  cat(myFile,file="~/Desktop/parallelDebuggingHell.txt",append=T,sep="\n")
  testH5.file <- H5File$new(myFile, mode = "r")
  testH5.file.dims <- testH5.file[["imagedata"]]$dims
  write.nrrd(aperm(testH5.file[["imagedata"]][1:10,,1,100:300,100:300],c(3,2,1)),file=file.path("~/Desktop/",paste(basename(myFile),"nrrd",sep=".")),dtype="short")
  testH5.file$close()
  testH5.file.dims
}
# 