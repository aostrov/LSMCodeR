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

registerDoParallel(cores = 8)
foreach((i=1:10)) %dopar% mean(myMatrix[,,i])
foreach(i=1:10) %dopar% mean(testH5.file[["fakeData"]][,,i]) # Works!!!!

foreach(myFile=physiologyFilesSP, .packages = "hdf5r") %dopar% {
  file.h5 <- H5File$new(file.path(myFile), mode = "r")
  test.dim <- file.h5[['imagedata']]$dims
  file.h5$close()
  test.dim
}
