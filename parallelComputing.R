myMatrix <- matrix(sample(c(1:50),1000,replace = TRUE))
testH5.file[["fakeData"]] <- myMatrix

# no_cores <- detectCores() - 1
no_cores <- 10
# Initiate cluster
cl <- makeCluster(no_cores)

clusterExport(cl,c("testH5.file","myMatrix"))
clusterEvalQ(cl,library("hdf5r"))
parSapply(cl,c(1:10),function(x) mean(myMatrix[,,x])) # works
parSapply(cl,c(1:10),function(x) mean(testH5.file[["fakeData"]][,,x])) # fails
stopCluster(cl)