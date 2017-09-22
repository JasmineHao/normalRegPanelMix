library(foreach)
library(parallel)
library(doParallel)
library(snow)
library(Rmpi)

cl <- makeMPIcluster(mpi.universe.size())
registerDoParallel(cl)

t <- Sys.time()
a <- foreach(i=1:24,.combine = rbind) %dopar% {
     i}

print(a)
print(Sys.time() - t)

mpi.close.Rslaves()
mpi.quit()