library(normalRegPanelMix)
library(foreach)
library(Rmpi)
library(stargazer)
library(parallel)

#Generate Data
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 1 #Number of X

set.seed(123456)
Nset <- c(100,500)
Tset <- c(2,5,10)
alphaset <- list(c(0.5,0.5),c(0.2,0.8))
muset <- list(c(-1,1),c(-0.5,0.5))
sigma <- c(0.8,1.2)
beta <- c(0,0)

GenerateSample <- function(phi,nrep){ 
  p = phi$p
  q = phi$q
  N = phi$N
  T = phi$T
  M = phi$M
  alpha = phi$alpha 
  mu = phi$mu
  gamma = phi$gamma
  beta = phi$beta
  Data <- replicate(nrep,generateData(alpha,mu,sigma,gamma,beta,N,T,M,p,q))
  return(list(phi=phi,Data=Data))
}



PerformEMtest <- function (data, an, m = 2, z = NULL, parallel) {
  # workers might need information
  library(normalRegPanelMix)# workers might need information
  # print(data)
  out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
  out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = data$Z,parlist=out.h0$parlist,an=an,update.alpha = 1,parallel = FALSE)
  return(2 * max(out.h1$penloglik - out.h0$loglik))
}



MPIgetEstimate <- function(Data,phi,nrep,an,m){
  lr.crit <- matrix(0.0,nr=nrep,ncol=3)
  lr.estimate <- matrix(0.0,nr=nrep,ncol=1)
  lr.size <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  ldata <- lapply(seq_len(ncol(Data)), function(i) Data[,i])
  parallel = FALSE
  lr.estimate <- cbind(mpi.applyLB(ldata, PerformEMtest, an = an, m=M, z = NULL,
                                   parallel = parallel))
  data <- ldata[[1]]
  out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
  crit <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z,cl=cl , parallel = FALSE)$crit
  for ( k in 1:nrep){
    lr.crit[k,] <- crit
    lr.size[k,] <- 1 * (lr.estimate[k,] > lr.crit[k,2])
  }
  
  return(list(est = lr.estimate , crit = lr.crit,nominal.size = apply(lr.size,2,mean)))
}


# Rmpi setup 
print("collecting workers..")
# mpi.spawn.Rslaves()
mpi.setup.rngstream()
mpi.bcast.Robj2slave(PerformEMtest, all=TRUE)
print("workers loaded.")


#GeneratePhiDataPairs
count <- 0
nrep <- 5
phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset)

NTset <- expand.grid(Nset,Tset)
Parset <- expand.grid(muset,alphaset)
nNT <- dim(NTset)[1]
nPar <- dim(Parset)[1]
result.f <- matrix(0,nr=(nNT),nc=nPar)
rownames(result.f) <- apply(NTset,1,paste,collapse = ",")
colnames(result.f) <- apply(Parset,1,paste,collapse = ",")

for (r in 1:nNT){
  N <-  NTset[r,1]
  T <-  NTset[r,2]
  print(paste(r,"/",nNT))
  count <- 0
  for (mu in muset){
    for (alpha in alphaset){
      
      t <- Sys.time()
      phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                 beta = beta, N = N, T = T, M = M, p = p, q = q)
      
      phi.data.pair <- GenerateSample(phi,nrep)
      count <- count + 1
      Data = phi.data.pair$Data
      phi = phi.data.pair$phi
      # phi.data[[count]] <- phi.data.pair
      an <- 0.5 *anFormula(phi,M,N,T,q=1) #The an function according the the empirical regression
      # an <- 0.0005
      result <- MPIgetEstimate(Data,phi,nrep,an,m=M)
      
      # print(result$nominal.size)
      result.f[r, count] <- result$nominal.size
      print(Sys.time() - t)
      
    }
    
  }
  
}


write.csv(result.f,file="sizeTestM2Regressor.csv")

mpi.bcast.cmd(mpi.quit(),caller.execute=TRUE)
mpi.close.Rslaves()
mpi.quit()