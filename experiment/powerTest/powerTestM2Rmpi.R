library(normalRegPanelMix)
library(doParallel)
library(Rmpi)

set.seed(123456)

#Generate Data
M <- 3 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X


Nset <- c(100,500)
Tset <- c(2,5)
alphaset <- list(c(1/3,1/3,1/3),c(0.25,0.5,0.25))
muset <- list(c(-1,0,1),c(-1.5,0,1.5),c(-1,0,2),c(-0.5,0,1.5))
sigmaset <- list(c(1,1,1),c(0.6,1.2,0.6),c(0.6,0.6,1.2))


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
  # if (q != 0){
  #   parlist <- list('alpha' = alpha,
  #                   'mubeta' = t(cbind(mu,beta)),
  #                   'sigma' = sigma, 'gam' = gamma)
  # }else{
  #   parlist <- list('alpha' = alpha, 'mubeta' = mu, 'sigma' = sigma, 'gam' = gamma)
  # }
  Data <- replicate(nrep,generateData(alpha,mu,sigma,gamma,beta,N,T,M,p,q))
  return(list(phi=phi,Data=Data))
}



PerformEMPowerTest <- function (data , m, z = NULL) {
  # workers might need information
  library(normalRegPanelMix)# workers might need information
  T <- dim(data$Y)[1]
  N <- dim(data$Y)[2]
  out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
  an <- anFormula(out.h0$parlist,m,N,T) * 0.5
  out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=an,update.alpha = 1,parallel = FALSE)
  return(2 * max(out.h1$penloglik - out.h0$loglik))
}


MPIgetEstimatePower <- function(Data,nrep,m){
  lr.crit <- matrix(0.0,nr=nrep,ncol=3)
  lr.estimate <- matrix(0.0,nr=nrep,ncol=1)
  lr.size <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  
  ldata <- lapply(seq_len(ncol(Data)), function(i) Data[,i])
  lr.estimate <- cbind(mpi.applyLB(ldata, PerformEMPowerTest, m=(M-1), z = NULL))
  
  #####################################################
  data <- ldata[[1]]
  out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
  crit <- regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z , parallel = FALSE,nrep=1000)$crit
  for ( k in 1:nrep){
    lr.crit[k,] <- crit
    lr.size[k,] <- 1 * (lr.estimate[k,] > lr.crit[k,2])
  }
  #####################################################
  
  return(list(est = lr.estimate , crit = lr.crit,nominal.size = apply(lr.size,2,mean)))
}

#GeneratePhiDataPairs
count <- 0
nrep <- 1000
phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(sigmaset)
power.data <- matrix(0,nr=(nset),nc=6)


# Rmpi setup 
print("collecting workers..")
mpi.spawn.Rslaves()
mpi.setup.rngstream()
mpi.bcast.Robj2slave(PerformEMtest, all=TRUE)
print("workers loaded.")

for (mu in muset){
  for (sigma in sigmaset){
    for (alpha in alphaset){
      for (N in Nset){
        for (T in Tset){
          count <- count + 1
          
          t <- Sys.time()
          phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                     beta = beta, N = N, T = T, M = M, p = p, q = q)
          
          phi.data.pair <- GenerateSample(phi,nrep)
          Data = phi.data.pair$Data
          phi = phi.data.pair$phi
          
          result <- MPIgetEstimatePower(Data,nrep,M)
          
          
          power.data[count, 1:3] <- 
            cbind(100*(result$nominal.size),phi$T,phi$N)
          power.data[count,4] <- paste(mu,collapse=",")
          power.data[count,5] <- paste(alpha,collapse=",")
          power.data[count,6] <- paste(sigma,collapse=",")
          
          print(Sys.time() - t)
        }
      }
    }
  }
}


write.csv(power.data,file="power.simulate.csv")

mpi.close.Rslaves()
mpi.quit()