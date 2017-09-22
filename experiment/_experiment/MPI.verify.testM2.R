library(normalRegPanelMix)
library(doParallel)
library(Rmpi)
#Generate Data
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X


Nset <- c(100,500)
Tset <- c(2,5,10)
alphaset <- list(c(0.5,0.5),c(0.2,0.8))
muset <- list(c(-1,1),c(-0.5,0.5))
sigma <- c(0.8,1.2)


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


PerformEMtest <- function (data, an, m = 2, z = NULL, parallel) {
  library(doParallel) # workers might need information
  library(normalRegPanelMix)# workers might need information
  
  out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
  out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=an,update.alpha = 1,parallel = FALSE)
  return(2 * max(out.h1$penloglik - out.h0$loglik))
}



MPIgetEstimate <- function(Data,phi,nrep,an,m){
  lr.crit <- matrix(0.0,nr=nrep,ncol=3)
  lr.estimate <- matrix(0.0,nr=nrep,ncol=1)
  lr.size <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  ptm <- proc.time()[1]
  parallel=FALSE
  ldata <- lapply(seq_len(ncol(Data)), function(i) Data[,i])
  lr.estimate <- cbind(mpi.applyLB(ldata, PerformEMtest, an = an, m=M, z = NULL,
                                   parallel = parallel))
  data <- ldata[[1]]
  out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
  crit <- regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = FALSE,nrep=1000)$crit
  for ( k in 1:nrep){
    lr.crit[k,] <- crit
    lr.size[k,] <- 1 * (lr.estimate[k,] > lr.crit[k,2])
  }

  print(proc.time()[1] - ptm)
  return(list(est = lr.estimate , crit = lr.crit,nominal.size = apply(lr.size,2,mean)))
}

#GeneratePhiDataPairs
count <- 0
nrep <- 200
phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset)
regression.data <- matrix(0,nr=(nset),nc=3)


#Initiate Rmpi slaves
mpi.spawn.Rslaves()
#Broadcast library to slaves

for (N in Nset){
  for (T in Tset){
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
        an <- anFormula(phi,M,N,T) * 0.5
        print(an)
        
        result <- MPIgetEstimate(Data,phi,nrep,an,m=M)
        
        # print(result$nominal.size)
        regression.data[count, ] <- 
          cbind(result$nominal.size,phi$T,phi$N)
        # write.csv(result$lr.estimate,file=paste("estimate.result", count))
        
        print(Sys.time() - t)
      }
    }
  }
}

mpi.close.Rslaves()
mpi.quit()
write.csv(regression.data,file="estimate.result")

