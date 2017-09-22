library(normalRegPanelMix)
library(foreach)
library(Rmpi)
library(stargazer)
library(parallel)
set.seed(123456)

Nset <- c(100,500)
Tset <- c(2,10)

anset <- c(0.05,0.1,0.15,0.3,0.4)

alphaset <- list(c(0.25, 0.25, 0.25, 0.25))
muset 		<- list(c(-4,-1,1,4), c(-5,-1,1,5), c(-6,-2,2,6), c(-6,-1,2,5), c(-5,0,2,4), c(-6,0,2,4))
sigmaset <- list(c(1, 1, 1, 1), c(1, 0.75, 0.5, 0.25))



#The parameters that are fixed
M <- 4 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X
gamma <- matrix(0)
beta <- matrix(0)

#Get the misclassfication
GetMisclTerm <- function(phi) {
  
  m <- phi$M 
  
  if (m == 2) 
  {
    omega.12  <- omega.12(phi)
    return (log(omega.12 /(0.5-omega.12)))
  }
  
  if (m == 3) # ln(omega_12 omega_23 / (0.5-omega_12)(0.5-omega_23))
  {
    omega.123 <- omega.123(phi)
    omega.12 <- omega.123[1]
    omega.23 <- omega.123[2]
    return (log(omega.12 * omega.23 / ((0.5-omega.12)*(0.5-omega.23))))
  }
  omega.1234 <- omega.1234(phi)
  omega.12 <- omega.1234[1]
  omega.23 <- omega.1234[2]
  omega.34 <- omega.1234[3]
  # (m == 4) # ln(omega_12 omega_23 omega_34 / (0.5-omega_12)(0.5-omega_23)(0.5-omega_34))
  return (log(omega.12 * omega.23 * omega.34 / 
                ((0.5-omega.12)*(0.5-omega.23)*(0.5-omega.34))))
  
}

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



PerformEMtest <- function (data, an, m = 2, parallel) {
  # workers might need information
  library(normalRegPanelMix)# workers might need information
  # print(data)
  out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
  out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=an,update.alpha = 1,parallel = FALSE)
  return(2 * max(out.h1$penloglik - out.h0$loglik))
}



MPIgetEstimate <- function(Data,phi,nrep,an,m){
  lr.crit <- matrix(0.0,nr=nrep,ncol=3)
  lr.estimate <- matrix(0.0,nr=nrep,ncol=1)
  lr.size <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  parallel=FALSE
  # cl <- makeCluster(7)
  # for (k in 1:nrep){
  #   data <- Data[,k]
  #   out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
  #   out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=an,update.alpha = 1,parallel = TRUE,cl=cl)
  #   lr.estimate[k] <- 2 * max(out.h1$penloglik - out.h0$loglik)
  # }
  # stopCluster(cl)
  
  
  
  ldata <- lapply(seq_len(ncol(Data)), function(i) Data[,i])
  lr.estimate <- cbind(mpi.applyLB(ldata, PerformEMtest, an = an, m=M,
                                   parallel = parallel))
  data <- ldata[[1]]
  
  out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=m,vcov.method = "none")
  crit <- try(regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z, parallel = FALSE,nrep=1000)$crit)
  if (class(crit) == "try-error"){
    crit <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, nbtsp = 199 ,parallel = FALSE)$crit
  } 
  for ( k in 1:nrep){
    lr.crit[k,] <- crit
    lr.size[k,] <- 1 * (lr.estimate[k,] > lr.crit[k,2])
  }
  
  return(list(est = lr.estimate , crit = lr.crit,nominal.size = apply(lr.size,2,mean)))
}

#GeneratePhiDataPairs
count <- 0
nrep <- 500
phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(sigmaset)
regression.data <- matrix(0,nr=(nset*length(anset)),nc=5)

# Rmpi setup 
print("collecting workers..")
mpi.spawn.Rslaves()
mpi.setup.rngstream()
mpi.bcast.Robj2slave(PerformEMtest, all=TRUE)
print("workers loaded.")

# ====== BEGIN EXPERIMENT ======
## 1. Initialization
# Case when m = 3
for (N in Nset){
  for (T in Tset){
    for (mu in muset){
      for (alpha in alphaset){
        for (sigma in sigmaset){
          for (an in anset){
            t <- Sys.time()
            phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                       beta = beta, N = N, T = T, M = M, p = p, q = q)
            
            count <- count + 1
            
            phi.data.pair <- GenerateSample(phi,nrep)
            Data = phi.data.pair$Data
            phi = phi.data.pair$phi
            
            # phi.data[[count]] <- phi.data.pair
            
            result <- MPIgetEstimate(Data,phi,nrep,an,m=M)
            
            omega <- GetMisclTerm(phi)
            print(result$nominal.size)
            print(Sys.time() - t)
            print(paste(count,"/",(nset*length(anset))) )
            
            regression.data[count, ] <-
              cbind(result$nominal.size,an, phi$N,phi$T,omega.term)
            
            
          }
        }
      }
    }
  }
}


colnames(regression.data) <- c("nom.size", "an" ,"N","T", "omega")
# Begin estimation

stargazer(regression.data)
#The colnames are nom.size, an, T, N
write.csv(regression.data,file="/work/haoyu2/orca/penaltyTestM4.csv",row.names=FALSE)

fit.data <- list(y = log(
  regression.data$nom.size/ (0.11 - regression.data$nom.size) ) ,
  x1 = 1 /  regression.data$T, 
  x2 = 1 /  regression.data$N, 
  x3 = log( regression.data$an/ (1 - regression.data$an )),
  x4 = regression.data$omega) 

fit.data$y <- replace(fit.data$y,fit.data$y == -Inf,-5)
an.model <- lm(y ~ x1 + x2 + x3 + x4 , data=fit.data,na.action = na.omit)
print(summary(an.model))
mpi.close.Rslaves()
mpi.quit()
