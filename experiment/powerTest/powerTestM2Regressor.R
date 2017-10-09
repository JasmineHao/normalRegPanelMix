library(normalRegPanelMix)
library(doParallel)
#Generate Data
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 1 #Number of X

set.seed(123456)
Nset <- c(100,500)
Tset <- c(2,10)

alphaset <- list(c(1/3,1/3,1/3),c(0.25,0.5,0.25))
muset <- list(c(-1,0,0.5),c(-0.5,0,0.5))
sigma <- c(1,1,1)
betaset <- list(c(0,0,0),c(-1,0,1))

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





getEstimate <- function(Data,nrep,an,cl){
  lr.crit <- matrix(0.0,nr=nrep,ncol=3)
  lr.estimate <- matrix(0.0,nr=nrep,ncol=1)
  lr.size <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  for (k in 1:nrep){
    data <- Data[,k]
    out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
    out.h1 <- regpanelmixMaxPhi(y=data$Y,x=data$X, z = data$Z,parlist=out.h0$parlist,an=an,update.alpha = 1,parallel = TRUE,cl=cl)
    lr.estimate[k,] <- 2 * max(out.h1$penloglik - out.h0$loglik)
    # lr.crit <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z,cl=cl ,an=an, parallel = TRUE)$crit
  }# out.h1 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M+1,vcov.method = "none")
  
  # crit <- regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z,cl=cl , parallel = TRUE,nrep=1000)$crit
  #THIS IS FOR TEST PURPOSE
  crit.no <- order(result$est)[round(nrep/2)]
  data <- Data[,crit.no]
  out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
  crit <- regpanelmixCritBoot(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z,cl=cl ,an=an, parallel = TRUE)$crit
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
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(betaset)

NTset <- expand.grid(Nset,Tset)
Parset <- expand.grid(muset,alphaset,betaset)
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
      for (beta in betaset){
        cl <- makeCluster(7)
        
        t <- Sys.time()
        phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                   beta = beta, N = N, T = T, M = (M+1), p = p, q = q)
        
        phi.data.pair <- GenerateSample(phi,nrep)
        count <- count + 1
        Data = phi.data.pair$Data
        phi = phi.data.pair$phi
        # phi.data[[count]] <- phi.data.pair
        an <- anFormula(phi,M,N,T,q=1) #The an function according the the empirical regression
        # an <- 1.0
        print(an)
        result <- getEstimate(Data,nrep,an,cl)
        
        # print(result$nominal.size)
        result.f[r, count] <- result$nominal.size
        print(Sys.time() - t)
      }  
    }
    
  }
  
}

result.f = 100 * (  result.f )

write.csv(result.f,file="powerTestM2Regressor.csv")

