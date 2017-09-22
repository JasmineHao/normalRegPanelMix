library(normalRegPanelMix)
library(doParallel)
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





getEstimate <- function(Data,nrep,an,cl){
  lr.crit <- matrix(0.0,nr=nrep,ncol=3)
  lr.estimate <- matrix(0.0,nr=nrep,ncol=1)
  lr.size <- matrix(0.0,nr=nrep,ncol=1) #Nomimal size
  ptm <- proc.time()[1]

  for (k in 1:nrep){
    
    data <- Data[,k]
    out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
    
    out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=an,update.alpha = 1,parallel = TRUE,cl=cl)
    
    lr.estimate[k,] <- 2 * max(out.h1$penloglik - out.h0$loglik)
    }# out.h1 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M+1,vcov.method = "none")
    
    
    
  
  crit <- regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z,cl=cl , parallel = TRUE,nrep=1000)$crit
  for ( k in 1:nrep){
    lr.crit[k,] <- crit
    lr.size[k,] <- 1 * (lr.estimate[k,] > lr.crit[k,2])
  }
  # lr.estimate <- 2 * lr.estimate
  # asymp.out <- regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z,parallel = FALSE,deep.check = 0)
  # crit <- asymp.out$crit
  print(proc.time()[1] - ptm)
  return(list(est = lr.estimate , crit = lr.crit,nominal.size = apply(lr.size,2,mean)))
}

#GeneratePhiDataPairs
count <- 0
nrep <- 2000
phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset)
regression.data <- matrix(0,nr=(nset),nc=3)
for (N in Nset){
  for (T in Tset){
    for (mu in muset){
      for (alpha in alphaset){
        cl <- makeCluster(6)
        
        t <- Sys.time()
        phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                   beta = beta, N = N, T = T, M = M, p = p, q = q)
        
        phi.data.pair <- GenerateSample(phi,nrep)
        count <- count + 1
        Data = phi.data.pair$Data
        phi = phi.data.pair$phi
        # phi.data[[count]] <- phi.data.pair
        an <- anFormula(phi,M,N,T) * 0.5#The an function according the the empirical regression
        print(an)
        
        result <- getEstimate(Data,nrep,an,cl)
        
        normimal.size <- result$nominal.size

        # print(result$nominal.size)
        regression.data[count, ] <- 
           cbind(result$nominal.size,phi$T,phi$N)
        # write.csv(result$lr.estimate,file=paste("estimate.result", count))
        
        print(Sys.time() - t)
      }
    }
  }
}





# for(i in 1:count){
#   t.out <- Sys.time()
#   phi <- phi.data[[i]]$phi
#   Data = phi.data.pair$Data
#   an <- anFormula(phi,M,phi$N,phi$T) 
#   result <- getEstimate(Data,nrep,an)
#   regression.data[i,] <- cbind(result$nominal.size,phi$N,phi$T)
#   print(cbind(result$nominal.size,phi$N,phi$T,phi$alpha,phi$mu,phi$sigma))
#   print(i)
#   print(Sys.time() - t.out)
#   
# }
write.csv(regression.data,file="estimate.result")

# registerDoParallel(detectCores())
# t <- Sys.time()
# foreach(i=1:2, .combine = rbind)%dopar%{
#   library(normalRegPanelMix)
#   phi <- phi.data[[i]]$phi
#   Data = phi.data.pair$Data
#   an <- anFormula(phi,M,phi$N,phi$T)
#   result <- getEstimate(Data,nrep,an)
#   cbind(result$nominal.size,phi$N,phi$T,phi$alpha,phi$mu,phi$sigma)
# }
# print(Sys.time() - t)
