library(normalRegPanelMix)
#Generate Data
N <- 100 #Number of people
T <- 2 #Time periods, it's important that T > 1
M <- 1 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X

alpha <- c(1) #Probability
mu <- c(0.5)

sigma <- c(1.2)
gamma <- matrix(0)
beta <- matrix(0)

if (q != 0){
  parlist <- list('alpha' = alpha, 'mubeta' = t(cbind(mu,beta)), 'sigma' = sigma, 'gam' = gamma)
}else{
  parlist <- list('alpha' = alpha, 'mubeta' = mu, 'sigma' = sigma, 'gam' = gamma)
}
nrep = 200

lr.crit <- matrix(0,nr=nrep,ncol=3)
lr.estimate <- rep(0,nrep)

ptm <- proc.time()[1]
Data <- replicate(nrep,generateData(alpha,mu,sigma,gamma,beta,N,T,M,p,q))



for (k in 1:nrep){
  ptm <- proc.time()[1]
  data <- Data[,k]
  
  out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
  

  out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=0.1)
  # out.h1 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M+1,vcov.method = "none")
  
  lr.estimate[k] <- min(out.h1$penloglik) - out.h0$loglik
  # lr.estimate[k] <- min(out.h1$loglik - out.h0$loglik)
  # lr.crit[k,] <- regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z,parallel = FALSE)$crit
  
  
  print(proc.time()[1] - ptm)
  
  
  # out.h0 <- normalmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
  # print(normalmixCrit(y=data$Y,  parlist=out.h0$parlist)$crit)
  # out.h1 <- normalmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=0.5)
}
 lr.estimate <- 2 * lr.estimate



#  
# print(mean(lr.estimate > lr.crit[,1]))
# 
# print(mean(lr.estimate > lr.crit[,2]))
# 
# print(mean(lr.estimate > lr.crit[,3]))
 crit <- regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z,parallel = FALSE)$crit
 print(mean(lr.estimate > crit[1]))
 print(mean(lr.estimate > crit[2]))
 print(mean(lr.estimate > crit[3]))
 