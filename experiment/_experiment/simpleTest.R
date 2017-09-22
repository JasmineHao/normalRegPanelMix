library(normalRegPanelMix)
#Generate Data
N <- 100 #Number of people
T <- 2 #Time periods, it's important that T > 1
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X

alpha <- c(0.5,0.5) #Probability
mu <- c(-0.5,0.5)

sigma <- c(1.2,0.8)
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
  lr.crit[k,] <- regpanelmixCrit(y=data$Y, x=data$X, parlist=parlist, z = data$Z,parallel = FALSE)$crit
  print(lr.crit[k,])
  out.h0 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
  # rNormalPanelMixMLE(out.h0$parlist$alpha,out.h0$parlist$mu,out.h0$parlist$sigma,out.h0$parlist$beta,out.h0$parlist$gamma,data$Y,data$X,data$Z,M,p,q)
  # out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=1/N)
  out.h1 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M+1,vcov.method = "none")
  lr.estimate[k] <- out.h1$penloglik - out.h0$loglik
  # lr.estimate.pen[k] <- out.h1$penloglik - out.h0$penloglik
  print(2 * lr.estimate[k])
  print(k)
  print(proc.time()[1] - ptm)
}
lr.estimate <- 2 * lr.estimate
lr.estimate.pen <- 2 * lr.estimate.pen


print(mean(lr.estimate > lr.crit[,1]))
print(mean(lr.estimate.pen > lr.crit[,1]))

print(mean(lr.estimate > lr.crit[,2]))
print(mean(lr.estimate.pen > lr.crit[,2]))

print(mean(lr.estimate > lr.crit[,3]))
print(mean(lr.estimate.pen > lr.crit[,3]))
