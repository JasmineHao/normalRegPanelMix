library(normalRegPanelMix)
#Generate Data
Nset <- c(100,500)
Tset <- c(2,5,10)
alphaset <- list(c(0.5,0.5),c(0.2,0.8))
muset <- list(c(-1,1),c(-0.5,0.5))
sigma <- c(0.8,1.2)
anset <- c(0.05,0.1,0.15,0.2,0.3,0.4)

#The parameters that are fixed
M <- 2 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X
gamma <- matrix(0)
beta <- matrix(0)

#The parameters that are not fixed
# N <- 100 #Number of people
# T <- 5 #Time periods, it's important that T > 1
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



ptm <- proc.time()[1]

getEstimate <- function(Data,nrep,anset){
  lr.crit <- matrix(0.0,nr=nrep,ncol=3)
  lr.estimate <- matrix(0.0,nr=nrep,ncol=length(anset))
  lr.size <- matrix(0.0,nr=nrep,ncol=length(anset)) #Nomimal size
  ptm <- proc.time()[1]
  for (k in 1:nrep){
    
    data <- Data[,k]
    
    out.h0 <- normalpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M,vcov.method = "none")
    
    for (j in 1 : length(anset) ){
      an = anset[j]
      out.h1 <- normalpanelmixMaxPhi(y=data$Y,parlist=out.h0$parlist,an=an,update.alpha = 1)
      lr.estimate[k,j] <- 2 * max(out.h1$penloglik - out.h0$loglik)
    }# out.h1 <- regpanelmixPMLE(y=data$Y,x=data$X, z = data$Z,m=M+1,vcov.method = "none")
    
    
    
  }
  crit <- regpanelmixCrit(y=data$Y, x=data$X, parlist=out.h0$parlist, z = data$Z,parallel = FALSE)$crit
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
nrep <- 100
phi.data <- list()
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset)
regression.data <- matrix(0,nr=nset*length(anset),nc=4)

for (N in Nset){
  for (T in Tset){
    for (mu in muset){
      for (alpha in alphaset){
        phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                   beta = beta, N = N, T = T, M = M, p = p, q = q)
        
        phi.data.pair <- GenerateSample(phi,nrep)
        # phi.data[[count]] <- phi.data.pair
        Data = phi.data.pair$Data
        # phi = phi.data.pair$phi
        result <- getEstimate(Data,nrep,anset)
        
        normimal.size <- result$nominal.size
        count <- count + 1
        print(count)
        regression.data[((count-1) *length(anset) + 1) : (count *length(anset) ), ] <- 
          cbind(result$nominal.size,anset,phi$T,phi$N)
        
        
        
        write.csv(result$lr.estimate,file=paste("estimate.result", count))
        
        
      }
    }
  }
}

# Begin estimation
# regression.data <- matrix(0,nr=count*length(anset),nc=4)
# for (i in 1:3){
#   Data <- phi.data[[i]]$Data
#   phi  <- phi.data[[i]]$phi
#   result <- getEstimate(Data,nrep,anset)
#   
#   normimal.size <- result$nominal.size
#   regression.data[((i-1) *length(anset) + 1) : (i *length(anset) ), ] <- 
#     cbind(result$nominal.size,anset,phi$T,phi$N)
#   write.csv(result$lr.estimate,file=paste("estimate.result", i))
#   # p.hat <- mean(result$est > result$crit[2])
# }

write.csv(regression.data,file="regression.csv")
fit.data <- list(y = log(regression.data$V1 / (0.1 - regression.data$V1) ) , x2 = log( regression.data$V2/ (1 - regression.data$V2 ))  , x1 =1 /  regression.data$V3,x3 = 1 / regression.data$V4) 
fit.data$y <- replace(fit.data$y,fit.data$y == -Inf,-2)
lfit <- lm(y ~ x1 + x2 + x3, data=fit.data,na.action = na.omit)
