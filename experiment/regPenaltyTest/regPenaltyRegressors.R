library(normalRegPanelMix)
M <- 2
p <- 0
q <- 0
Nset <- c(100,500)
Tset <- c(2,5,10)
alphaset <- list(c(0.5,0.5),c(0.2,0.8))
muset <- list(c(-1,1),c(-0.5,0.5))
sigmaset <- list(c(1, 1), c(1.5, 0.75))
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(sigmaset)
an.data <- matrix(0,nr=nset,nc=3)

count <- 0
for (N in Nset){
  for (T in Tset){
    for (mu in muset){
      for (alpha in alphaset){
        for (sigma in sigmaset){
            count <- count + 1
            phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                       beta = beta, N = N, T = T, M = M, p = p, q = q)
            
            an <- anFormula(phi,M,N,T) * 0.5 
            an.data[count, ] <-
              cbind(an, phi$N,phi$T)

        }
      }
    }
  }
}

m_2.an <- mean(an.data[,1])

###########################################################
M <- 3
alphaset <- list(c(1/3,1/3,1/3),c(0.25,0.5,0.25))
muset 		<- list(c(-4, 0, 4), c(-4, 0, 5), c(-5, 0, 5), c(-4, 0, 6), c(-5, 0, 6), c(-6, 0, 6))
sigmaset <- list(c(1, 1, 1), c(0.75, 1.5, 0.75))

nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(sigmaset)
an.data <- matrix(0,nr=nset,nc=3)

count <- 0
for (N in Nset){
  for (T in Tset){
    for (mu in muset){
      for (alpha in alphaset){
        for (sigma in sigmaset){
          count <- count + 1
          phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                     beta = beta, N = N, T = T, M = M, p = p, q = q)
          
          an <- anFormula(phi,M,N,T) * 0.5 
          an.data[count, ] <-
            cbind(an, phi$N,phi$T)
          
        }
      }
    }
  }
}

m_3.an <- mean(an.data[,1])


###########################################################

M <- 4
alphaset <- list(c(0.25, 0.25, 0.25, 0.25))
muset 		<- list(c(-4,-1,1,4), c(-5,-1,1,5), c(-6,-2,2,6), c(-6,-1,2,5), c(-5,0,2,4), c(-6,0,2,4))
sigmaset <- list(c(1, 1, 1, 1), c(1, 0.75, 0.5, 0.25))

nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(sigmaset)
an.data <- matrix(0,nr=nset,nc=3)

count <- 0
for (N in Nset){
  for (T in Tset){
    for (mu in muset){
      for (alpha in alphaset){
        for (sigma in sigmaset){
          count <- count + 1
          phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                     beta = beta, N = N, T = T, M = M, p = p, q = q)
          
          an <- anFormula(phi,M,N,T) * 0.5 
          an.data[count, ] <-
            cbind(an, phi$N,phi$T)
          
        }
      }
    }
  }
}

m_4.an <- mean(an.data[,1])
