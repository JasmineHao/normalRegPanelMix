library(normalRegPanelMix)
library(foreach)
library(Rmpi)
#Generate Data
set.seed(123456)

Nset <- c(100,500)
Tset <- c(2,10)

anset <- c(0.05,0.1,0.15,0.2,0.3,0.4)

alphaset <- list(c(1/3,1/3,1/3),c(0.25,0.5,0.25))
muset 		<- list(c(-4, 0, 4), c(-4, 0, 5), c(-5, 0, 5), c(-4, 0, 6), c(-5, 0, 6), c(-6, 0, 6))
sigmaset <- list(c(1, 1, 1), c(0.75, 1.5, 0.75))


#The parameters that are fixed
M <- 3 #Number of Type
p <- 0 #Number of Z
q <- 0 #Number of X
gamma <- matrix(0)
beta <- matrix(0)

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


#GeneratePhiDataPairs
count <- 0
nset <- length(Nset) * length(Tset) * length(muset) * length(alphaset) * length(sigmaset)
omega.data <- matrix(0,nr=(nset*length(anset)),nc=1)


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
            count <- count+1
            phi = list(alpha = alpha,mu = mu,sigma = sigma, gamma = gamma,
                       beta = beta, N = N, T = T, M = M, p = p, q = q)
            
            print(an)
            omega <- GetMisclTerm(phi)
            omega.data[count,1] <- omega
            print(Sys.time() - t)
            
          }
        }
      }
    }
  }
}


colnames(regression.data) <- c("nom.size", "an" ,"N","T", "omega")
# Begin estimation

stargazer(regression.data)
