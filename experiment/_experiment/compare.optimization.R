library(nloptr)

q <- 0
n_lam <- (q+3)*(q+2) / 2
df <- q + 2
nrep <- 3
m <- 1
ndraw <- 200

# Get
K <- matrix(rnorm(n_lam*n_lam),nc=n_lam)
I <- t(K) %*% K
e <- eigen(I, symmetric=TRUE)
u <- t(e$vec %*% (t(e$vec) * sqrt(e$val)) %*% matrix(rnorm(nrep*n_lam*m), nrow=n_lam*m))
Z <- u %*% solve(I)

#Define the function
obj.func <- function(x,z,I){
  t(z - x) %*% I %*% (z - x) 
}

#Define the gradient
f.grad <- function(x,z,I){
  2*t(x - z) %*% I
}

#Check derivative
check.derivatives(.x=rep(1,3),func=obj.func,
                  func_grad = f.grad,check_derivatives_print='all',
                  z = Z[1,], I = I)

#slsqp, direct maximization

slsqp.sol <- matrix(0,nr=nrep,ncol=n_lam)
slsqp.time <- matrix(0,nr=nrep,ncol=1)
slsqp.draw <- matrix(0,nr=ndraw,ncol=n_lam)
slsqp.max.time <- matrix(0,nr=nrep,ncol=1)

hin <- function(x){
  c(x[1],x[2],x[1]*x[2] - x[3]^2)
}

heq <- function(x){
  x[3]^2 - x[1]*x[2]
}

for (i in 1:nrep){
  ptm <- proc.time()[3]
  z <- Z[i,]
  t <- slsqp(rnorm(3),obj.func,gr = f.grad,z=z,I=I,hin = hin)$par
  print(norm(as.matrix(hin(x))))
  
  for(j in 1:ndraw){
    
    
    slsqp.draw[j,] <- t
  }
  
  slsqp.draw.LR <- rowSums(slsqp.draw %*% I *slsqp.draw)
  max.ind <- which.max(slsqp.draw.LR)
  max.count <- abs(slsqp.draw.LR - slsqp.draw.LR[max.ind]) < 0.001
  slsqp.max.time[i] <- mean(max.count)
  slsqp.sol[i,] <- slsqp.draw[max.ind,]
  slsqp.time[i] <- proc.time()[3] - ptm
  print(slsqp.time[i])
}
