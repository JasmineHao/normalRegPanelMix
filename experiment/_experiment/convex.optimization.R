library(nloptr)

q <- 0
n_lam <- (q+3)*(q+2) / 2
df <- q + 2
nrep <- 30
m <- 1
ndraw <- 5

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

hin <- function(x){
  c(x[1],x[2],x[1]*x[2] - x[3]^2)
}

hin_2 <- function(x){
  c(x[1],x[2])
}

heq <- function(x){
  x[3]^2 - x[1]*x[2]
}

count <- 0
for (i in 1:nrep){
  ptm <- proc.time()[3]
  z <- Z[i,]
  t <- slsqp(rnorm(3),obj.func,gr = f.grad,z=z,I=I,hin = hin)$par
  cond <- hin(t)
  
  if (!(abs(cond[3]) < 1e-5)){
    count <- count + 1
    
    cseq <- lapply( t , function(x){seq(-x,x,length=ndraw)})
    draws <- as.matrix(expand.grid(cseq))
    
    eq.optimization <-  function(x){slsqp(x,obj.func,gr = f.grad,z=z,I=I,hin = hin_2,heq = heq)$par}
    eq.result <- t(apply(draws,1,eq.optimization))
    draw.LR <- rowSums(eq.result %*% I *eq.result)
    t <- eq.result[which.max(draw.LR),]
  }
  
  
  slsqp.sol[i,] <- t
  slsqp.time[i] <- proc.time()[3] - ptm
  
}
#############################################################
#Lambda optimization
#############################################################

slsqp.LR <- rowSums(slsqp.sol %*% I *slsqp.sol)

mean(slsqp.LR)

mean(rowSums(Z %*% I * Z))

mean(slsqp.time)

