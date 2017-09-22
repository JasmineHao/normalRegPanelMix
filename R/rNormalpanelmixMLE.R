#' @export
#' @description Estimates parameters of a finite panel mixture of univariate normals by
#' penalized maximum log-likelhood functions.
#' @title rnormalpanelmixlikelihood
#' @name rnormalpanelmixlikelihood
#' @param y nt by 1 vector of data
#' @param x nt by q matrix of data for x (if exists)
#' @param z nt by p matrix of regressor associated with gamma
#' @param b (p+q+3)*m by 1 vector, b <- (alpha,mu,beta,sigma,q)
#' @param m The number of components in the mixture
#' @param q The number of variable in x
#' @param p The number of variable in z
#' @return  log-likelihood


rNormalPanelMixMLE <- function (alpha,mu,sigma,beta,gamma,y,x,z,M,p,q){
  t <- nrow(y)
  n <- ncol(y)
  y <- as.vector(y)
  
  alpha <- alpha / sum(alpha)
  ll <- - n*t*log(2*pi)/2
  l_j <- matrix(0,nr=M,nc=n)
  r <- matrix(0,nr=M,nc=n)
  mubeta <- cbind(mu,beta)
  x1 <- cbind(1,x)

  for (mm in 1:M){
    if (p != 0){
      ytilde <- y - z %*% gamma
    }else{
      ytilde <- y
    }
    
    if (q != 0){
      ytilde <- ytilde - x1 %*% mubeta[mm,]
      
    }else{
      ytilde <- ytilde - mu[mm]
    }
    r_mm <- 0.5*(ytilde/sigma[mm])^2
    
    r_mm <- matrix(r_mm,nr=t)
    r[mm,] <- colSums(r_mm) + t*log(sigma[mm])
    #Avoid prod of small numbers, use log sum
  }
  
  r_min <- apply(r,2,min)
  for (mm in 1:M){
    r[mm,] <- r[mm,] - r_min
    l_j[mm,] <- alpha[mm]*exp(-r[mm,])
    
    #l_j[mm,] <- exp(log(l_j[mm,]) - r_min)
  }
  # Likelihood
  
  ll <- ll + sum( log(colSums(l_j)) - r_min )
  return(ll)
}
