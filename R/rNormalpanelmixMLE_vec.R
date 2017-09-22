#' @export
#' @description Estimates parameters of a finite panel mixture of univariate normals by
#' penalized maximum log-likelhood functions.
#' @title rNormalPanelMixMLE_vec
#' @name rNormalPanelMixMLE_vec
#' @param alpha m by 1 vector
#' @param mu m by 1 vector
#' @param sigma m by 1 vector
#' @param gamma p by 1 vector
#' @param beta m by q matrix
#' @param y nt by 1 vector of data
#' @param x nt by q matrix of data for x (if exists)
#' @param z nt by p matrix of regressor associated with gamma
#' @param m The number of components in the mixture
#' @param p The number of variable in z
#' @param q The number of variable in x
#' @return  log-likelihood


rNormalPanelMixMLE_vec <- function (b0,y,x,z,M,p,q){
  t <- nrow(y)
  n <- ncol(y)
  y <- as.vector(y)
  
  alpha <- b0[1:M]
  mu <- b0[(M+1):(2*M)]
  sigma <- b0[(2*M+1):(3*M)]
  if (q != 0){
  beta <- b0[(3*M+1):((3+q)*M)]
  }else{
    beta <- NULL}
  if (p != 0){
  gamma <- b0[((3+q)*M+1):((3+q)*M+p)]
  }else{
  gamma <- NULL  
  }
  alpha <- alpha / sum(alpha)
  ll <- - n*t*log(2*pi)/2
  l_j <- matrix(0,nr=M,nc=n)
  r <- matrix(0,nr=M,nc=n)
  if (q != 0){
    mubeta <- cbind(mu,matrix(beta,nr=M))
    x1 <- cbind(1,x)
  }else{
    mubeta <- mu
  }
  
  # y <- as.vector(y)
  
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
