#' Add together two numbers.
#'
#' @export
#' @param alpha m vector.
#' @param mu m vector.
#' @param sigma m vector.
#' @param gamma p vector.
#' @param beta q vector
#' @param N int
#' @param T int
#' @param M number of types
#' @param p columns of Z
#' @param q columns of X
#' @return list of items
#' \item{Y}
#' \item{X}
#' \item{Z}
#' @examples
#'generateData(alpha,mu,sigma,gamma,beta,N,T,M,p,q)
# The example of X, Z

# X = matrix(rnorm(N*T*q),nc=q)
# Z = matrix(rnorm(N*T*p),nc=p)
generateData <- function(alpha,mu,sigma,gamma,beta,N,T,M,p,q,X=NULL,Z=NULL){
  sprintf("N = %d",N)
  sprintf("T = %d",T)

  R <- matrix(0,nrow=N,ncol=M)
  if (sum(alpha) != 1){
    alpha <- alpha / sum(alpha)
  }
  
  if (length(alpha) != M | length(mu) != M){
    stop('M must be the size of alpha')
  }
  
  prior <- runif(N)
  alpha.cum <- c(0,cumsum(alpha))
   
  if (M > 1){
    for (m in 1:M){
      lb <-  alpha.cum[m]
      ub <- alpha.cum[m+1] 
      
      R[,m] <- 1 * ((prior <= ub) & (prior > lb))
      }
  }else{
    R <- matrix(1,nrow=N,ncol=M)
  }
  Y = matrix(0,T,N)

  if ((q != 0) && (is.null(X))){
    X = matrix(rnorm(N*T*q),nc=q)
  }
  if ((p != 0) && (is.null(Z))){
    X = matrix(rnorm(N*T*p),nc=p)
  }
  

  mu_R <- R %*% mu
  sigma_R <- R %*% sigma
  u <- matrix(rnorm(N*T),nrow=T,ncol=N)

  for (nn in 1:N) {
    y_nn <- rep(0,T)
    y_nn <- mu_R[nn] + sigma_R[nn] * u[,nn]
    if (q > 1 ){
      beta_R <- R %*% beta
      y_nn = y_nn + X[(T*(nn-1)+1) : (T*nn),] %*% beta_R[nn,]
     }else if (q == 1){
       beta_R <- R %*% beta
       y_nn = y_nn + X[(T*(nn-1)+1) : (T*nn),] * beta_R[nn]
     }else{}

    if (p > 1){
      y_nn = y_nn + Z[(T*(nn-1)+1) : (T*nn),] %*% gamma
    }else if (p ==1){
      y_nn = y_nn + Z[(T*(nn-1)+1) : (T*nn),] * gamma
    }else{ }
    Y[,nn] <- y_nn
  }
  if (p == 0){Z = NULL}
  if (q == 0){X = NULL}
  
  return(list(Y=Y,Z=Z,X=X))
  }
