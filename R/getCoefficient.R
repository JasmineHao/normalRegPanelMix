#' @title get coefficients name
#' @export
#' @description get a list of coeficients from output
#' @param N int
#' @param T int
#' @param M number of types
#' @param p columns of Z
#' @param q columns of X
#' @return list of item
#' \item{alpha}
#' \item{mu}
#' \item{beta}
#' \item{sigma}
#' \item{gamma}
getCoefficient <- function(coefficients, N,T,M,p,q){
  alpha_hat <- as.vector(coefficients[1:M])
  mubeta_hat <- as.vector(coefficients[M+1:((q+1)*M)])
  mu_hat <- mubeta_hat[(0:(M-1))*(q+1)+1]
  if (q != 0){
    ind <- rep(0,0)
    for (mm in 0:(M-1)){
      ind <- append(ind,(2:(q+1))+mm*(q+1))
    }
    beta_hat <- t(matrix(mubeta_hat[ind],nr=q))
  }else{
    beta_hat = matrix(0)
  }
  sigma_hat <- coefficients[((q+2)*M+1):((q+3)*M)]
  gamma_hat <- coefficients[((q+3)*M+1):((q+3)*M+p)]
  
  return(list('alpha' = alpha_hat, 'mu' =mu, 
              'beta' = beta_hat, 'sigma' = sigma_hat, 'gamma' = gamma_hat ))
}