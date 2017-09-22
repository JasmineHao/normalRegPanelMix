

#  get.t.from.lambda
# lambda Take in lambda, list of mu, sigma, beta
# @return the second order expansion
get.t.from.lambda <- function(lambda){
  mu_j <- lambda[1]
  sigma_j <- lambda[2]
  beta_j <- lambda[length(lambda)-2:length(lambda)]
  q_j <- length(lambda)
  q <- q_j - 2
  t <- rep(0,q_j * (q_j - 1) /2 + q_j)
  t[1] <- mu_j^2
  t[2] <- mu_j * (sigma_j^2)
  t[3] <- sigma_j^4
  if (q > 0){
    t[4:(3+q)] <- mu_j * beta_j
    t[(4+q):(3+2*q)] <- (sigma_j^2) * beta_j
    t[(4+2*q):(3+3*q)] <- beta_j^2
  }
  # :(3+3*q+(q-1)*q/2)
  if (q > 1) {
    qq <-  1
    for (j in 1:(q-1)) {
      for (i in (j+1):q) {
        t[(3+3*q+qq)] <- 2*beta_j[j]*beta_j[i]
        qq <- qq+1
      }
    }
  }
  return(t)
}

hermite <- function(Z, sigma)
# Computes the normalized Hermite polynomial for computing
# the critical values and p-values of the modified EM test
# Input
#   Z (n by m) : normalized data
#  sigma (m by 1): parameters
# Output
#   H (n by m by 4): normalized Hermite polynomials
{
n <- nrow(Z)
m    <- length(sigma)

H <- array(0,dim=c(n,m,4))
H[,,1] <- Z/sigma
H[,,2] <- t(t(Z^2-1)/2/sigma^2)
H[,,3] <- t(t(Z^3-3*Z)/6/sigma^3)
H[,,4] <- t(t(Z^4-6*Z^2+3)/24/sigma^4)
return(H)
}  # end function hermite


tKR <- function (a,b) {
  # Computes the transpose of the Khatri-Rao product of t(a) and t(b)
  n <- nrow(a)
  k <- ncol(a)
  KR <- matrix(unlist(lapply(1:k, function(i) (a[,i]*b))),nrow=n)
  KR
}


coef.to.list <- function(coefficients, z = NULL) {
# 　Convert coefficients to list
len     <- length(coefficients)
p       <- 0
gam   <- NULL

if (!is.null(z)) {
  z <- as.matrix(z)
  p <- ncol(z)
  gam <- coefficients[(len-p+1):len]
}

m <- (len-p)/3
if (round(m) != m) {
  stop("The dimension of the coefficients is incompatible with z. Please check the data.")
}

param   <- matrix(coefficients[1:(len-p)], nrow=m, ncol=3)
alpha   <- param[, 1]
mu      <- param[, 2]
sigma   <- param[, 3]

a = list(alpha = alpha, mu = mu, sigma = sigma, gam = gam)

a

}

#' @description Generates mixed normal random variables with regressor x
#' @title rnormregpanelmix
#' @name rnormregpanelmix
#' @param n The number of observations
#' @param t Number of time periods
#' @param x n by k-1 matrix that does NOT include a constant
#' @param alpha m by 1 vector that represents proportions of components
#' @param mubeta k by m matrix that represents (mu times k regression coefficients) on x for m components
#' @param sigma m by 1 vector that represents sd of components
#' @return n by 1 vector that is formed by regressor x
rnormregpanelmix <- function (n, t, x = NULL, alpha, mubeta, sigma) {
  # Generates mixed normal random variables with regressor x
  # Input
  #  n : number of observations
  #   x : (n by k-1) matrix NOT including a constant
  #   alpha  : m-vector
  #  mubeta  : k by m matrix
  #  sigma  : m-vector
  # Output
  #  y : n by 1 vector
  # if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
  #   set.seed(normalregMix.test.seed)
  
  m     <- length(alpha)
  nt    <- n*t
  mubeta   <- matrix(mubeta, ncol=m)
  y <- rep(0,nt)
  if (!is.null(x)){
    x <- as.matrix(x)
    if (nrow(x) != nt) { stop("y and x must have the same number of rows.") }
    x1   <- cbind(1,x)
    ii   <- sample(m, n, replace=TRUE, prob=alpha)
    
    for (nn in 1:n){
      y[(t*(nn-1)+1):(t*nn)]  <- rnorm(t, mean = x1[(t*(nn-1)+1):(t*nn),]%*% mubeta[, ii[nn]] , sd = sigma[ii[nn]])
    }
  } else {
    ii   <- sample(m, n, replace=TRUE, prob=alpha)
    for (nn in 1:n){
      y[(t*(nn-1)+1):(t*nn)]  <- rnorm(t, mean = mubeta[, ii[nn]] , sd = sigma[ii[nn]])
    }
  }
  
  y
  
} 



#' @description Computes omega_{j|i} defined in (2.1) of Maitra and Melnykov (2010)
#' @export
#' @title omega.ji
#' @name omega.ji
#' @param phi_i 3 by 1 column consisting of alpha, mu, sigma of ith component
#' @param phi_j 3 by 1 column consisting of alpha, mu, sigma of jth component
#' @return omega_{j|i}
#' @references Maitra, R., and Melnykov, V. (2010)
#' Simulating Data to Study Performance of Finite Mixture Modeling and Model-Based Clustering Algorithms,
#' \emph{Journal of Computational and Graphical Statistica},
#' \bold{19}, 354--376.
# Returns a misclassification rate omega_ji given two components i, j,
# i.e. the probability of choosing component j where
# the true model is ith component.
omega.ji <- function(phi_i, phi_j) {
  alpha_i <- phi_i[1]
  alpha_j <- phi_j[1]
  mu_i <- phi_i[2]
  mu_j <- phi_j[2]
  sigma_i <- phi_i[3]
  sigma_j <- phi_j[3]

  a <- (1/sigma_j^2 - 1/sigma_i^2)
  b <- mu_i / sigma_i^2 - mu_j / sigma_j^2
  c <- mu_j^2 / sigma_j^2 - mu_i^2 / sigma_i^2

  if (sigma_i == sigma_j)
    if (mu_i > mu_j)
      omega_ji = pnorm((2 * log(alpha_j/alpha_i) - c)/(2*b),
                       mean = mu_i, sd = sigma_i)
  else
    omega_ji = 1 - pnorm((2 * log(alpha_j/alpha_i) - c)/(2*b),
                         mean = mu_i, sd = sigma_i)
  else {
    d <- 2 * log(alpha_j * sigma_i / (alpha_i * sigma_j)) - c + (b^2 / a)
    da <- max(d/a, 0)
    if (sigma_i > sigma_j)
      omega_ji = pnorm(sqrt(da)-b/a, mean = mu_i, sd = sigma_i) -
      pnorm(-sqrt(da)-b/a, mean = mu_i, sd = sigma_i)
    else
      omega_ji = 1 +
      pnorm(-sqrt(da)-b/a, mean = mu_i, sd = sigma_i) -
      pnorm(sqrt(da)-b/a, mean = mu_i, sd = sigma_i)
  }
  return (omega_ji)
}

#' @description Computes omega_{12} defined in Maitra and Melnykov (2010)
#' @export
#' @title omega.12
#' @name omega.12
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gam
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gam_1, ..., gam_m))
#' @return The misclassification rate omega_ij
#' @references Maitra, R., and Melnykov, V. (2010)
#' Simulating Data to Study Performance of Finite Mixture Modeling and Model-Based Clustering Algorithms,
#' \emph{Journal of Computational and Graphical Statistica},
#' \bold{19}, 354--376.
omega.12 <- function(parlist)
  # Computes omega_{12} for testing H_0:m=2 against H_1:m=3
{
  phi1 <- c(alpha = parlist$alpha[1], mu = parlist$mu[1], sigma = parlist$sigma[1])
  phi2 <- c(alpha = parlist$alpha[2], mu = parlist$mu[2], sigma = parlist$sigma[2])

  part1 <- omega.ji(phi1, phi2)
  part2 <- omega.ji(phi2, phi1)

  return((part1 + part2) / 2)
}  # end function omega.12


#' Computes omega_{12} and omega_{23} defined in Maitra and Melnykov (2010)
#' @export
#' @title omega.123
#' @name omega.123
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @return A 2 by 1 vector whose first element is omega_12 and second element is omega_23
#' @references Maitra, R., and Melnykov, V. (2010)
#' Simulating Data to Study Performance of Finite Mixture Modeling and Model-Based Clustering Algorithms,
#' \emph{Journal of Computational and Graphical Statistica},
#' \bold{19}, 354--376.
omega.123 <- function(parlist)
{
  phi1 <- c(alpha = parlist$alpha[1], mu = parlist$mu[1], sigma = parlist$sigma[1])
  phi2 <- c(alpha = parlist$alpha[2], mu = parlist$mu[2], sigma = parlist$sigma[2])
  phi3 <- c(alpha = parlist$alpha[3], mu = parlist$mu[3], sigma = parlist$sigma[3])

  part1 <- omega.ji(phi1, phi2)
  part2 <- omega.ji(phi2, phi1)
  w12 <- (part1 + part2)/2

  part3 <- omega.ji(phi2, phi3)
  part4 <- omega.ji(phi3, phi2)
  w23 <- (part3 + part4)/2

  return(c(w12, w23))

}  # end function omega.123

#' @description Computes omega_{12}, omega_{23}, and omega_{34} defined in Maitra and Melnykov (2010)
#' @export
#' @title omega.1234
#' @name omega.1234
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @return A 3 by 1 vector consisting of omega_12, omega_23, and omega_34
#' @references Maitra, R., and Melnykov, V. (2010)
#' Simulating Data to Study Performance of Finite Mixture Modeling and Model-Based Clustering Algorithms,
#' \emph{Journal of Computational and Graphical Statistica},
#' \bold{19}, 354--376.
omega.1234 <- function(parlist)
{
  phi1 <- c(alpha = parlist$alpha[1], mu = parlist$mu[1], sigma = parlist$sigma[1])
  phi2 <- c(alpha = parlist$alpha[2], mu = parlist$mu[2], sigma = parlist$sigma[2])
  phi3 <- c(alpha = parlist$alpha[3], mu = parlist$mu[3], sigma = parlist$sigma[3])
  phi4 <- c(alpha = parlist$alpha[4], mu = parlist$mu[4], sigma = parlist$sigma[4])

  part1 <- omega.ji(phi1, phi2)
  part2 <- omega.ji(phi2, phi1)
  w12 <- (part1 + part2)/2

  part3 <- omega.ji(phi2, phi3)
  part4 <- omega.ji(phi3, phi2)
  w23 <- (part3 + part4)/2

  part5 <- omega.ji(phi3, phi4)
  part6 <- omega.ji(phi4, phi3)
  w34 <- (part5 + part6)/2

  return(c(w12, w23, w34))

}  # end function omega.1234

coef.to.list <- function(coefficients, z = NULL) {
  # Convert coefficients to list
  len     <- length(coefficients)
  p       <- 0
  gam   <- NULL

  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    gam <- coefficients[(len-p+1):len]
  }

  m <- (len-p)/3
  if (round(m) != m) {
    stop("The dimension of the coefficients is incompatible with z. Please check the data.")
  }

  param   <- matrix(coefficients[1:(len-p)], nrow=m, ncol=3)
  alpha   <- param[, 1]
  mu      <- param[, 2]
  sigma   <- param[, 3]

  a = list(alpha = alpha, mu = mu, sigma = sigma, gam = gam)

  a

}
