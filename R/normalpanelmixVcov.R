#' @description Computes the variance-covariance matrix of the MLE of
#' m-component normal mixture.
#' @title normalpanelmixVcov
#' @name normalpanelmixVcov
#' @param y n by 1 vector of data
#' @param coefficients (alpha_1, ..., alpha_m, mu_1, ..., mu_m, sigma_1, ..., sigma_m, gam)
#' @param z n by p matrix of regressor associated with gamma
#' @param vcov.method Method used to compute the variance-covariance matrix,
#' one of \code{"Hessian"} and \code{"OPG"}. #' The default option is \code{"Hessian"}.
#' When \code{method = "Hessian"}, the variance-covarince matrix is
#' estimated by the Hessian using the formula given in Boldea and Magnus (2009).
#' When \code{method = "OPG"}, the outer product of gradients is used.
#' @return The variance-covariance matrix of the MLE of
#' m-component normal mixture given the data and coefficients.
#' @references   Boldea, O. and Magnus, J. R. (2009)
#' Maximum Likelihood Estimation of the Multivariate Normal Mixture Model,
#' \emph{Journal of the American Statistical Association},
#' \bold{104}, 1539--1549.
normalpanelmixVcov <- function(y, coefficients, z = NULL, vcov.method = c("Hessian", "OPG"))
{
t  <- nrow(y)
n  <- ncol(y)
y <- as.vector(y)
len <- length(coefficients)
p <- 0
gam  <- NULL
vcov.method <- match.arg(vcov.method)

if (!is.null(z)) {
  z <- as.matrix(z)
  p <- ncol(z)
  gam <- coefficients[(len-p+1):len]
}

m <- (len-p)/3
if (round(m) != m) {
  stop("The dimension of the coefficients is incompatible with z. Please check the data.")
}

alpha   <- coefficients[1:m]
mu      <- coefficients[(m+1):(2*m)]
sigma   <- coefficients[(2*m+1):(3*m)]

if (m == 1) {
  if (is.null(z)) {
    I <- n*diag(c(1/sigma^2, 2/sigma^2))  # information matrix
  } else {
    z1 <- cbind(1, z)
    I <- matrix(0, nrow=p+2, ncol=p+2)
    I[1:(p+1), 1:(p+1)] <- t(z1) %*% z1/sigma^2
    I[(p+2), (p+2)] <- n*2/sigma^2
    s.1 <- c(1,(p+2), 2:(p+1))
    I <- I[s.1, s.1]
  }
  vcov  <- solve(I)
  # Because the variance is parameterized as sigma^2, we convert it to sigma
  c.mat.vec <- c(1, (1/sigma^(1/2))/2, rep(1, p))
  vcov <- diag(c.mat.vec) %*% vcov %*% diag(c.mat.vec)

} else { # end of if (m == 1)
  # m >= 2

  # Compute posterior probabilities, and adjust y if z is present
  sigma0  <- rep(1, m)  # dummy
  mu0     <- double(m)  # dummy
  an      <- 1/n  # penalty term for variance
  h       <- 0
  tau     <- 0.5
  k       <- 0
  epsilon <- 1e-08
  maxit = 2
  ninits = 1

  b <- matrix( rep( coefficients, ninits), ncol = ninits)
  if (is.null(z)) {
    out.p <- cppnormalpanelmixPMLE(b, y, matrix(0),  mu0, sigma0, m, p ,t , an, maxit, ninits, epsilon, tau, h, k)

  } else {
    out.p <- cppnormalpanelmixPMLE(b, y, z,  mu0, sigma0, m, p, t, an, maxit, ninits, epsilon, tau, h, k)
    # Adjust y
    y <- y - z %*% gam
  }
  post <- matrix(out.p$post, nrow=n)

  p1 <- seq(1, (2*m-1), by=2) # sequence of odd numbers, 1,3,...,2*m-1
  p2 <- seq(2, (2*m), by=2)    # sequence of even numbers, 2,4,...,2*m

  # Matrices used in computing vcov

  a <- diag(1/alpha[-m], nrow=m-1, ncol=m-1)
  a <- cbind(a, -1/alpha[m])  # m-1 by m matrix of a_i's
  abar <- a %*% t(post) # m-1 by n

  Z0 <- t((t(matrix(rep.int(y, m), ncol=m))-mu)/sigma)  # normalized data, n by m
  f <- t(t(exp(-Z0^2/2)/sqrt(2*pi))/sigma)      # pdf, n by m
  phi <- t(t(f)*alpha)                # n by m
  f0 <- rowSums(phi)                  # data pdf, n by 1

  vinv <- 1/(sigma*sigma)

  b <- t(t(Z0)/sigma)  # n by m
  B <- t(vinv - t(b*b))  # n by m

  c0 <- array(0, dim=c(n, m, 2))
  c0[, , 1] <- b
  c0[, , 2] <- -B/2

  # Computes Hessian-based I
  if (vcov.method == "Hessian")  {
    other.method = "OPG"
    C0 <- array(0, dim=c(n, m, 2, 2))
    C0[, , 1, 1] <- t(matrix(vinv, nrow=m, ncol=n))  # n by m
    C0[, , 2, 1] <- C0[, , 1, 2] <- t(t(b)*vinv)    # n by m
    C0[, , 2, 2] <- t((vinv -2*t(B))*vinv)/2       # n by m

    Q.pi <- - abar %*% t(abar)  # m-1 by m-1
    Q.pi.theta <- matrix(0, nrow=m-1, ncol=2*m)  # m-1 by 2m
    for (i in 1:m){
      zi <- a[, i] - abar  # m-1 by n
      wi <- c0[, i, ]*post[, i]  # n by 2
      Q.i <- colSums(tKR(wi, t(zi)))  # 2*(m-1) vector
      # first m-1 elements correspond to mu x pi
      # second m-1 elements correspond to sigma x pi
      Q.pi.theta[, i] <- Q.i[1:(m-1)]
      Q.pi.theta[, m+i] <- Q.i[m:(2*(m-1))]
    }

    Q.theta <- matrix(0, nrow=2*m, ncol=2*m)
    for (i in 2:m){ # off-diagonal blocks
      for (j in 1:(i-1)){
        wi  <- c0[, i, ]*post[, i]
        wj  <- c0[, j, ]*post[, j]
        Q.ij <- - colSums(tKR(wi, wj))
        Q.theta[(2*i-1):(2*i), (2*j-1):(2*j)] = t(matrix(Q.ij, nrow=2, ncol=2))
      }
    }

    Q.theta <- Q.theta + t(Q.theta)
    for (i in 1:m){ # diagonal blocks
      C.ii <- array(C0[, i, , ], dim=c(n, 2, 2))
      Q.ii.1 <- apply(C.ii*post[, i], c(2, 3), sum)
      w.ii <- tKR(c0[, i, ], c0[, i, ])*post[, i]*(1-post[, i])
      Q.ii.2 <- matrix(colSums(w.ii), nrow=2, ncol=2)
      Q.theta[(2*i-1):(2*i), (2*i-1):(2*i)] <- -Q.ii.1 + Q.ii.2
    }
    # odd rows and columns of Q.theta correspond to mu
    # even rows and columns of Q.theta correspond to sigma
    Q.theta <- Q.theta[c(p1,p2), c(p1,p2)] # first block = wrt mu, second blosk = wrt sigma

    dimI <- m-1+2*m
    I <- matrix(0, nrow=dimI, ncol=dimI)
    I[1:(m-1), 1:(m-1)] <- - Q.pi
    I[1:(m-1), m:dimI]  <- - Q.pi.theta
    I[m:dimI, 1:(m-1)]  <- - t(Q.pi.theta)
    I[m:dimI, m:dimI]   <- - Q.theta

    if (!is.null(z)) {
      dbar  <-  z*rowSums(post*b) # n by p
      Q.gam.theta <- matrix(0, nrow=p, ncol=2*m)  # p by 2*m matrix
      for (i in 1:m) {
        C.i <- array(C0[, i, 1, ], dim=c(n, 2))  # n by 2
        Q.i.1 <- colSums(tKR(-C.i+b[, i]*c0[, i, ], z*post[, i]))
        Q.i.2 <- colSums(tKR(c0[, i, ]*post[, i], dbar))  # p*2 vector
        Q.gam.theta[, (2*i-1):(2*i)] <- matrix(Q.i.1-Q.i.2, nrow=p, ncol=2)
      }

      Q.gam.theta <- Q.gam.theta[, c(p1, p2),  drop=FALSE]  # p by 2*m
      w1 <- (post*b)%*%t(a) - rowSums(post*b)*t(abar)  # n by m-1
      Q.pi.gam.0 <- colSums(tKR(w1, z))  # (m-1)*p vector
      Q.pi.gam  <- matrix(Q.pi.gam.0, nrow=m-1, ncol=p)
      Q.gam     <- - t(z) %*% (z*rowSums(post*B)) -
        matrix(colSums(tKR(dbar, dbar)), nrow=p, ncol=p)
      I <- cbind(I, -rbind(Q.pi.gam, t(Q.gam.theta)))
      I <- rbind(I, -cbind(t(Q.pi.gam), Q.gam.theta, Q.gam))
    }  # end if (!is.null(z))

  } else  { # compute I with (vcov.method == "OPG")
    other.method = "Hessian"
    score <- t(abar)
    for (j in 1:m) { score <- cbind(score, c0[, j, ]*post[, j]) }

    ind <- c(c(1:(m-1)), p1+m-1, p2+m-1)
    score <- score[, ind]
    I <- t(score) %*% score

    if (!is.null(z)) {
      dbar  <-  z*rowSums(post*b) # n by p
      score <- cbind(score, dbar)
      I <- t(score) %*% score
    }

  } # end if (vcov.method == "OPG")

  vcov <- try(solve(I))
  if (class(vcov) == "try-error" || any(diag(vcov) <0) ) {
    vcov <- matrix(NaN, nrow = 3*m-1+p, ncol = 3*m-1+p)
    warning("Fisher information matrix is singular and/or the
            variance is estimated to be negative. Consider using vcov.method=\"",other.method,"\".")
  }

  # Because the variance is parameterized as sigma^2, we convert is to sigma

  c.mat.vec <- c(rep(1, m-1+m), (1/sigma^(1/2))/2, rep(1, p))
  vcov <- diag(c.mat.vec) %*% vcov %*% diag(c.mat.vec)

  # Add the variance of alpha_m
  len   <- length(coefficients)
  M.mat <- diag(len-1)
  M.mat <- rbind(M.mat[1:(m-1), ], c(rep(-1, m-1), rep(0, len-m)),
                 M.mat[m:(len-1), ])

  vcov     <- M.mat %*% vcov %*% t(M.mat)

}   # end else (i.e., m >= 2)

vcov

}  # end function normalpanelmixVcov
