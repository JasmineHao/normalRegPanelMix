#' Computes the variance-covariance matrix of the MLE of m-component normal mixture.
#' @export
#' @title regpanelmixVcov
#' @name regpanelmixVcov
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param coefficients (alpha_1, ..., alpha_m, mu_1, ..., mu_m, sigma_1, ..., sigma_m, gam)
#' @param z n by p matrix of regressor associated with gam
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
regpanelmixVcov <- function(y, x, coefficients, z = NULL, vcov.method = c("Hessian", "OPG")) {
  # Computes the variance-covariance matrix of the MLE of m-component normal regression mixture
  # Input
  #  y  : n by 1 vector of dependent variable
  #  x  : n by q1-1 matrix of regressor NOT including an intercept
  #  coefficients : (alpha_1,...,alpha_m,mubeta_1^T, ...,mubeta_m^T,sigma_1, ..., sigma_m,gam^T)
  #  z  : n by p matrix of regressor associated with gam
  # Output
  #  vcov: variance-covariance matrix
  t  <- nrow(y)
  n  <- ncol(y)
  y     <- as.vector(y)

  len   <- length(coefficients)
  p     <- 0
  gam <- NULL
  vcov.method <- match.arg(vcov.method)

  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    gam <- coefficients[(len-p+1):len]
  }

  x   <- as.matrix(x)
  x1  <- cbind(1, x)
  q1   <- ncol(x1)
  q   <- ncol(x)

  m  <- (len-p)/(3+q)
  if (round(m) != m)
    stop("The dimension of the coefficients is incompatible with x and z. Please check the data.")

  alpha   <- coefficients[1:m]  # m-vector
  mubeta  <- matrix(coefficients[(m+1):((2+q)*m)], nrow=q+1, ncol=m)  # q+1 by m
  sigma   <- coefficients[((2+q)*m+1):((3+q)*m)]  # m-vector

  if (m == 1) {
    xz1 <- cbind(x1,z)
    I <- matrix(0, nrow=q1+p+1, ncol=q1+p+1)
    I[1:(q1+p), 1:(q1+p)] <- t(xz1) %*% xz1/sigma^2
    I[(q1+p+1), (q1+p+1)] <- n*2/sigma^2
    if (p != 0){
      s.1 <- c(1:q1, (q1+p+1), (q1+1):(q1+p))
      I <- I[s.1, s.1]
    }

    vcov <- solve(I)
    # Because the variance is parameterized as sigma^2, we convert it to sigma
    c.mat.vec <- c(rep(1,q1),(1/sigma^(1/2))/2,rep(1,p))
    vcov <- diag(c.mat.vec) %*% vcov %*% diag(c.mat.vec)

  } else {  # end of if (m == 1)
    # m >= 2
    # Compute posterior probabilities, and adjust y if z is present
    sigma0  <- rep(1, m)  # dummy
    mu0     <- double(m)  # dummy
    an      <- 1/n  # penalty term for variance
    h       <- 0
    tau     <- 0.5
    k <- 0
    epsilon <- 1e-08
    maxit = 2
    ninits = 1
    b <- matrix( rep( coefficients, ninits), ncol = ninits)

    if (is.null(z))
      out.p <- cppRegPanelmixPMLE(b, y, x, matrix(0), mu0, sigma0, m, p, t, an, maxit, ninits, epsilon, tau, h, k)
    else
    {
      out.p <- cppRegPanelmixPMLE(b, y, x, z, mu0, sigma0, m, p, t, an, maxit, ninits, epsilon, tau, h, k)
      # Adjust y
      y <- as.vector(y - z %*% gam)
    }

    post <- matrix(out.p$post, nrow=n)

    p2 <- seq(q1+1, (q1+1)*m, by=q1+1)  # sequence of q1+1, (q1+1)*2, ... , (q1+1)*m
    p1 <- (1:((q1+1)*m))[-p2]        # other values from 1, ..., (q1+1)*m

    a <- diag(1/alpha[-m], nrow=m-1, ncol=m-1)
    a <- cbind(a, -1/alpha[m])  # m-1 by m matrix of a_i's
    abar <- a %*% t(post)  # m-1 by n

    xtheta <- x1 %*% mubeta  # n by m

    Z0 <- t(t(y-xtheta)/sigma)          # normalized data, n by m
    f <- t(t(exp(-Z0^2/2)/sqrt(2*pi))/sigma)  # pdf, n by m
    phi <- t(t(f)*alpha)            # n by m
    f0 <- rowSums(phi)              # data pdf, n by 1

    vinv <- 1/(sigma*sigma)  # m-vector

    b <- t(t(Z0)/sigma)  # n by m
    B <- t(vinv - t(b*b))  # n by m

    c0 <- array(0,dim=c(n, m, q1+1))
    c0[, , (1:q1)] <- array(tKR(x1, b), dim=c(n, m, q1))
    c0[, , (q1+1)] <- -B/2

    # Compute Hessian-based I
    if (vcov.method == "Hessian") {
      other.method = "OPG"
      C0 <- array(0, dim=c(n, m, q1+1, q1+1))
      x11 <- array(tKR(x1, x1), dim = c(n, q1, q1))
      for (i in 1:m) {
        C0[, i, (1:q1), (1:q1)] <- x11*vinv[i]
        C0[, i, (1:q1), q1+1]   <- C0[, i, q1+1, (1:q1)] <- x1*b[, i]*vinv[i] # n by q1
      }
      C0[, , q1+1, q1+1] <- t((vinv - 2*t(B))*vinv)/2      # n by m

      Q.pi <- - abar %*% t(abar)  # m-1 by m-1

      Q.pi.theta <- matrix(0,nrow=m-1,ncol=(q1+1)*m)  # m-1 by (q1+1)*m
      for (i in 1:m) {
        zi <- a[, i] - abar  # m-1 by n
        wi <- c0[, i, ]*post[, i]  # n by q1+1
        Q.i <- colSums(tKR(wi, t(zi)))  # (q1+1)*(m-1) vector
        # first q1*(m-1) elements correspond to mubeta x pi,
        # last m-1 elements correspond to sigma x pi,
        Q.pi.theta[,(q1*(i-1)+1):(q1*i)] <- matrix(Q.i[1:(q1*(m-1))],ncol=q1)  # m-1 by q1 matrix
        Q.pi.theta[, q1*m+i] <- Q.i[(q1*(m-1)+1):((q1+1)*(m-1))]  # m-1 vector
      }

      Q.theta <- matrix(0, nrow=(q1+1)*m, ncol=(q1+1)*m)
      for (i in 2:m) {  # off-diagonal blocks
        for (j in 1:(i-1)) {
          wi  <- c0[, i, ]*post[, i] # n by q1+1
          wj  <- c0[, j, ]*post[, j] # n by q1+1
          Q.ij <- - colSums(tKR(wi, wj))  # (q1+1)*(q1+1) vector
          Q.theta[((q1+1)*(i-1)+1):((q1+1)*i), ((q1+1)*(j-1)+1):((q1+1)*j)] = t(matrix(Q.ij, nrow=q1+1, ncol=q1+1))
        }
      }

      Q.theta <- Q.theta + t(Q.theta)
      for (i in 1:m) {  # diagonal blocks
        C.ii   <- array(C0[, i, , ], dim=c(n, q1+1, q1+1))
        Q.ii.1   <- apply(C.ii*post[,i], c(2, 3), sum)
        w.ii   <- tKR(c0[, i, ], c0[, i, ])*post[, i]*(1-post[, i])
        Q.ii.2   <- matrix(colSums(w.ii), nrow=q1+1, ncol=q1+1)
        Q.theta[((q1+1)*(i-1)+1):((q1+1)*i), ((q1+1)*(i-1)+1):((q1+1)*i)] <- -Q.ii.1 + Q.ii.2
      }

      # q1+1,2*(q1+1),...,m*(q1+1)th rows and columns = sigma
      # other rows and columns = mubeta
      Q.theta <- Q.theta[c(p1, p2), c(p1, p2)]  # first block = wrt mubeta, second blosk = wrt sigma

      dimI <- m-1+(q1+1)*m
      I <- matrix(0, nrow=dimI, ncol=dimI)
      I[1:(m-1), 1:(m-1)] <- - Q.pi
      I[1:(m-1), m:dimI]  <- - Q.pi.theta
      I[m:dimI, 1:(m-1)]  <- - t(Q.pi.theta)
      I[m:dimI, m:dimI]   <- - Q.theta

      if (!is.null(z)) {
        dbar <-  z*rowSums(post*b)  # n by p
        Q.gam.theta <- matrix(0, nrow=p, ncol=(q1+1)*m)  # p by (q1+1)*m matrix
        for (i in 1:m) {
          C.i <- array(C0[, i, 1, ], dim=c(n, q1+1))  # n by q1+1
          Q.i.1 <- colSums(tKR(-C.i+b[, i]*c0[, i, ], z*post[, i])) # p*(q1+1) vector
          Q.i.2 <- colSums(tKR(c0[, i, ]*post[, i], dbar))  # p*(q1+1) vector
          Q.gam.theta[, ((q1+1)*(i-1)+1):((q1+1)*i)] <- matrix(Q.i.1+Q.i.2, nrow=p, ncol=q1+1)
        }

        Q.gam.theta <- Q.gam.theta[, c(p1, p2), drop=FALSE]  # p by (q1+1)*m
        w1 <- (post*b)%*%t(a) - rowSums(post*b)*t(abar)  # n by m-1
        Q.pi.gam.0 <- colSums(tKR(w1, z))  # (m-1)*p vector
        Q.pi.gam  <- matrix(Q.pi.gam.0, nrow=m-1, ncol=p)
        Q.gam     <- - t(z)%*%(z*rowSums(post*B)) -
          matrix(colSums(tKR(dbar, dbar)), nrow=p, ncol=p)

        I <- cbind(I, -rbind(Q.pi.gam, t(Q.gam.theta)))
        I <- rbind(I, -cbind(t(Q.pi.gam), Q.gam.theta, Q.gam))
      }  # end if (!is.null(z))

    }  else {  # compute I with (method == "OPG")
      other.method = "Hessian"
      c0.a <- array(0, dim=c(n, m, 2))
      c0.a[, , 1] <- b  # n by m
      c0.a[, , 2] <- -B/2  # n by m

      score <- t(abar)

      for (j in 1:m) {
        # score.o <- cbind(score.o, c0[, j, ]*post[, j])
        score <- cbind(score, x1*c0.a[, j, 1]*post[, j], c0.a[, j, 2]*post[, j])
        # print(all.equal(score.o, score))
      }

      ind <- c(1:(m-1), p1+m-1, p2+m-1)
      score <- score[, ind]
      I <- t(score) %*% score

      if (!is.null(z))  {
        dbar <-  z*rowSums(post*b)  # n by p
        score <- cbind(score, dbar)
        I <- t(score) %*% score
      }

    }  # end if (method=="OPG")

    vcov <- try(solve(I))
    if (class(vcov) == "try-error" || any(diag(vcov) <0) ) {
      vcov <- matrix(NaN, nrow = (2+q1)*m-1+p, ncol = (2+q1)*m-1+p)
      warning("Fisher information matrix is singular and/or the
              variance is estimated to be negative. Consider using vcov.method=\"",other.method,"\".")
    }

    # Because the variance is parameterized as sigma^2, we convert it to sigma

    c.mat.vec <- c(rep(1, m-1+m*q1), (1/sigma^(1/2))/2, rep(1, p))
    vcov <- diag(c.mat.vec) %*% vcov %*% diag(c.mat.vec)
    # vcov.opg <- diag(c.mat.vec) %*% vcov.opg %*% diag(c.mat.vec)

    # Add the variance of alpha_m
    M.mat <- diag(len-1)
    M.mat <- rbind(M.mat[1:(m-1),], c(rep(-1,m-1),rep(0,len-m)), M.mat[m:(len-1),])

    vcov <- M.mat %*% vcov %*% t(M.mat)
    # vcov.opg <- M.mat %*% vcov.opg %*% t(M.mat)

  }   # end else (i.e., m >= 2)

  vcov

}  # end function regpanelmixVcov
