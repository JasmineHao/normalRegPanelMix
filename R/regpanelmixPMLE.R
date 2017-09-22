
#' @description Estimates parameters of a finite mixture of univariate normals by
#' the method of penalized maximum likelhood. Using this function is equivalent to
#' calling normalmixPMLE with regressors specified by x as a parameter.
#' @export
#' @title regpanelmixPMLE
#' @name regpanelmixPMLE
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param m The number of components in the mixture
#' @param z n by p matrix of regressor associated with gam
#' @param vcov.method Method used to compute the variance-covariance matrix, one of \code{"Hessian"} and \code{"OPG"}.
#' The default option is \code{"Hessian"}. When \code{method = "Hessian"}, the variance-covarince matrix is
#' estimated by the Hessian using the formula given in Boldea and Magnus (2009).
#' When \code{method = "OPG"}, the outer product of gradients is used.
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit The maximum number of iterations.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param binit The initial value of parameter vector that is included as a candidate parameter vector
#' @return  A list of class \code{normalMix} with items:
#' \item{coefficients}{A vector of parameter estimates. Ordered as \eqn{\alpha_1,\ldots,\alpha_m,\mu_1,\ldots,\mu_m,\sigma_1,\ldots,\sigma_m,\gam}.}
#' \item{parlist}{The parameter estimates as a list containing alpha, mu, and sigma (and gam if z is included in the model).}
#' \item{vcov}{The estimated variance-covariance matrix.}
#' \item{loglik}{The maximized value of the log-likelihood.}
#' \item{penloglik}{The maximized value of the penalized log-likelihood.}
#' \item{aic}{Akaike Information Criterion of the fitted model.}
#' \item{bic}{Bayesian Information Criterion of the fitted model.}
#' \item{postprobs}{n by m matrix of posterior probabilities for observations}
#' \item{components}{n by 1 vector of integers that indicates the indices of components
#' each observation belongs to based on computed posterior probabilities}
#' \item{call}{The matched call.}
#' \item{m}{The number of components in the mixture.}
#' @note \code{regpanelmixPMLE} maximizes the penalized log-likelihood function
#' using the EM algorithm with combining short and long runs of EM steps as in Biernacki et al. (2003).
#' \code{regpanelmixPMLE} first runs the EM algorithm from \code{ninits}\eqn{* 4m(1 + p)} initial values
#' with the convertence criterion \code{epsilon.short} and \code{maxit.short}.
#' Then, \code{regpanelmixPMLE} uses \code{ninits} best initial values to run the EM algorithm
#' with the convertence criterion \code{epsilon} and \code{maxit}.
#' @references     Biernacki, C., Celeux, G. and Govaert, G. (2003)
#' Choosing Starting Values for the EM Algorithm for Getting the
#' Highest Likelihood in Multivariate Gaussian Mixture Models,
#' \emph{Computational Statistics and Data Analysis}, \bold{41}, 561--575.
#'
#' Boldea, O. and Magnus, J. R. (2009)
#' Maximum Likelihood Estimation of the Multivariate Normal Mixture Model,
#' \emph{Journal of the American Statistical Association},
#' \bold{104}, 1539--1549.
#'
#' Chen, J., Tan, X. and Zhang, R. (2008)
#' Inference for Normal Mixtures in Mean and Variance,
#' \emph{Statistica Sinica}, \bold{18}, 443--465.
#'
#' McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley \& Sons, Inc.
#' @examples
#' data(faithful)
#' attach(faithful)
#' regpanelmixPMLE(y = eruptions, x = waiting, m = 1)
#' regpanelmixPMLE(y = eruptions, x = waiting, m = 2)
regpanelmixPMLE <- function (y, x, m = 2, z = NULL, vcov.method = c("Hessian", "OPG", "none"),
                        ninits = 10, epsilon = 1e-08, maxit = 2000,
                        epsilon.short = 1e-02, maxit.short = 500, binit = NULL) {
  if (is.null(x)){
    return(normalpanelmixPMLE(y = y, x = x, m = m, z = z, vcov.method = vcov.method,
           ninits= ninits, epsilon = epsilon, maxit = maxit,
           epsilon.short = epsilon.short, maxit.short = maxit.short,
           binit = binit))
  }
  t  <- nrow(y)
  n  <- ncol(y)
  nt <- n*t
  
  y   <- as.vector(y)
  x   <- as.matrix(x)   # n by (q1-1) matrix
  # n   <- length(y)
  if (nrow(x) != nt) { stop("y and x must have the same number of rows.") }
  x1  <- cbind(1, x)
  q1   <- ncol(x1)

  p       <- 0
  gam   <- NULL
  ninits.short <- ninits*10*(q1+p)*m
  vcov.method <- match.arg(vcov.method)

  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    if (nrow(z) != n) { stop("y and z must have the same number of rows.") }
  }

  npar    <- m-1 + (q1+1)*m + p  # number of parameters
  xz      <- cbind(x, z)
  ls.out  <- lsfit(xz, y)
  sd0     <- sqrt(mean(ls.out$residuals^2))

  if (m == 1) {
    mubeta <- as.matrix(unname(ls.out$coeff[1:q1]))
    if (!is.null(z)) {gam <- unname(ls.out$coeff[(q1+1):(q1+p)])}
    res     <- ls.out$residuals
    sigma   <- sqrt(mean(res*res))
    loglik  <- - (n/2)*(1 + log(2*pi) + 2*log(sigma))
    aic     <- -2*loglik + 2*npar
    bic     <- -2*loglik + log(n)*npar
    penloglik <- loglik

    parlist <- list(alpha = 1, mubeta = mubeta, sigma = sigma, gam = gam)
    coefficients <- c(alpha = 1, mubeta = mubeta, sigma = sigma, gam = gam)
    postprobs <- rep(1, n)

  } else {  # m >= 2

    # generate initial values
    tmp <- regpanelmixPMLEinit(y = y, x = x, z = z, ninits = ninits.short, m = m)

    h       <- 0    # setting h=0 gives PMLE
    tau     <- 0.5  # setting tau=0.5 gives PMLE
    k <- 0 # setting k=0 gives PMLE

    sigma0  <- rep(sd0, m)
    mu0     <- double(m)    # dummy
    an      <- 1/n  # penalty term for variance

    if (is.null(z))
      ztilde <- matrix(0) # dummy
    else
      ztilde <- z

    # short EM
    b0 <- rbind( tmp$alpha, tmp$mubeta, tmp$sigma, tmp$gam )
    if (!is.null(binit)) {
      b0[ , 1] <- binit
    }
    out.short <- cppRegPanelmixPMLE(b0, y, x, ztilde, mu0, sigma0, m, p, t, an, maxit.short,
                               ninits.short, epsilon.short)
    # long EM
    components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
    b1 <- b0[ ,components] # b0 has been updated
    out <- cppRegPanelmixPMLE(b1, y, x, ztilde, mu0, sigma0, m, p, t, an, maxit, ninits, epsilon)
    index     <- which.max(out$penloglikset)
    alpha <- b1[1:m,index] # b0 has been updated

    mubeta <- matrix(b1[(1+m):((q1+1)*m),index],nrow=q1,ncol=m)
    # mubeta <- matrix(b1[3:10,index],nrow=q1,ncol=m)
    sigma <- b1[(1+(q1+1)*m):((q1+2)*m),index]
    if (!is.null(z)) {
      gam     <- b1[((q1+2)*m+1):((q1+2)*m+p),index]
    }
    penloglik <- out$penloglikset[index]
    loglik    <- out$loglikset[index]

    postprobs <- matrix(out$post[,index], nrow=n)

    aic <- -2*loglik + 2*npar
    bic <- -2*loglik + log(n)*npar

    mu.order  <- order(mubeta[1,])
    alpha     <- alpha[mu.order]
    mubeta    <- mubeta[,mu.order]
    sigma     <- sigma[mu.order]

    postprobs <- postprobs[, mu.order]
    colnames(postprobs) <- c(paste("comp", ".", 1:m, sep = ""))

    mubeta.name <- matrix(0,nrow = q1, ncol = m)
    mubeta.name[1,] <- paste("mu", 1:m, sep = "")

    if (q1 == 2) {
      mubeta.name[2,] <- paste("beta", 1:m,  sep = "")
    } else {
      for (i in 1:(q1-1)) {
        for (j in 1:m) {
          mubeta.name[i+1,j] <- paste("beta", j, i, sep = "")
        }
      }
    }

    parlist <- list(alpha = alpha, mubeta = mubeta, sigma = sigma, gam = gam)
    coefficients <- unlist(parlist)
    names(coefficients)[(m+1):((q1+1)*m)] <- c(mubeta.name)
  }  # end m >= 2

  if (vcov.method == "none") {
    vcov <- NULL
  } else {
    vcov <- regpanelmixVcov(y = y, x = x, coefficients = coefficients, z = z , vcov.method = vcov.method)
  }

  a <- list(coefficients = coefficients, parlist = parlist, vcov = vcov, loglik = loglik,
            penloglik = penloglik, aic = aic, bic = bic, postprobs = postprobs,
            components = getComponentcomponents(postprobs),
            call = match.call(), m = m, label = "PMLE")

  class(a) <- "normalregMix"

  a

}  # end function regpanelmixPMLE
