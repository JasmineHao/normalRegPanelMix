#' @description Estimates parameters of a finite panel mixture of univariate normals by
#' penalized maximum log-likelhood functions.
#' @export
#' @title normalpanelmixPMLE
#' @name normalpanelmixPMLE
#' @param y n by 1 vector of data
#' @param x n by q matrix of data for x (if exists)
#' @param m The number of components in the mixture
#' @param z n by p matrix of regressor associated with gamma
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
#' @return  A list of class \code{normalpanelmix} with items:
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
#' @note \code{normalpanelmixPMLE} maximizes the penalized log-likelihood function
#' using the EM algorithm with combining short and long runs of EM steps as in Biernacki et al. (2003).
#' \code{normalpanelmixPMLE} first runs the EM algorithm from \code{ninits}\eqn{* 4m(1 + p)} initial values
#' with the convertence criterion \code{epsilon.short} and \code{maxit.short}.
#' Then, \code{normalpanelmixPMLE} uses \code{ninits} best initial values to run the EM algorithm
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
#'
#' normalpanelmixPMLE(y = eruptions, m = 1)
#' normalpanelmixPMLE(y = eruptions, m = 2)
#'
#' out <- normalpanelmixPMLE(y = eruptions, m = 2)
#' summary(out)
normalpanelmixPMLE <- function (y, x = NULL, m = 2, z = NULL, vcov.method = c("Hessian", "OPG", "none"),
                           ninits = 25, epsilon = 1e-08, maxit = 2000,
                           epsilon.short = 1e-02, maxit.short = 500, binit = NULL) {
  t <- dim(y)[1] #Number of year
  n <- dim(y)[2] #Number of firms
  nt <- n * t
  p <- 0
  y <- as.vector(y)
  if (!is.null(x))
    #TODO(Jasmine): fix this
    return (regpanelmixPMLE(y = y, x = x, m = m, z = z, vcov.method = vcov.method,
                       ninits= ninits, epsilon = epsilon, maxit = maxit,
                       epsilon.short = epsilon.short, maxit.short = maxit.short,
                       binit = binit))

  gam <- NULL
  ninits.short <- ninits*10*(1 + p)*m
  vcov.method <- match.arg(vcov.method)

  if (!is.null(z)) {
    z     <- as.matrix(z)
    p     <- ncol(z)
    if (nrow(z) != nt) { stop("y and z must have the same number of rows.") }
    ls.out   <- lsfit(z, y)
    sd0    <- sqrt(mean(ls.out$residuals^2))
  } else {
    sd0   <- sd(y) * sqrt((n - 1) / n)
  }

  # if (m == 1) {
  #   if (!is.null(z)) {
  #     mu     <- unname(ls.out$coeff[1])
  #     gam <- unname(ls.out$coeff[2:(1 + p)])
  #     res    <- ls.out$residuals
  #     sigma <- sqrt(mean(res*res))
  #   } else {
  #     mu     <- mean(y)
  #     sigma  <- sd0
  #   }
  # 
  #   loglik   <- - (n/2) *(1 + log(2*pi) + 2*log(sigma))
  #   aic      <- -2*loglik + 2*(m-1+2*m+p)
  #   bic      <- -2*loglik + log(n)*(m-1+2*m+p)
  #   penloglik <- loglik
  # 
  #   parlist <- list(alpha = 1, mu = mu, sigma = sigma, gam = gam)
  #   coefficients <- c(alpha = 1, mu = mu, sigma = sigma, gam = gam)
  #   postprobs <- rep(1, n)
  # 
  # } else {  # m >= 2

    # generate initial values
    tmp <- normalpanelmixPMLEinit(y = y, z = z, ninits = ninits.short, m = m)

    # the following values for (h, k, tau, an) are given by default
    # h       <- 0  # setting h=0 gives PMLE
    # k       <- 0  # k is set to 0 because this is PMLE
    # tau     <- 0.5  # tau is set to 0.5 because this is PMLE
    an      <- 1/n  # penalty term for variance
    sigma0  <- rep(sd0, m)
    mu0     <- double(m+1) # dummy

    if (is.null(z)) {
      ztilde <- matrix(0) # dummy
    } else {
      ztilde <- z
    }
    # short EM
    b0 <- as.matrix(rbind( tmp$alpha, tmp$mu, tmp$sigma, tmp$gam ))
    if (!is.null(binit)) {
      b0[ , 1] <- binit
    }
    
    out.short <- cppnormalpanelmixPMLE(b0, y, ztilde, mu0, sigma0, m, p, t , an, maxit.short,
                                  ninits.short, epsilon.short)

    # long EM
    components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
    b1 <- b0[ ,components] # b0 has been updated
    out <- cppnormalpanelmixPMLE(b1, y, ztilde, mu0, sigma0, m, p,t, an, maxit, ninits, epsilon)

    index     <- which.max(out$penloglikset)
    alpha <- b1[1:m,index] # b0 has been updated
    mu <- b1[(1+m):(2*m),index]
    sigma <- b1[(1+2*m):(3*m),index]
    if (!is.null(z)) {
      gam     <- b1[(3*m+1):(3*m+p),index]
    }
    penloglik <- out$penloglikset[index]
    loglik    <- out$loglikset[index]
    postprobs <- matrix(out$post[,index], nrow=n)

    aic     <- -2*loglik + 2*(m-1+2*m+p)
    bic     <- -2*loglik + log(n)*(m-1+2*m+p)

    mu.order  <- order(mu)
    alpha     <- alpha[mu.order]
    mu        <- mu[mu.order]
    sigma     <- sigma[mu.order]

    postprobs <- postprobs[, mu.order]
    if (m > 1){
      colnames(postprobs) <- c(paste("comp", ".", 1:m, sep = ""))
    }
    parlist <- list(alpha = alpha, mubeta = mu, sigma = sigma, gam = gam)
    coefficients <- unlist(parlist)

  # } # end m >= 2

  if (vcov.method == "none") {
    vcov <- NULL
  } else {
    vcov <- normalpanelmixVcov(y = y, coefficients = coefficients, z = z , vcov.method = vcov.method)
  }

  a <- list(coefficients = coefficients, parlist = parlist, vcov = vcov, loglik = loglik,
            penloglik = penloglik, aic = aic, bic = bic, postprobs = postprobs,
            components = getComponentcomponents(postprobs),
            call = match.call(), m = m, label = "PMLE")

  class(a) <- "normalregMix"

  a
}  # end function normalpanelmixPMLE
