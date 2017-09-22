#' @description Given a pair of h and tau and data, compute ordinary &
#' penalized log-likelihood ratio resulting from MEM algorithm at k=1,2,3,
#' tailored for parallelization.
#' @export
#' @title normalpanelmixMaxPhiStep
#' @name normalpanelmixMaxPhiStep
#' @param htaupair A set of h and tau
#' @param y n by 1 vector of data
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param p Dimension of z
#' @param an a term used for penalty function
#' @param ninits The number of randomly drawn initial values.
#' @param ninits.short The number of candidates used to generate an initial phi, in short MEM
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @param verb Determines whether to print a message if an error occurs.
#' @return A list of phi, log-likelihood, and penalized log-likelihood resulting from MEM algorithm.
normalpanelmixMaxPhiStep <- function (htaupair, y, parlist, z = NULL, p,
                                 an,
                                 ninits, ninits.short,
                                 epsilon.short, epsilon,
                                 maxit.short, maxit,
                                 verb,update.alpha)
{
  alpha0 <- parlist$alpha

  m      <- length(alpha0)
  m1     <- m+1
  k      <- 1
  t  <- nrow(y)
  n  <- ncol(y)
  h      <- as.numeric(htaupair[1])
  tau    <- as.numeric(htaupair[2])

  mu0    <- parlist$mu
  mu0h   <- c(-1e+10,mu0,1e+10)        # m+2 by 1
  sigma0 <- parlist$sigma
  sigma0h<- c(sigma0[1:h],sigma0[h:m]) # m+1 by 1
  gam0 <- parlist$gam
  if (is.null(z)) {
    ztilde <- matrix(0) # dummy
    gam <- NULL
  }else{
    ztilde <- as.matrix(z)
  }
  # generate initial values
  tmp <- normalpanelmixPhiInit(y = y, parlist = parlist, z = z, h=h, tau = tau, ninits = ninits.short)

  # short EM
  b0 <- as.matrix(rbind(tmp$alpha, tmp$mu, tmp$sigma, tmp$gam))
  out.short <- cppnormalpanelmixPMLE(b0, as.vector(y), ztilde, mu0h, sigma0h, m1, p, t, an, maxit.short, ninits.short,
                                epsilon.short, tau, h, k, update.alpha)
  components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
  if (verb && any(out.short$notcg)) {
    cat(sprintf("non-convergence rate at short-EM = %.3f\n",mean(out.short$notcg)))
  }
  # long EM
  b1 <- as.matrix(b0[ ,components])
  out <- cppnormalpanelmixPMLE(b1, as.vector(y), ztilde, mu0h, sigma0h, m1, p, t, an, maxit, ninits, epsilon, tau, h, k, update_alpha = update.alpha)

  index     <- which.max(out$penloglikset)
  alpha <- b1[1:m1,index]
  mu <- b1[(1+m1):(2*m1),index]
  sigma <- b1[(1+2*m1):(3*m1),index]
  if (!is.null(z)) {
    gam     <- b1[(3*m1+1):(3*m1+p),index]
  }
  mu.order  <- order(mu)
  alpha     <- alpha[mu.order]
  mu        <- mu[mu.order]
  sigma     <- sigma[mu.order]
  sigma0h <- sigma0h[mu.order]
  b <- as.matrix( c(alpha, mu, sigma, gam) )

  # initilization
  loglik <-  vector("double", 3)
  penloglik <-  vector("double", 3)
  coefficient <- vector("double", length(b))

  penloglik[1] <- out$penloglikset[[index]]
  loglik[1]    <- out$loglikset[[index]]
  for (k in 2:3) {
    ninits <- 1
    maxit <- 2
    # Two EM steps
    out <- cppnormalpanelmixPMLE(b, as.vector(y), ztilde, mu0h, sigma0h, m1, p,t, an, maxit, ninits, epsilon, tau, h, k, update_alpha = update.alpha)
    alpha <- b[1:m1,1] # b has been updated
    mu <- b[(1+m1):(2*m1),1]
    sigma <- b[(1+2*m1):(3*m1),1]
    if (!is.null(z)) {
      gam     <- b[(3*m1+1):(3*m1+p),1]
    }
    loglik[k]    <- out$loglikset[[1]]
    penloglik[k]   <- out$penloglikset[[1]]

    # Check singularity: if singular, break from the loop
    if ( any(sigma < 1e-06) || any(alpha < 1e-06) || is.na(sum(alpha)) ) {
      loglik[k]    <- -Inf
      penloglik[k]   <- -Inf
      break
    }

    mu.order  <- order(mu)
    alpha     <- alpha[mu.order]
    mu        <- mu[mu.order]
    sigma     <- sigma[mu.order]
    sigma0h <- sigma0h[mu.order]
  }
  coefficient <- as.matrix( c(alpha, mu, sigma, gam) ) # at k=3

  return (list(coefficient = coefficient, loglik = loglik, penloglik = penloglik))
}
