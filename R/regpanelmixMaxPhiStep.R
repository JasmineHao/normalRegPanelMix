#' @description Given a pair of h and tau and data, compute ordinary &
#' penalized log-likelihood ratio resulting from MEM algorithm at k=1,2,3, 
#' tailored for parallelization.
#' @export
#' @title regpanelmixPhiStep
#' @name regpanelmixPhiStep
#' @param htaupair A set of h and tau
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gam
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gam_1, ..., gam_m))
#' @param z n by p matrix of regressor associated with gam
#' @param p Dimension of z
#' @param an a term used for penalty function
#' @param ninits The number of randomly drawn initial values.
#' @param ninits.short The number of candidates used to generate an initial phi, in short MEM
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @param verb Determines whether to print a message if an error occurs.
#' @return A list of coefficients, log-likelihood, and penalized log-likelihood resulting from MEM algorithm.
regpanelmixPhiStep <- function (htaupair, y, x, parlist, z = NULL, p,
                           an,
                           ninits, ninits.short,
                           epsilon.short, epsilon,
                           maxit.short, maxit,
                           verb,update.alpha)
{
  alpha0 <- parlist$alpha
  
  x1    <- cbind(1,x)
  q1     <- ncol(x1)
  m      <- length(alpha0)
  m1     <- m+1
  t  <- nrow(y)
  n  <- ncol(y)
  nt <- n*t
  y   <- as.vector(y)
  x   <- as.matrix(x) 
  # n      <- length(y)
  h      <- as.numeric(htaupair[1])
  tau    <- as.numeric(htaupair[2])
  k <- 1
  
  mubeta0 <- parlist$mubeta
  mu0h <- c(-1e+10,mubeta0[1,],1e+10)        # m+2 by 1
  sigma0 <- parlist$sigma
  sigma0h <- c(sigma0[1:h],sigma0[h:m])  # m+1 by 1
  gam0 <- parlist$gam
  
  if (is.null(z)) {
    ztilde <- matrix(0) # dummy
    gam <- NULL
  } else{
    ztilde <- as.matrix(z)
  }
  # generate initial values
  # print(parlist)
  tmp <- regpanelmixPhiInit(y = y, x = x, z = z, parlist=parlist, h=h, tau, ninits = ninits.short)
  
  # short EM
  b0 <- as.matrix(rbind( tmp$alpha, tmp$mubeta, tmp$sigma, tmp$gam ))
  out.short <- cppRegPanelmixPMLE(b0, y, x, ztilde, mu0h, sigma0h, m1, p, t, an, maxit.short,
                             ninits.short, epsilon.short, tau, h, k)
  # long EM
  components <- order(out.short$penloglikset, decreasing = TRUE)[1:ninits]
  b1 <- as.matrix(b0[ ,components]) # b0 has been updated
  out <- cppRegPanelmixPMLE(b1, y, x, ztilde, mu0h, sigma0h, m1, p, t, an, maxit, ninits, epsilon, tau, h, k, update_alpha = update.alpha)
  
  index     <- which.max(out$penloglikset)
  alpha <- b1[1:m1,index] # b0 has been updated
  mubeta <- matrix(b1[(1+m1):((q1+1)*m1),index],nrow=q1,ncol=m1)
  sigma <- b1[(1+(q1+1)*m1):((q1+2)*m1),index]
  if (!is.null(z)) {
    gam     <- b1[((q1+2)*m1+1):((q1+2)*m1+p),index]
  }
  mu.order  <- order(mubeta[1,])
  alpha     <- alpha[mu.order]
  mubeta    <- mubeta[ ,mu.order]
  sigma     <- sigma[mu.order]
  sigma0h <- sigma0h[mu.order]
  b <- as.matrix(c(alpha,as.vector(mubeta),sigma,gam))
  
  # initilization
  penloglik <-  vector("double", 3)
  loglik <-  vector("double", 3)
  coefficient <- vector("double", length(b))
  
  penloglik[1] <- out$penloglikset[[index]]
  loglik[1]    <- out$loglikset[[index]]
  for (k in 2:3) {
    ninits <- 1
    maxit <- 2
    # Two EM steps
    out <- cppRegPanelmixPMLE(b, y, x, ztilde, mu0h, sigma0h, m1, p, t, an, maxit, ninits, epsilon, tau, h, k,update_alpha = update.alpha)
    alpha <- b[1:m1,1] # b0 has been updated
    mubeta <- matrix(b[(1+m1):((q1+1)*m1),1],nrow=q1,ncol=m1)
    sigma <- b[(1+(q1+1)*m1):((q1+2)*m1),1]
    if (!is.null(z)) {
      gam     <- b[((q1+2)*m1+1):((q1+2)*m1+p),1]
    }
    loglik[k]    <- out$loglik[[1]]
    penloglik[k]   <- out$penloglik[[1]]
    
    # Check singularity: if singular, break from the loop
    if ( any(sigma < 1e-06) || any(alpha < 1e-06) || is.na(sum(alpha)) ) {
      for (i in k:3) {
        loglik[k]    <- -Inf
        penloglik[k]   <- -Inf
      }
      break
    }
    mu.order  <- order(mubeta[1,])
    alpha     <- alpha[mu.order]
    mubeta    <- mubeta[ ,mu.order]
    sigma     <- sigma[mu.order]
    sigma0h <- sigma0h[mu.order]
  }
  coefficient <- as.matrix(c(alpha,as.vector(mubeta),sigma,gam)) # at k=3
  
  return (list(coefficient = coefficient, loglik = loglik, penloglik = penloglik))
}
