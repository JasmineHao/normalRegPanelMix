#' @description Generates lists of parameters for initial candidates used by
#' the modified EM test for mixture of normals.
#' @title normalpanelmixPhiInit
#' @name normalpanelmixPhiInit
#' @param y n by 1 vector of data
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma,
#' and gamma in the form of (alpha = (alpha_1, ..., alpha_m),
#' mu = (mu_1, ..., mu_m), sigma = (sigma_1, ..., sigma_m),
#' gam = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param h h used as index for pivoting
#' @param tau Tau used to split the h-th component
#' @param ninits number of initial values to be generated
#' @return A list with the following items:
#' \item{alpha}{m+1 by ninits matrix for alpha}
#' \item{mu}{m+1 by ninits matrix for mu}
#' \item{sigma}{m+1 by ninits matrix for sigma}
#' \item{gam}{m+1 by ninits matrix for gamma}
normalpanelmixPhiInit <- function (y, parlist, z = NULL, h, tau, ninits = 1)
{
  # if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
  #   set.seed(normalregMix.test.seed)

  t  <- nrow(y)
  n  <- ncol(y)
  p     <- ncol(z)

  mu0      <- parlist$mu
  sigma0   <- parlist$sigma
  alpha0   <- parlist$alpha
  m       <- length(alpha0)

  if (is.null(z))
    gam <- NULL
  else {
    gam0  <- parlist$gam
    y    <- y - z %*% gam0
    gam <- matrix(runif(p*ninits,min=0.5,max=1.5),nrow=p)*gam0
  }

  if (m>=2){
    mid <- (mu0[1:(m-1)]+mu0[2:m])/2  # m-1 by 1
    lb0 <- c(min(y),mid)        # m by 1
    lb  <- c(lb0[1:h],lb0[h:m])      # m+1 by 1
    ub0 <- c(mid,max(y))        # m by 1
    ub  <- c(ub0[1:h],ub0[h:m])      # m+1 by 1
  } else {
    lb  <- c(min(y),min(y))
    ub  <- c(max(y),max(y))
  }

  mu <- matrix(runif((m+1)*ninits,min=lb,max=ub),nrow=m+1)

  sigma.hyp <- c(sigma0[1:h],sigma0[h:m])  # m+1 by 1
  sigma <- matrix(runif((m+1)*ninits,min=sigma.hyp*0.25,max=sigma.hyp*2),nrow=m+1)

  alpha.hyp <- c(alpha0[1:h],alpha0[h:m])  # m+1 by 1
  alpha.hyp[h:(h+1)] <- c(alpha.hyp[h]*tau,alpha.hyp[h+1]*(1-tau))
  alpha <- matrix(rep.int(alpha.hyp,ninits),nrow=m+1)

  list(alpha = alpha, mu = mu, sigma = sigma, gam = gam)

}  # end function normalpanelmixPhiInit
