#' Generate initial values used by the PMLE of mixture of normals
#' @export
#' @title regpanelmixPMLEinit
#' @name regpanelmixPMLEinit
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param z n by p matrix of regressor associated with gam
#' @param ninits number of initial values to be generated
#' @param m The number of components in the mixture
#' @return A list with the following items:
#' \item{alpha}{m by ninits matrix for alpha}
#' \item{mubeta}{q+1 by m by ninits array for mu and beta}
#' \item{sigma}{m by ninits matrix for sigma}
#' \item{gam}{m by ninits matrix for gam}
regpanelmixPMLEinit <- function (y, x, z = NULL, ninits = 1, m = 2)
{
  # if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
  #   set.seed(normalregMix.test.seed)
  t  <- nrow(y)
  n  <- ncol(y)
  y <- as.vector(y)
  q1  <- ncol(x)+1
  p  <- ncol(z)

  gam <- NULL
  if (!is.null(z)) {
    out     <- lsfit(cbind(x, z), y)
    gam0  <- out$coef[(q1+1):(q1+p)]
    gam   <- matrix(runif(p*ninits, min=0.5, max=1.5), nrow=p)*gam0
    mubeta_hat <- out$coef[1:q1]
    y     <- y - z %*% gam0
    r     <- out$residuals
    stdR  <- sd(r)
  } else {
    out         <- lsfit(x, y)
    mubeta_hat  <- out$coef
    r           <- out$residuals
    stdR        <- sd(r)
  }

  alpha <- matrix(runif(m*ninits), nrow=m)
  alpha <- t(t(alpha)/colSums(alpha))


  minMU <- min(y - x %*% mubeta_hat[-1])
  maxMU <- max(y - x %*% mubeta_hat[-1])
  mubeta <- matrix(0, nrow=q1*m, ncol=ninits)
  for (j in 1:m) {
    mubeta[(q1*(j-1)+1), ] <- runif(ninits, min=minMU, max=maxMU)
    for (i in 2:q1) {
      mubeta[(q1*(j-1)+i), ] <- mubeta_hat[i]*runif(ninits, min=-2, max=2)
    }
  }
  sigma <- matrix(runif(m*ninits, min=0.01, max=1), nrow=m)*stdR

  list(alpha = alpha, mubeta = mubeta, sigma = sigma, gam = gam)

}  # end function regpanelmixPMLEinit
