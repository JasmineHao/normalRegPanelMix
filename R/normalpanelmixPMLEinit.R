
#' @description Generate initial values used by the PMLE of mixture of normals
#' @title normalpanelmixPMLEinit
#' @name normalpanelmixPMLEinit
#' @param y n by 1 vector of data
#' @param z n by p matrix of regressor associated with gamma
#' @param ninits number of initial values to be generated
#' @param m The number of components in the mixture
#' @return A list with the following items:
#' \item{alpha}{m by ninits matrix for alpha}
#' \item{mu}{m by ninits matrix for mu}
#' \item{sigma}{m by ninits matrix for sigma}
#' \item{gam}{m by ninits matrix for gam}
normalpanelmixPMLEinit <- function (y, z = NULL, ninits = 1, m = 2)
{
  # if (normalregMix.test.on) # initial values controlled by normalregMix.test.on
  #   set.seed(normalregMix.test.seed)

  t  <- nrow(y)
  n  <- ncol(y)
  y  <- as.vector(y)

  p <- ncol(z)
  gam <- NULL
  if (!is.null(z)){
    out     <- lsfit(z,y)
    gam0  <- out$coef[-1]
    gam   <- matrix(runif(p*ninits, min=0.5, max=1.5), nrow=p) * gam0
    y       <- out$residuals + out$coef[1]
  }

  alpha <- matrix(runif(m * ninits), nrow=m)
  alpha <- t(t(alpha) / colSums(alpha))
  mu    <- matrix(runif(m * ninits, min = min(y), max = max(y)), nrow=m)
  sigma <- matrix(runif(m * ninits, min = 0.01, max = 2)*sd(y), nrow=m)

  list(alpha = alpha, mu = mu, sigma = sigma, gam = gam)

}  # end function normalpanelmixPMLEinit
