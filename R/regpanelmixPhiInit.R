#' @description Generates lists of parameters for initial candidates used by
#' the modified EM test for mixture of normals.
#' @export
#' @title regpanelmixPhiInit
#' @name regpanelmixPhiInit
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gam
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gam_1, ..., gam_m))
#' @param z n by p matrix of regressor associated with gam
#' @param h h used as index for pivoting
#' @param tau Tau used to split the h-th component
#' @param ninits number of initial values to be generated
#' @return A list with the following items:
#' \item{alpha}{m+1 by ninits matrix for alpha}
#' \item{mubeta}{q+1 by m+1 by ninits array for mu and beta}
#' \item{sigma}{m+1 by ninits matrix for sigma}
#' \item{gam}{m+1 by ninits matrix for gam}
regpanelmixPhiInit <- function (y, x, z = NULL, parlist, h, tau, ninits = 1)
{
  # if (normalregpanelmix.test.on) # initial values controlled by normalregpanelmix.test.on
  #   set.seed(normalregpanelmix.test.seed)
  
  y  <- as.vector(y)
  n   <- length(y)
  x   <- matrix(x,nrow=n)
  x1   <- cbind(1,x)
  q1   <- ncol(x1)
  
  alpha0   <- parlist$alpha
  mubeta0  <- parlist$mubeta
  sigma0   <- parlist$sigma
  mu0    <- mubeta0[1,]
  beta0  <- mubeta0[-1,,drop=FALSE]
  gam  <- NULL
  # print(parlist)
  if (!is.null(z)) {
    gam0  <- parlist$gam
    p     <- ncol(z)
    y     <- as.vector(y - z %*% gam0)
    gam1   <- c(rep(1,p),runif(p*(ninits-1),min=-2,max=2))*gam0
    gam   <- matrix(gam1,nrow=p)
  }

  minMU  <- min(y - x%*%beta0)
  maxMU  <- max(y - x%*%beta0)
  
  m     <- length(alpha0)
  
  if (m>=2) {
    mid <- (mu0[1:(m-1)]+mu0[2:m])/2  # m-1 by 1
    lb0 <- c(minMU,mid)          # m by 1
    lb  <- c(lb0[1:h],lb0[h:m])      # m+1 by 1
    ub0 <- c(mid,maxMU)          # m by 1
    ub  <- c(ub0[1:h],ub0[h:m])      # m+1 by 1
    beta.hyp <- cbind(beta0[,1:h,drop=FALSE],beta0[,h:m,drop=FALSE])    # q1-1 by m+1
    mu.hyp <- c(mu0[1:h],mu0[h:m])  # m+1 by 1 vector
  } else {  # m=1
    lb  <- c(minMU,minMU)
    ub  <- c(maxMU,maxMU)
    beta.hyp <- cbind(beta0,beta0)   # q1-1 by 2 matirx
    mu.hyp <- c(mu0,mu0)    # 2 by 1 vector
  }
  
  mubeta <- matrix(0, nrow=q1*(m+1), ncol=ninits)
  for (j in 1:(m+1)) {
    mubeta[(q1*(j-1)+1), ] <- runif(ninits, min=lb, max=ub)
    mubeta[(q1*(j-1)+1), 1] <- mu.hyp[j]
    for (i in 2:q1) {
      mubeta[(q1*(j-1)+i), ] <- beta.hyp[i]*runif(ninits, min=-2, max=2)
      mubeta[(q1*(j-1)+i), 1] <- beta.hyp[i]
    }
  }
  
  sigma.hyp <- c(sigma0[1:h],sigma0[h:m])  # m+1 by 1
  sigma1 <- c(rep(1,m+1),runif((m+1)*(ninits-1),min=0.25,max=2))*sigma.hyp
  sigma <- matrix(sigma1,nrow=m+1)
  
  alpha.hyp <- c(alpha0[1:h],alpha0[h:m])  # m+1 by 1
  alpha.hyp[h:(h+1)] <- c(alpha.hyp[h]*tau,alpha.hyp[h+1]*(1-tau))
  alpha <- matrix(rep.int(alpha.hyp,ninits),nrow=m+1)
  
  list(alpha = alpha, mubeta = mubeta, sigma = sigma, gam = gam)
  
}  # end function regpanelmixPhiInit
