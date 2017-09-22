#' @description  Sequentially performs MEM test given the data for y and x
#' on the null hypothesis H_0: m = m_0 where m_0 is in {1, 2, ..., maxm}
#' @export
#' @title normalpanelmixMEMtestSeq
#' @name normalpanelmixMEMtestSeq
#' @param y n by 1 vector of data for y
#' @param x n by q matrix of data for x (if exists)
#' @param z n by p matrix of regressor associated with gamma
#' @param maxm The maximum number of components set as null hypothesis in the mixture
#' @param ninits The number of randomly drawn initial values.
#' @param maxit The maximum number of iterations.
#' @param nbtsp The number of bootstrap observations; by default, it is set to be 199
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' create a new one for computation accordingly.
#' @param crit.bootstrap.from The minimum m in null hypothesis to have critical values
#' calculated from bootstrap for the test statistics
#' @return A list of with the following items:
#' \item{alpha}{maxm by maxm matrix, whose i-th column is a vector of alphas estimated given the null hypothesis m_0 = i}
#' \item{mu}{maxm by maxm matrix, whose i-th column is a vector of mus estimated given the null hypothesis m_0 = i}
#' \item{sigma}{maxm by maxm matrix, whose i-th column is a vector of sigmas estimated given the null hypothesis m_0 = i}
#' \item{beta}{A list of length maxm, whose i-th element is a q times i matrix of betas estimated given the null hypothesis m_0 = i}
#' \item{gam}{maxm by maxm matrix, whose i-th column is a vector of gammas estimated given the null hypothesis m_0 = i}
#' \item{emstat}{A maxm vector of values of modified EM statistics of the model at m_0 = 1, 2, ..., maxm}
#' \item{pvals}{A maxm by 3 matrix whose i-th row indicates a vector of p-values at k = 1, 2, 3}
#' \item{aic}{A maxm vector of Akaike Information Criterion of the fitted model at m_0 = 1, 2, ..., maxm}
#' \item{bic}{A maxm vector of Bayesian Information Criterion of the fitted model at m_0 = 1, 2, ..., maxm}
#' \item{loglik}{A maxm vector of log-likelihood values of the model at m_0 = 1, 2, ..., maxm}
#' \item{penloglik}{A maxm vector of penalized log-likelihood values of the model at m_0 = 1, 2, ..., maxm}
#' \item{pmle.result}{A list of output from normalpanelmixPMLE under the number of components selected by sequantial hypothesis testing}
#' @examples
#' data(faithful)
#' attach(faithful)
#' normalpanelmixMEMtestSeq(y = eruptions)
normalpanelmixMEMtestSeq <- function (y, x = NULL, z = NULL,  maxm = 3, ninits = 10, maxit = 2000,
                                 nbtsp = 199, parallel = FALSE, cl = NULL,
                                 crit.bootstrap.from = 3) {
  # Compute the modified EM test statistic for testing H_0 of m components
  # against H_1 of m+1 components for a univariate finite mixture of normals


  if (!is.null(x))
    return (regmixMEMtestSeq(y = y, x = x, z = z, maxm = maxm, ninits = ninits, maxit = maxit,
                             nbtsp = nbtsp, parallel = parallel, cl = cl,
                             crit.bootstrap.from = crit.bootstrap.from))

  t <- dim(y)[1] #Number of year
  n <- dim(y)[2] #Number of firms
  # y   <- as.vector(y)
  p   <- 0
  if (!is.null(z)) {
    z <- as.matrix(z)
    p <- ncol(z)
    gam <- matrix(0, nrow = p, ncol = maxm)
  }
  else
    gam <- NULL

  out   <- vector('list', length = maxm)
  aic    <- bic <- double(maxm)
  pvals   <- emstat <- matrix(0, nrow = maxm, ncol = 3)
  loglik  <- penloglik <- double(maxm)

  alpha   <- mu <- sigma <- matrix(0, nrow = maxm, ncol = maxm)

  # Test H_0:m=1, H_0:m=2, ...
  binit <- NULL
  for (m in 1:maxm){
    pmle.result   <- normalpanelmixPMLE(y = y, m = m, z = z, vcov.method = "none",
                                   ninits = ninits, maxit = maxit, binit = binit)
    loglik[m] <- loglik0 <- pmle.result$loglik
    penloglik[m] <- penloglik0 <- pmle.result$penloglik
    aic[m]  <- pmle.result$aic
    bic[m]  <- pmle.result$bic

    parlist <- pmle.result$parlist
    alpha0  <- parlist$alpha
    mu0     <- parlist$mu
    sigma0  <- parlist$sigma
    gam0  <- parlist$gam

    alpha[,m] <- c(alpha0, double(maxm - m))
    mu[,m]     <- c(mu0, double(maxm - m))
    sigma[,m] <- c(sigma0, double(maxm - m))

    cat(sprintf("%d-component model estimate:\n",m))
    tab = as.matrix(rbind(alpha0, mu0, sigma0))
    rownames(tab) <- c("alpha", "mu", "sigma")
    colnames(tab) <- c(paste("comp", ".", 1:m, sep = ""))
    print(tab, digits = 4)

    if (!is.null(z)){
      gam[, m] <- gam0
      cat("gam =", gam0,"\n")
    }
    cat(sprintf("\nAIC, BIC, and log-likelihood of 1 to %.i", m), "component models \n")
    cat(c("AIC    =", sprintf(' %.2f', aic[1:m])), "\n")
    cat(c("BIC    =", sprintf(' %.2f', bic[1:m])), "\n")
    cat(c("loglik =", sprintf('%.2f', loglik[1:m])), "\n\n")

    if (m <= maxm){

      cat(sprintf("Testing the null hypothesis of %d components\n", m))

      an    <- anFormula(parlist = parlist, m = m, n = n)
      par1  <- normalpanelmixMaxPhi(y = y, parlist = parlist, z = z, an = an,
                               ninits = ninits, maxit = maxit, parallel = parallel)
      emstat.m  <- 2*(par1$penloglik - loglik0)

      # use the estimate of b as one of the initial values
      binit <- par1$coefficient

      cat(c("modified EM-test statitic ", sprintf('%.3f',emstat.m)),"\n")
      if (m <= crit.bootstrap.from) {
        em.out <- normalpanelmixCrit(y=y, parlist=parlist, z=z, values = emstat.m)
        cat(c("asymptotic p-value       ", sprintf('%.3f',em.out$pvals)),"\n \n")
      } else {
        em.out <- normalpanelmixCritBoot(y=y, parlist=parlist, z=z, values = emstat.m,
                                    ninits = ninits, nbtsp = nbtsp, parallel = parallel, cl = cl)
        cat(c("bootstrap p-value        ", sprintf('%.3f',em.out$pvals)),"\n \n")
      }
      pvals[m,]     <- em.out$pvals
      emstat[m,]    <- emstat.m
    }
  }

  for (m in 1:maxm)
    if ( pvals[m,2] >= 0.05 ) {
      cat(sprintf("\nThe number of components selected by Sequential Hypothesis Testing (alpha=0.05) = %.i", m), " \n")
      cat(sprintf("The number of components selected by AIC = %.i", which.min(aic)), " \n")
      cat(sprintf("The number of components selected by BIC = %.i", which.min(bic)), " \n")
      binit <- as.vector(c(alpha[1:m,m], mu[1:m,m], sigma[1:m,m],  gam[,m]))
      pmle.result   <- normalpanelmixPMLE(y = y, m = m, z = z,
                                     ninits = 2, maxit = maxit, binit = binit)
      cat(sprintf("\nThe summary of the estimated %.i", m), "component model: \n")
      print(summary(pmle.result))
      break
    }

  a = list(alpha = alpha, mu = mu, sigma = sigma, gam = gam, emstat = emstat, pvals = pvals, aic = aic, bic = bic, loglik = loglik, penloglik = penloglik, pmle.result = pmle.result)

  a
}  # end normalpanelmixMEMtestSeq
