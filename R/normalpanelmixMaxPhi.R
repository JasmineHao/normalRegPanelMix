#' @description Compute ordinary & penalized log-likelihood ratio resulting from
#' MEM algorithm at k=1,2,3.
#' @export
#' @title normalpanelmixMaxPhi
#' @name normalpanelmixMaxPhi
#' @param y n by 1 vector of data
#' @param parlist The parameter estimates as a list containing alpha, mu, sigma, and gamma
#' in the form of (alpha = (alpha_1, ..., alpha_m), mu = (mu_1, ..., mu_m),
#' sigma = (sigma_1, ..., sigma_m), gam = (gamma_1, ..., gamma_m))
#' @param z n by p matrix of regressor associated with gamma
#' @param an a term used for penalty function
#' @param tauset A set of initial tau value candidates
#' @param ninits The number of randomly drawn initial values.
#' @param epsilon.short The convergence criterion in short EM. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon.short}.
#' @param epsilon The convergence criterion. Convergence is declared when the penalized log-likelihood increases by less than \code{epsilon}.
#' @param maxit.short The maximum number of iterations in short EM.
#' @param maxit The maximum number of iterations.
#' @param verb Determines whether to print a message if an error occurs.
#' @param parallel Determines whether package \code{doParallel} is used for calculation
#' @param cl Cluster used for parallelization; if it is \code{NULL}, the system will automatically
#' create a new one for computation accordingly.
#' @return A list with items:
#' \item{loglik}{Log-likelihood resulting from MEM algorithm at k=1,2,3.}
#' \item{penloglik}{Penalized log-likelihood resulting from MEM algorithm at k=1,2,3.}
normalpanelmixMaxPhi <- function (y, parlist, z = NULL, an, tauset = c(0.1,0.3,0.5),
                             ninits = 10, epsilon.short = 1e-02, epsilon = 1e-08,
                             maxit.short = 500, maxit = 2000,
                             verb = FALSE,
                             parallel = FALSE,
                             cl = NULL,update.alpha=1) {
  # Given a parameter estimate of an m component model and tuning paramter an,
  # maximize the objective function for computing the modified EM test statistic
  # for testing H_0 of m components against H_1 of m+1 for a univariate normal finite mixture

  warn  <- options(warn=-1) # Turn off warnings

  if (!is.null(z)) {
    z     <- as.matrix(z)
    p     <- ncol(z)
  } else {
    p <- 0
  }
  m <- length(parlist$alpha)

  ninits.short <- ninits*10*(1+p)*m

  loglik.all <- matrix(0,nrow=m*length(tauset),ncol=3)
  penloglik.all <- matrix(0,nrow=m*length(tauset),ncol=3)
  coefficient.all <- matrix(0,nrow=m*length(tauset),ncol=(3*(m+1)+p))

  if (parallel) {
    if (is.null(cl)){
      cl <- makeCluster(detectCores())
      # print(detectCores())
    }
    registerDoParallel(cl)
    results <- foreach (t = 1:length(tauset),
                        .export = 'normalpanelmixMaxPhiStep', .combine = c)  %:%
      foreach (h = 1:m) %dopar% {
        normalpanelmixMaxPhiStep (c(h, tauset[t]), y, parlist, z, p,
                             an,
                             ninits, ninits.short,
                             epsilon.short, epsilon,
                             maxit.short, maxit,
                             verb,update.alpha) }
    on.exit(cl)
    loglik.all <- t(sapply(results, "[[", "loglik"))
    penloglik.all <- t(sapply(results, "[[", "penloglik"))
    coefficient.all <- t(sapply(results, "[[", "coefficient"))
  }
  else
    for (h in 1:m)
      for (t in 1:length(tauset)) {
        rowindex <- (t-1)*m + h
        tau <- tauset[t]
        result <- normalpanelmixMaxPhiStep(c(h, tau), y, parlist, z, p,
                                      an,
                                      ninits, ninits.short,
                                      epsilon.short, epsilon,
                                      maxit.short, maxit,
                                      verb,update.alpha)
        loglik.all[rowindex,] <- result$loglik
        penloglik.all[rowindex,] <- result$penloglik
        coefficient.all[rowindex,] <- result$coefficient
      }

  loglik <- apply(loglik.all, 2, max)  # 3 by 1 vector
  penloglik <- apply(penloglik.all, 2, max)  # 3 by 1 vector
  index <- which.max(loglik.all[ ,3]) # a par (h,m) that gives the highest likelihood at k=3
  coefficient <- as.vector(coefficient.all[index,])

  out <- list(coefficient = coefficient, loglik = loglik, penloglik = penloglik)

  out


}  # end normalpanelmixMaxPhi
