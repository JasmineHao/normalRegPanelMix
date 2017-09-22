#'
#'@export
#'@description return panel data likelihood
#' @title cppregpanelmixMLE
#' @name cppregpanelmixMLE
#' @param alpha m by 1 vector
#' @param mu m by 1 vector
#' @param sigma m by 1 vector
#' @param gamma p by 1 vector
#' @param beta m by q matrix
#' @param y nt by 1 vector of data
#' @param x nt by q matrix of data for x (if exists)
#' @param z nt by p matrix of regressor associated with gamma
#' @param m The number of components in the mixture
#' @param p The number of variable in z
#' @param q The number of variable in x
#' @return  log-likelihood


cppregpanelmixMLE <- function(ys, xs, zs, alpha0s, mu0s, sigma0s, beta0s, gamma0s, m, q, p, t) {
  .Call('normalRegPanelMix_cppregpanelmixMLE', PACKAGE = 'normalRegPanelMix', ys, xs, zs, alpha0s, mu0s, sigma0s, beta0s, gamma0s, m, q, p, t)
}