#' @export

cppregpanelmixMLE <- function(ys, xs, zs, alpha0s, mu0s, sigma0s, beta0s, gamma0s, m, q, p, t) {
    .Call('normalRegPanelMix_cppregpanelmixMLE', PACKAGE = 'normalRegPanelMix', ys, xs, zs, alpha0s, mu0s, sigma0s, beta0s, gamma0s, m, q, p, t)
}
#' @export
cppnormalpanelmixPMLE <- function(bs, ys, zs, mu0s, sigma0s, m, p, t, an, maxit = 2000L, ninits = 10L, tol = 1e-8, tau = 0.5, h = 0L, k = 0L) {
    .Call('normalRegPanelMix_cppnormalpanelmixPMLE', PACKAGE = 'normalRegPanelMix', bs, ys, zs, mu0s, sigma0s, m, p, t, an, maxit, ninits, tol, tau, h, k)
}
