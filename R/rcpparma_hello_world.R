#'@export
#'@importFrom Rcpp evalCpp 
#'@useDynLib normalRegPanelMix
#'
rcpparma_hello_world <- function() {
  .Call('normalRegPanelMix_rcpparma_hello_world', PACKAGE = 'normalRegPanelMix')
}