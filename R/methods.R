normalregMix.test.on <- FALSE
normalregMix.test.seed <- 8888577

#' @description Turns on/off the test mode.
#' @export
#' @title testMode
#' @name testMode
#' @description When the modified EM-algorithm is run, initial values are randomly created
#' based on the data given. If the test mode is turned on, these initial values
#' are going to be created with the random seed provided. This method would be suitable
#' for users who would like to replicate experiments. By default, the test mode is turned off.
#' @param on Option to turn on the test mode
#' @param seed The random seed to be used for initialization
#' @param hide.message Determines whether to print the current seed and status
testMode <- function(on = FALSE, seed = 8888577, hide.message = TRUE)
{
  unlockBinding("normalregMix.test.on", getNamespace("normalregMix"))
  unlockBinding("normalregMix.test.seed", getNamespace("normalregMix"))
  assign("normalregMix.test.on", on, getNamespace("normalregMix"))
  assign("normalregMix.test.seed", seed, getNamespace("normalregMix"))
  lockBinding("normalregMix.test.on", getNamespace("normalregMix"))
  lockBinding("normalregMix.test.seed", getNamespace("normalregMix"))

  if (!hide.message)
    print(paste("The test mode is currently",
                switch(as.character(normalregMix.test.on), "TRUE" = "ON", "FALSE" = "OFF"),
                "with seed",
                as.character(normalregMix.test.seed)))
}

#' @description Prints the multiple plots for diagnosis. 
#' @export
#' @title plotDiag
#' @name plotDiag
#' @param components n vector of the indices of components for each observation
#' @param y n vector of data that represents dependent variables
#' @param x n by q matrix of data that represent(s) covariates
#' @param m The number of components
#' @examples
#' data(faithful)
#' attach(faithful)
#' regmix.model.m2 <- regmixPMLE(y = eruptions, x = waiting, m = 2)
#' \dontrun{plotDiag(regmix.model.m2$components, y = eruptions, x = waiting, m  = 2)}
plotDiag <- function(components, y = y, x = x, m = 2)
{
  dimx <- dim(as.matrix(x))[2]
  if ((dimx <= 1) || (is.null(dimx)))
    return (NULL)

  pivot.names <- as.character(seq(1, dimx))
  if (!is.null(colnames(x)))
    pivot.names <- colnames(x)

  ivs <- as.matrix(x)

  for (j in 1:m)
  {
    ivs.component <- as.matrix(ivs[components == j,])
    ys.component <- y[components == j]
    for (pivot in 1:dimx)
    {
      pivot.name <- pivot.names[pivot]
      ivs.pivot <- ivs.component[,pivot]
      ivs.others <- ivs.component[,-pivot]
      lm.y.other <- lm(ys.component ~ ivs.others)
      lm.pivot.other <- lm(ivs.pivot ~ ivs.others)
      plot.df <- data.frame(y.on.others = lm.y.other$residuals,
                            pivot.on.others = lm.pivot.other$residuals)
      plot <- ggplot(plot.df, aes(x=pivot.on.others, y=y.on.others))
      plot <- plot + geom_point(shape=1) + geom_smooth(method=lm) +
        xlab(paste("pivot.on.others (pivot on ", pivot.name, ", component ",
                   as.character(j), ")", sep = ""))
      print(plot)
    }
  }
}

#' @description Generates a vector that indicates which component each observation belongs to,
#' based on its posterior probability
#' @export
#' @title getComponentcomponents
#' @name getComponentcomponents
#' @param postprobs n by m matrix of posterior probabilities for
#' m-component model on n observations
#' @return n by 1 vector of components that indicate which component each observation belongs to
#' based on its posterior probability
getComponentcomponents <- function(postprobs)
{
  postprobs.mat <- as.matrix(postprobs)
  apply(postprobs.mat, 1, function(i) (which(i==max(i))))
}

#' @description Returns the summary of a normalregMix instance
#' @export
#' @title summary.normalregMix
#' @name summary.normalregMix
#' @param object normalregMix instance
#' @param reorder Determines whether components are reordered in summary
#' @param digits Digits used for reporting
#' @param ... Other arguments that do not affect the function
summary.normalregMix <- function(object, reorder = FALSE, digits = 3, ...) {

if (object$label == "PMLE") {
  coef <- object$coefficients
  vcov <- object$vcov
  m    <- object$m
  # When m=1, the first element of coef is alpha (=1), so we drop it
  if (object$m == 1) { coef <- coef[-1] }

  if (reorder && m != 1){
    len1 <- length(object$coefficients)
    k <- nrow(object$parlist$mubeta)
    if (is.null(k)) { k <- 1 }
    sel.vec <- NULL
    for (j in 1:m){
      sel <- c(j, (m+(j-1)*k+1):(m+j*k), j+m*(k+1))
      sel.vec <- c(sel.vec,sel)
    }
    if (!is.null(object$parlist$gam)) {
      p <- length(object$parlist$gam)
      sel <- c((len1-p+1):len1)
      sel.vec <- c(sel.vec,sel)
    }
    reorder.mat <- matrix(0, nrow=len1, ncol=len1)
    sel.vec.2 <- cbind(1:len1, sel.vec)
    reorder.mat[sel.vec.2] <- 1
    coef <- coef[sel.vec]
    vcov <- reorder.mat %*% vcov %*% t(reorder.mat)
  }

  se    <- sqrt(diag(vcov))
  tval  <- coef / se
  TAB   <- cbind(Estimate = coef, StdErr =se, t.value = tval, p.value = 2*pnorm(-abs(tval)))
  res   <- list(coefficients = TAB, parlist = object$parlist, vcov = vcov,
              loglik = object$loglik, penloglik = object$penloglik,
              aic = object$aic, bic = object$bic, call = object$call,
              m = object$m, reorder = reorder,
              digits = digits)
  class(res) <- "summary.normalregMix"
  res

} else {
  stop("The object is not an output of normalmixPMLE/regmixPMLE.")
}

}  # end function summary.normalregMix

#' @description Prints the summary of a normalregMix instance
#' @export
#' @title print.summary.normalregMix
#' @name print.summary.normalregMix
#' @param x normalregMix instance
#' @param ... Other arguments that do not affect the function
print.summary.normalregMix <- function(x, ...) {

cat("\nCall:\n")
print(x$call)
cat("\n")

m   <- x$m
tab <- x$coefficients
reorder <- x$reorder
k <- nrow(x$parlist$mubeta)
if (is.null(k)) { k <- 1 }

if (reorder) {
#   cat(sprintf("Number of components: %d\n",m))
  for (j in 1:m){
    cat(sprintf("Component %i\n",j))
    coef.j = tab[c(((k+2)*(j-1)+1):((k+2)*j)), ]
    # rownames(coef.j) <- c("alpha","mu","sigma")
    printCoefmat(coef.j, P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE, digits=x$digits)
  }
  if (!is.null(x$parlist$gam)) {
    p <- length(x$parlist$gam)
    gam <- tab[(nrow(tab)-p+1):nrow(tab), , drop = FALSE]
    cat("gam\n")
    printCoefmat(gam, P.values = TRUE, has.Pvalue = TRUE, signif.legend = FALSE, digits=x$digits)
  }
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
} else {
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, digits = x$digits)
}

cat(sprintf("\nloglik at estimate:  %.3f\n", x$loglik))
cat(sprintf("penloglik at estimate: %.3f\n", x$penloglik))
cat(sprintf("AIC: %.3f\n", x$aic))
cat(sprintf("BIC: %.3f\n", x$bic))

}

#' @description Prints the description of a normalregMix instance
#' @export
#' @title print.normalregMix
#' @name print.normalregMix
#' @param x normalregMix instance
#' @param ... Other arguments that do not affect the function
print.normalregMix <- function(x, ...) {

cat("\nCall:\n")
print(x$call)

if (x$label == "MEMtest") {
  cat(sprintf("\nTesting the null hypothesis of %d components\n", x$m))
  cat("                            k = 1  k = 2  k = 3 \n")
  cat(c("modified EM-test statistic ",sprintf('%.3f ', x$emstat)),"\n")
  if (x$crit.method == "asy") {
  cat(c("asymptotic p-value         ",sprintf('%.3f ', x$pvals)),"\n")
  } else if (x$crit.method == "boot") {
    cat(c("bootstrap p-value          ",sprintf('%.3f ', x$pvals)),"\n")
  }
} else if (x$label == "PMLE") {
  cat("\nCoefficients:\n")
  print(x$coefficients, digits=4)
  cat(sprintf("\nloglik at estimate: %.3f\n", x$loglik))
  cat(sprintf("penloglik at estimate: %.3f\n", x$penloglik))
  cat(sprintf("AIC: %.3f\n", x$aic))
  cat(sprintf("BIC: %.3f\n", x$bic))
} else {
  stop("The object is not a valid normalMix object.")
}
}
