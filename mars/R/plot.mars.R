#' Plot function for MARS
#'
#' @description Provides a plot of information regarding the passed mars object.
#' @usage Input the mars object from the mars function into the plot function.
#'
#' @param x A mars object.
#' @param ... Additional arguments - not used in current function.
#'
#' @family methods
#'
#' @return A plot of the passed mars object
#' @export
#'
#' @examples plot(mars(formula=y~x1+x2, data = marstestdata))
#'
#'
plot.mars <- function(x,...) {
  data <- eval(x$call$data)
  tt <- terms(x$formula, data=data)
  tt <- delete.response(tt)
  mf <- model.frame(tt,data)
  mt <- attr(mf,"terms")
  X <- model.matrix(mt,mf)[-1] #remove intercept
  Bf <- x$Bfuncs
}

