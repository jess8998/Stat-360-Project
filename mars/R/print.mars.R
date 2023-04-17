#' Print function for MARS
#'
#' @description Provides brief information regarding the passed mars object.
#' @usage Call the mars function.
#'
#' @param x A mars object.
#' @param ... Additional arguments - not used in current function.
#'
#' @family methods
#'
#' @return The mars object output which include Call and Coefficients.
#' @export
#'
#' @examples mars(formula=y~x1+x2, data = marstestdata)
#'
#'
print.mars <- function(x,...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}
