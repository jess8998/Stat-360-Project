#' Summary function for MARS
#'
#' @description Provides a summary of information regarding the passed mars object.
#' @usage Input the mars object from the mars function into the summary function.
#'
#' @param object A mars object.
#' @param ... Additional arguments - not used in current function.
#'
#' @family methods
#'
#' @return A summary of the passed mars object
#' Inlcudes the Call, Coefficients, Residuals and Basis Function Components.
#' @export
#'
#' @examples summary(mars(formula=y~x1+x2, data = marstestdata))
#'

summary.mars <- function(object,...) {
  cat("Call:\n")
  print(object$call)
  cat("\nCoefficients:\n")
  print(object$coefficients)
  cat("\nResiduals:\n")
  Min <- min(object$residuals, na.rm=TRUE)
  Q1 <- quantile(object$residuals,0.25, na.rm=TRUE)
  Median <- median(object$residuals, na.rm=TRUE)
  Q3 <- quantile(object$residuals,0.75,na.rm=TRUE)
  Max <- max(object$residuals, na.rm=TRUE)
  resid <- cbind(Min, Q1, Median, Q3, Max)
  row.names(resid) <- "Residuals"
  print(resid)
  cat("\nBasis Function Components:\n")
  for(i in 1:length(object$Bfuncs)) {
    if(i == 1) {
      cat(paste0("B",i-1, ":"))
      cat("\nIntercept\n")
    } else {
      cat(paste0("\nB",i-1, ":"))
      for(j in 1:nrow(object$Bfuncs[[i]])) {
        cat(paste0("\nComponent ",j,": variable x",object$Bfuncs[[i]][j,]["v"],
                   "; sign ", object$Bfuncs[[i]][j,]["s"],
                   "; knot at ",object$Bfuncs[[i]][j,]["t"], "\n"))
      }
    }
  }
}
