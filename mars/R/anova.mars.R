#' anova function for MARS
#'
#' @description Provides a ANOVA of information regarding the passed mars object.
#' @usage Input the mars object from the mars function into the anova function.
#'
#' @param object A mars object.
#' @param ... Additional arguments - not used in current function.
#'
#' @family methods
#'
#' @return ANOVA of the passed mars object
#' @export
#'
#' @examples anova(mars(formula=y~x1+x2, data = marstestdata))
#'
#'

anova.mars <- function(object,...) {
  n1 <- nrow(object$model)
  n2 <- nrow(other[[1]]$model)
  df1 <- ncol(object$Bfuncs) + 1
  df2 <- ncol(other[[1]]$Bfuncs) + 1
  RSS1 <- sum(object$residuals^2)
  RSS2 <- sum(other[[1]]$residuals^2)
  F <- ((RSS1 - RSS2) / (df1 - df2)) / (RSS2 / (n2 - df2))
  pval <- pf(F, df1 - df2, n2 - df2, lower.tail = FALSE)
  return(list(statistic = F, parameter = df1 - df2, p.value = pval))
}
