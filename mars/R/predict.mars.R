#' Predict function for MARS
#'
#' @description Produces predictions based on the passed mars object and an new data set.
#' @usage Input a mars object and a new data set into the predict.mars function.
#'
#' @param object A mars object
#' @param newdata A new data set to create predictions
#' @param ... Additional arguments - not used in current function
#'
#' @family methods
#'
#' @return A matrix of predicted Basis Functions.
#' @export
#'
#' @examples marstest <- mars(y~x1+x2, data=marstestdata)
#' predict.mars(marstest, newdata=data.frame(x1=rnorm(250), x2=rnorm(250)))
#'

predict.mars <- function(object,newdata,...) {
  if(missing(newdata) || is.null(newdata)) {
    B <- as.matrix(object$B)
  }
  else {
    tt <- terms(object$formula,data=newdata)
    tt <- delete.response(tt)
    mf <- model.frame(tt,newdata)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)[,-1] # remove intercept
    B <- make_B(X,object$Bfuncs)
  }
  beta <- object$coefficients
  drop(B %*% beta)
}


### make_B() function
make_B <- function(X,object) {
  B <- init_B(nrow(X),length(object)-1)
  Bfuncs <- object
  for(i in 2:length(Bfuncs)) {
    for(j in 1:nrow(Bfuncs[[i]])) {
      v <- Bfuncs[[i]][j,]["v"]
      s <- Bfuncs[[i]][j,]["s"]
      t <- Bfuncs[[i]][j,]["t"]
      if (j == 1){
        B[,i] <- h(X[,v],s,t)
      }
      else if(j != 1){
        B[,i] <- B[,i] * h(X[,v],s,t)
      }
    }
  }
  B <- as.matrix(B)
  return(B)
}
