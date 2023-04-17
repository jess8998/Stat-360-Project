#' Multivariate Adaptive Regression Splines (MARS)
#'
#' @description The formula for regression modeling analysis utilizing Friedman's MARS approach to determine the optimal curve of fit for nonlinear functions.
#' @usage mars(formula, data,control)

#' @param formula An R formula
#' @param data A data frame containing the data
#' @param control An object of class `mars.control`
#' @details This function uses Friedman's MARS technique to fit a nonlinear regression model to the data. It uses a forward-backward greedy search algorithm to select basis functions, and cuts away the basis functions using cross-validation. The resulting model is returned as an object of class `mars`, which can be further analyzed using various functions such as `plot`, `anova`, `predict`, `summary`, and `print`.
#'
#' @return A ‘mars’ object containing the final regression result as well as details about the basis functions used.
#' The anova, plot, predict, print and summary methods may be used for better understanding of the mars output.
#' @export marstestdata
#'
#' @seealso [mars.control] for creating a control object to be passed as an argument to the mars function.
#' @seealso [anova.mars] for comparing the performance of mars models, and selecting the best mars model from the set of models.
#' @seealso [plot.mars] for creating plots that show the results of the mars function.
#' @seealso [predict.mars] to generate predictions based on results obtained from fitting a mars model using the mars function.
#' @seealso [print.mars]  to display the ‘Call’ and ‘Coefficients’ arguments of the mars function.
#' @seealso [summary.mars] to generate a summary of the results acquired from the mars function.
#'
#' @examples marstest <- mars(y~., data=marstestdata)
#'
#' @import stats
#' @import ISLR
#' @author Adriel Stanley Vidianto 301407418
#' Sumin Kang 301376418
#' Jesslyn Devina Anwar 301468731
#'
#' @references Jerome H. Friedman. "Multivariate Adaptive Regression Splines". The Annals of Statistics, 19(1). 1-67, March,1991.
#' \url{https://doi.org/10.1214/aos/1176347963}


mars <- function(formula,data,control=NULL) {
  cc <- match.call() # save the call
  mf <- model.frame(formula,data)
  y <- model.response(mf)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf) [,-1,drop=FALSE]
  x_names <- colnames(x)
  control <- validate_mars.control(control)
  fwd <- fwd_stepwise(y,x,control)
  bwd <- bwd_stepwise(fwd,control)
  #fit the final best model and output results
  fit <- lm(y~.-1,data=data.frame(y=y,bwd$B))
  out <- c(list(call=cc,formula=formula,y=y,B=bwd$B,Bfuncs=bwd$Bfuncs,
                x_names=x_names),fit)
  class(out) <- c("mars",class(fit))
  out
}

fwd_stepwise <-function(y,x,control=mars.control()){
  #initialize N,n,B, and Bfuncs
  N <- length(y) # sample size
  n <- ncol(x) # number of predictors
  B <- init_B(N,control$Mmax)
  Bfuncs <- vector(mode="list", length = control$Mmax+1)

  for(i in 1:(control$Mmax/2)) { # contrast to indexing 2...Mmax in Friedman
    M <- (2*i) -1
    lof_best <- Inf
    for(m in 1:M) { # choose a basis function to split
      svars <- setdiff(1:n, Bfuncs[[m]][,"v"])
      for(v in svars) { #select a variable to split on
        tt <- split_points(x[,v],B[,m])
        for(t in tt) {
          Bnew <- data.frame(B[,1:M],
                             Btem1=B[,m]*h(x[,v],+1,t),
                             Btem2=B[,m]*h(x[,v],-1,t))
          gdat <- data.frame(y=y,Bnew)
          lof <- LOF(y~.,gdat, control)
          if(lof < lof_best) {
            lof_best <- lof
            split_best <- c(m=m,v=v,t=t)
          } # end if
        } #end for loop over t (the splits)
      } #end for loop over v (the variables)
    } #end for loop over m (the basis function)
    mstar <- split_best["m"]
    vstar <- split_best["v"]
    tstar <- split_best["t"]
    B[,M+1] <- B[,mstar]*h(x[,vstar],-1,tstar)
    B[,M+2] <- B[,mstar]*h(x[,vstar],+1,tstar)
    Bfuncs[[M+1]] <- rbind(Bfuncs[[mstar]],c(s=-1,vstar,tstar))
    Bfuncs[[M+2]] <- rbind(Bfuncs[[mstar]],c(s=1,vstar,tstar))
  } # end loop over M
  colnames(B) <- paste0("B", (0:(ncol(B)-1)))
  names(y) <- c(1:length(y))
  return(list(y=y,B=B,Bfuncs=Bfuncs))
}

init_B <- function(N,Mmax) {
  B <- data.frame(matrix(NA,nrow=N,ncol=(Mmax+1)))
  B[,1] <- 1
  names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)
}

bwd_stepwise <- function(fwd,control) {
  #fwd is a list with elements y, B and B funcs
  Mmax <- (ncol(fwd$B)-1)
  Jstar <- 2:(Mmax+1)
  Kstar <- Jstar
  dat <- data.frame(y=fwd$y,fwd$B)
  lofstar <- LOF(y~.-1,dat,control)
  for(M in (Mmax+1):2)  {
    b <- Inf
    L <- Kstar
    if(control$trace) cat("L:", L, "\n")
    for(m in L) {
      K <- setdiff(L,m)
      dat <- data.frame(y=fwd$y,fwd$B[,K])
      lof <- LOF(y~.,dat,control)
      if(lof < b) {
        b <- lof
        Kstar <- K
      }
      if(lof < lofstar) {
        lofstar <- lof
        Jstar <- K
      }
    }
  }
  Jstar <- c(1,Jstar)
  return(list(y=fwd$y,B=fwd$B[,Jstar],Bfuncs=fwd$Bfuncs[Jstar]))
  #output from bwd: y, B, Bfuncs
}

LOF <- function(formula, data, control) {
  mod <- lm(formula,data)
  RSS <- sum((mod$res)^2)
  d <- control$d
  N <- nrow(data)
  M <- length(coefficients(mod))-1
  cM <- sum(hatvalues(mod))
  cTilde <- sum(diag(hatvalues(mod)))+d*M
  gcv <- RSS*N/((N-cTilde)^2)
  return(gcv)
}

h <- function(x, s, t) {
  ifelse(x>t, pmax(0, x - t), pmax(0, t - x))
}

split_points <- function(xv,Bm) {
  out <- sort(unique(xv[Bm>0]))
  (out <- out[-length(out)])
}



#------------------------------------------------------------------------
# constructor, validator and helper for class mars.control
#------------------------------------------------------------------------
#

new_mars.control <- function(control) {
  structure(control,class="mars.control")
}

validate_mars.control <- function(control) {
  stopifnot(is.integer(control$Mmax),
            is.numeric(control$d),
            is.logical(control$trace))
  if(control$Mmax < 2) {
    warning("Mmax must be >= 2; Reset it to 2")
    control$Mmax <- 2}
  if(control$Mmax %% 2 > 0) {
    control$Mmax <- 2*ceiling(control$Mmax/2)
    warning("Mmax should be an even integer. Reset it to ",control$Mmax)}
  control
}

#' 'Contructor for `mars.control` objects
#'
#' This function constructs a `mars.control` object that is specifies
#' the parameters used in the model fitting procedure.
#'
#' @param Mmax Maximum number of basis functions. Should be an even integer. Default value is 2.
#' @param d The coefficient in the penalty term of the generalized cross validation measure. Default is 3.
#' @param trace Whether we should print status information about the fitting. Default is `FALSE`.
#'
#' @return A `mars.control` object
#' @export marstestdata
#'
#' @examples mc <- mars.control(Mmax=10)

mars.control <- function(Mmax=2,d=3,trace=FALSE) {
  Mmax <- as.integer(Mmax)
  control <- list(Mmax=Mmax,d=d,trace=trace)
  control <- validate_mars.control(control)
  new_mars.control(control)
}
