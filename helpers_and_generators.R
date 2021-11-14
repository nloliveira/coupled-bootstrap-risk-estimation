##############################
######## data generators
##############################

library(MASS)
library(metR)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(foreach)
library(doSNOW)
library(glmnet)
library(doParallel)
library(rpart)
#library(bestsubset)
library(gridExtra)
library(igraph)
library(tvR)
library(flsa)


## got this function online, looking for the source to give proper credit!
withSeed <- local({
  function(expr, seed, ..., substitute=TRUE, envir=parent.frame()) {
    # Argument 'expr':
    if (substitute) expr <- substitute(expr)
    
    # Argument 'envir':
    if (!is.environment(envir))
      throw("Argument 'envir' is not a list: ", class(envir)[1L])
    
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Record entry seed
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    env <- globalenv()
    oseed <- env$.Random.seed
    # Restore on exit
    on.exit({
      if (is.null(oseed)) {
        rm(list=".Random.seed", envir=env, inherits=FALSE)
      } else {
        assign(".Random.seed", value=oseed, envir=env, inherits=FALSE)
      }
    })
    
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Set temporary seed
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    set.seed(seed=seed, ...)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Evaluate expression
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    eval(expr, envir=envir)
  } # withSeed()
})

#' Generates multivariate gaussian
#' n observations of a p-variate gaussian distribution with mean 0, variance s2, and correlation rho
#' 
#' @param n number of observations
#' @param p dimension of multivariate gaussian
#' @param rho correlation between variables
#' @param s2 variance
#' @return n x p matrix
gen_x <- function(n, p, rho, s2 = 1){
  Sig <- matrix(rho, p, p)
  diag(Sig) <- rep(1, p)
  Sig <- diag(sqrt(s2), p, p)%*%Sig%*%diag(sqrt(s2), p, p)
  X <- rmvnorm(n, rep(0, p), sigma = Sig)
  return(X)
}

#' Generates response vector - independent and normally distributed
#' 
#' @param X matrix of features
#' @param sig2 variance of noise
#' @param g_true true mean function. Maps x to E[Y]
#' @return vector in R^n 
gen_y_gaussian <- function(X, sig2, g_true){
  n <- nrow(X)
  meany <- g_true(X)
  y <- meany + rnorm(n, 0, sqrt(sig2))
  return(y)
}

##############################
######## simulation helpers, real error, plotting functions
##############################

#' Approximates via MC real risk and covariance term
#' 
#' @param X matrix of features
#' @param alphavec vector of values of alpha
#' @param sig2 variance of noise
#' @param g_true true mean function. Maps x to E[Y]
#' @param g model for which we are estimating the risk_alpha
#' @param nrep number of MC samples to approximate risk_alpha
#' @param gen_y function used to generate samples of Y
#' @return dataframe wit approximated Cov(Y, g(Y^+)), Cov(Y^+, g(Y^+)), Risk_alpha, Err_alpha for all values of alpha provided in alphavec
real_values_X <- function(X, alphavec, sig2, g_true, g, nrep, gen_y = gen_y_gaussian){
  n <- nrow(X)
  # generates many samples of Y
  y <- replicate(nrep, gen_y(X, sig2, g_true))
  if(n == 1) y <- matrix(y, 1, nrep)
  # generates perturbations omega
  w <- matrix(rnorm(n*nrep, 0, sqrt(sig2)), nrow = n, ncol = nrep)
  muvec <- g_true(X)
  
  out <- data.frame(alpha = alphavec, 
                    covygyplus = 0, 
                    covyplusgyplus = 0, 
                    erralpha = 0, 
                    riskalpha = 0)
  for(i in 1:length(alphavec)){
    ## create vectors of noisier data and compute predicted values 
    yplus <- y + sqrt(alphavec[i])*w
    gyplus <- apply(yplus, 2, g, X)
    if(n == 1) gyplus <- matrix(gyplus, 1, nrep)
    ## approximate parameters of interest
    out$riskalpha[i] <- mean(apply(gyplus, 2, function(gyplusi){mean((gyplusi - muvec)^2)}))
    out$erralpha[i] <- out$riskalpha[i] + (1+alphavec[i])*sig2
    for(j in 1:n){
      out$covygyplus[i] <- out$covygyplus[i] + cov(y[j,],gyplus[j,])/n
      out$covyplusgyplus[i] <- out$covyplusgyplus[i] + cov(yplus[j,],gyplus[j,])/n
    }
  }
  return(out)
}


