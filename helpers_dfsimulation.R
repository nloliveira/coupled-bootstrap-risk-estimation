######################################
#######  helper functions for DF simulation
######################################

#' fits Lasso regression model and predicts Y for a vector of lambdas
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param lambda sequence of lambdas
#' @return matrix of \hat Y
g_lasso_alllambdas <- function(y, xtrain, lambda){
  mod <- glmnet(xtrain, y, lambda = lambda)
  return(list(yhat = predict(mod, newx = xtrain), size_support = mod$df))
}

#' fits FS regression and predicts Y for all possible number of steps
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param steps maximum number of steps
#' @return matrix of \hat Y
g_subset_allsteps <- function(y, xtrain, steps = 100){
  mod <- fs(xtrain, y, maxsteps = steps)
  return(list(yhat = predict(mod, newx = xtrain)[,-1], size_support = mod$df))
}

#' approximates real support size and degrees of freedom for different lambdas and alphas in Lasso regression
#' 
#' @param y response vector
#' @param alphavec vector of alphas
#' @param lambda sequence of lambdas
#' @param sig2 variance
#' @param g_true true mean function
#' @param g mean function being evaluated
#' @param nrep number of MC samples
#' @param gen_y function that generates Y
#' @return matrix of \hat Y
real_values_X_alllambdas <- function(X, alphavec, lambda, sig2, g_true, g, nrep, gen_y = gen_y_gaussian){
  n <- nrow(X)
  y <- replicate(nrep, gen_y(X, sig2, g_true))
  if(n == 1) y <- matrix(y, 1, nrep)
  w <- matrix(rnorm(n*nrep, 0, sqrt(sig2)), nrow = n, ncol = nrep)
  muvec <- g_true(X)
  
  out <- matrix(0, length(alphavec), length(lambda))
  support <- matrix(0, length(alphavec), length(lambda))
  for(i in 1:length(alphavec)){
    yplus <- y + sqrt(alphavec[i])*w
    gyplus <- plyr::aaply(yplus, 2, function(yplusi){g(yplusi, X, lambda)$yhat})
    support[i,] <- plyr::aaply(yplus, 2, function(yplusi){g(yplusi, X, lambda)$size_support}) %>% apply(2, mean)
    for(j in 1:n){
      out[i,] <- out[i,] + apply(gyplus[,j,], 2, function(gyplusi_onelamb){cov(yplus[j,],gyplusi_onelamb)/(sig2*(1+alphavec[i]))})
    }
  }
  
  out <- data.frame(out)
  out$alpha <- alphavec
  names(out) <- c(as.character(lambda), "alpha")
  out <- reshape2::melt(out, id.vars = "alpha")
  names(out) <- c("alpha", "lambda", "df")
  
  support <- data.frame(support)
  support$alpha <- alphavec
  names(support) <- c(as.character(lambda), "alpha")
  support <- reshape2::melt(support, id.vars = "alpha")
  names(support) <- c("alpha", "lambda", "support")
  return(out %>% inner_join(support, by = c("alpha", "lambda")))
}

#' approximates real support size and degrees of freedom for different step sizes and alphas in FS regression
#' 
#' @param y response vector
#' @param alphavec vector of alphas
#' @param steps maximum number of steps
#' @param sig2 variance
#' @param g_true true mean function
#' @param g mean function being evaluated
#' @param nrep number of MC samples
#' @param gen_y function that generates Y
#' @return matrix of \hat Y
real_values_X_allsteps <- function(X, alphavec, steps = 100, sig2, g_true, g, nrep, gen_y = gen_y_gaussian){
  n <- nrow(X)
  y <- replicate(nrep, gen_y(X, sig2, g_true))
  if(n == 1) y <- matrix(y, 1, nrep)
  w <- matrix(rnorm(n*nrep, 0, sqrt(sig2)), nrow = n, ncol = nrep)
  muvec <- g_true(X)
  
  out <- matrix(0, length(alphavec), steps+1)
  support <- matrix(0, length(alphavec), steps+1)
  for(i in 1:length(alphavec)){
    yplus <- y + sqrt(alphavec[i])*w
    gyplus <- plyr::aaply(yplus, 2, function(yplusi){g(yplusi, X, steps)$yhat})
    support[i,] <- plyr::aaply(yplus, 2, function(yplusi){g(yplusi, X, steps)$size_support}) %>% apply(2, mean)
    for(j in 1:n){
      out[i,] <- out[i,] + apply(gyplus[,j,], 2, function(gyplusi_onelamb){cov(yplus[j,],gyplusi_onelamb)/(sig2*(1+alphavec[i]))})
    }
  }
  
  out <- data.frame(out)
  out$alpha <- alphavec
  names(out) <- c(as.character(1:(steps+1)), "alpha")
  out <- reshape2::melt(out, id.vars = "alpha")
  names(out) <- c("alpha", "steps", "df")
  
  support <- data.frame(support)
  support$alpha <- alphavec
  names(support) <- c(as.character(1:(steps+1)), "alpha")
  support <- reshape2::melt(support, id.vars = "alpha")
  names(support) <- c("alpha", "steps", "support")
  return(out %>% inner_join(support, by = c("alpha", "steps")))
}