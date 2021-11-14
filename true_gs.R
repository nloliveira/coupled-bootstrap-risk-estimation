#########
## library of functions that can be used as true mean functions. 
## Takes matrix X as input and outputs E[Y|X] = g(X)
#########

#' true mean of Y|X as a linear regression, all coefficients different from 0
#' 
#' @param X matrix of features
#' @return vector of E[Y]
g_true_lin <- function(X){
  withSeed({
    beta <- runif(ncol(X),-.5,.5)
  }, 57)
  if(is.vector(X)){
    return(beta%*%X)
  } else{
    return(as.vector(X%*%beta))
  }
}

#' true mean of Y|X is a linear function where 1/3 of the coefficients are equal to 0
#' 
#' @param X matrix of features
#' @return vector of E[Y]
g_true_lin_sparse <- function(X){
  withSeed({
    beta <- runif(ncol(X),-5,5)
  }, 57)
  beta[1:round(ncol(X)/3)] <- 0
  if(is.vector(X)){
    return(beta%*%X)
  } else{
    return(as.vector(X%*%beta))
  }
}

#' true mean of Y|X as a sparse linear regression with at most 5 features with coef diff than 0
#' 
#' @param X matrix of features
#' @return vector of E[Y]
g_true_lin_sparse_5 <- function(X){
  withSeed({
    beta <- runif(ncol(X),-.5,.5) ## withseed ensures that the generated vector beta is always the same
  }, 57)
  if(ncol(X)>5){
    beta[6:ncol(X)] <- 0
  } else{
    beta[1:round(ncol(X)/3)] <- 0
  }
  if(is.vector(X)){
    return(beta%*%X)
  } else{
    return(as.vector(X%*%beta))
  }
}

#' true mean of Y|X as a regression tree
#' 
#' @param X matrix of features. needs to have at least 7 features
#' @return vector of E[Y]
g_true_tree <- function(X){
  if(is.vector(X)){
    1.5*(X[1]>=1) +  0.5*(X[1]<1) -2*(X[2]>=1)*(X[1]<0) + 6*(X[3]>7) - 1.5*(X[1]<=-3) + 3*(X[5]>=0) + (X[5]<=1.5)*(X[7]>=1) - (X[1]<=0)*(X[2]<=0)*(X[4]<=1)
  } else{
    1.5*(X[,1]>=1) +  0.5*(X[,1]<1) -2*(X[,2]>=1)*(X[,1]<0) + 6*(X[,3]>7) - 1.5*(X[,1]<=-3) + 3*(X[,5]>=0) + (X[,5]<=1.5)*(X[,7]>=1) - (X[,1]<=0)*(X[,2]<=0)*(X[,4]<=1)
  }
}

#' true mean of Y|X as a regression tree in a 'checkboard' pattern
#' 
#' @param X matrix of features. needs to have at least 2 features
#' @return vector of E[Y]
g_true_tree_notlin <- function(X){
  if(is.vector(X)){
    2*(X[1]>=0)*(X[2]>=0) + 1*(X[1]<0)*(X[2]<0) - 1.5*(X[1]>0)*(X[2]<0) -1.5*(X[1]<0)*(X[2]>0)  
  } else{
    2*(X[,1]>=0)*(X[,2]>=0) + 1*(X[,1]<0)*(X[,2]<0) - 1.5*(X[,1]>0)*(X[,2]<0) -1.5*(X[,1]<0)*(X[,2]>0)  
  }
}

