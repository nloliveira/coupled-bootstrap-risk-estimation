##############################
######## different methods to compute error
##############################

#' Computes \hatErr_\alpha 
#' 
#' @param y response vector
#' @param w matrix of B noise vectors
#' @param gyplus matrix of g fitted on y+sqrt alpha w[,j]
#' @param alpha noise size
#' @param sig2 variance 
#' @param B_seq either vector or scalar. If scalar, number of bootstrap samples B. If vector, sequence of number of bootstrap samples that are used to produce estimates
#' @return data frame where for each value of B, returns estimated training error, out of sample error, risk, and covariance term
riskAlpha <- function(y, w, gyplus, alpha, sig2, B_seq){
  B <- ncol(w)
  errminusplus <- numeric(B)
  errtraining <- numeric(B)
  varterms <- numeric(B)
  out <- data.frame(B = B_seq, method = "CB", alpha = alpha, training = NA, error = NA, risk = NA, cov = NA)
  j = 1
  for(i in 1:B){
    yminus <- y - w[,i]/sqrt(alpha)
    yplus <- y + w[,i]*sqrt(alpha)
    errminusplus[i] <- mean((c(yminus) - gyplus[,i])^2)
    errtraining[i] <- mean((c(yplus) - gyplus[,i])^2)
    varterms[i] <- sig2*alpha - mean(w[,i]^2)/alpha
    if(i == B_seq[j]){
      out$training[j] <- mean(errtraining[1:i])
      out$error[j] <- mean(errminusplus[1:i]) + mean(varterms[1:i])
      out$risk[j] <- out$error[j] - sig2*(1+alpha)
      out$cov[j] <- (out$error[j] - out$training[j])/2
      j <- j+1
    }
  }
  return(out)
}



#' Computes \hatCov_\alpha 
#' 
#' @param y response vector
#' @param yhat g(y, X)
#' @param w matrix of B noise vectors
#' @param gyplus matrix of g fitted on y+sqrt alpha w[,j]
#' @param alpha noise size
#' @param sig2 variance 
#' @param B_seq either vector or scalar. If scalar, number of bootstrap samples B. If vector, sequence of number of bootstrap samples that are used to produce estimates
#' @return data frame where for each value of B, returns estimated training error \|Y-g(Y)\|^2_2, noisier training error \|Y^+-g(Y^+)\|^2_2, and BY covariance term
onlycov_uncentered <- function(y, yhat, w, gyplus, alpha, sig2, B_seq){
  B <- ncol(w)
  n <- length(y)
  ynew <- matrix(NA, n, B)
  train_noisy <- numeric(B)
  out <- data.frame(B = B_seq, method = "uncentered_cov", training_noisier = NA, training = NA, cov = NA)
  j = 1
  for(i in 1:B){
    train_noisy[i] <- mean((y+sqrt(alpha)*w[,i] - gyplus[,i])^2)
    if(i == B_seq[j]){
      covii <- numeric(n)
      for(k in 1:n){ 
        covii[k] <- mean(sqrt(alpha)*w[k,1:i]*(gyplus[k,1:i] - yhat[k]))
      }
      out$cov[j] <- mean(covii)
      out$training[j] <- mean((y - yhat)^2)
      out$training_noisier[j] <- mean(train_noisy[1:i])
      j <- j+1
    }
  }
  return(out)
}

