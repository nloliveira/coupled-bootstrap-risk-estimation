######################################
#######  helper functions for image denoising application
######################################

normalize <- function(img, up = 1, low = 0){
  out = (img - min(img))/(max(img) - min(img))
  out = out*(up-low) + low
  return(out)
}

##' estimates df fused lasso used for image denoising
##' \hat df = number of 'blocks' of the same color in an image
##' idea: transform image matrix into a graph. each pixel is a node, pixels are connected if they are the same color
##' uses igraph to count the number of communities in the graph
##' 
##' @param mat matrix equivalent to image. each entry is a number representing the pixel color
##' @return number. unbiased estimate of df
count_df_image <- function(mat){
  ## creates graph
  edges <- NULL
  mat_indices <- matrix(1:(nrow(mat)*ncol(mat)), nrow(mat), byrow = T)
  for(i in 1:(nrow(mat)-1)){
    for(j in 1:(ncol(mat)-1)){
      if(mat[i,j] == mat[i,j+1]) edges <- c(edges, mat_indices[i,j], mat_indices[i,j+1])
      if(mat[i,j] == mat[i+1,j]) edges <- c(edges, mat_indices[i,j], mat_indices[i+1,j])
    }
  }
  i = nrow(mat)
  for(j in 1:(ncol(mat)-1)){
    if(mat[i,j] == mat[i,j+1]) edges <- c(edges, mat_indices[i,j], mat_indices[i,j+1])
  }
  j = ncol(mat)
  for(i in 1:(nrow(mat)-1)){
    if(mat[i,j] == mat[i+1,j]) edges <- c(edges, mat_indices[i,j], mat_indices[i+1,j])
  }
  
  img_graph <- make_graph(edges, directed = FALSE, n = nrow(mat)*ncol(mat))
  return(count_components(img_graph))
}

##' adds gaussian noise to image
##' 
##' @param img matrix equivalent to image. each entry is a number representing the pixel color
##' @param sd standard deviation of the noise to be added to the entire image
##' @return matrix. noised image
add_gauss_noise <- function(img, sd){
  img + array(rnorm(nrow(img)*ncol(img), sd=sd), dim(img))
}

##' real risk
##' 
##' @param img original image
##' @param lambda_seq sequence of lambda values to be used for fused lasso
##' @param s2 variance
##' @param nsim number of replications for MC approximation
##' @return data frame with columns lambda, risk
true_risk <- function(img, lambda_seq, s2, nsim = 100, add_noise = add_gauss_noise){
  
  normalize <- function(img, up = 1, low = 0){
    out = (img - min(img))/(max(img) - min(img))
    out = out*(up-low) + low
    return(out)
  }
  
  iterations = nsim
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- makeCluster(detectCores())
  registerDoSNOW(cl)
  results <- foreach(i = 1:iterations, .packages=c("flsa"), .options.snow = opts) %dopar% {
    noised <- normalize(add_noise(img, sd = sqrt(s2)))
    sol <- flsa(noised, lambda1 = 0, lambda2 = lambda_seq)
    apply(sol, 1, function(el){
      mean((el - img)^2)
    })
  }
  stopCluster(cl)
  rm(cl)
  out <- plyr::ldply(results, function(el){el}) %>% 
    reshape2::melt(value.name = "risk", variable.name = "lambda") %>% 
    dplyr::group_by(lambda) %>% dplyr::summarise(risk = mean(risk)) %>% dplyr::ungroup()
  return(out)
}

##' unbiasedly estimates risk using unbiased estimator of df
##' 
##' @param noised noised image
##' @param lambda_seq sequence of lambda values to be used for fused lasso
##' @param s2 provided variance of the noise to used in the risk estimation
##' @return data frame with columns lambda, df = unbiased estimate of df, risk = unbiased estimate of risk
unbiased_risk <- function(noised, lambda_seq, s2 = NULL, sol = NULL){
  
  if(is.null(sol)){
    #t0 <- Sys.time()
    sol <- flsa(noised, lambda1 = 0, lambda2 = lambda_seq)
    #Sys.time()-t0
  }
  
  ## TODO change this to chosing lambda that makes an elbow in training error!
  if(is.null(s2)){
    ## estimating s2
    s2 <- sd(noised - sol[5,,])
  }
  
  df <- numeric(length(lambda_seq))
  risk <- numeric(length(lambda_seq))
  for(i in 1:dim(sol)[1]){
    df[i] <- count_df_image(sol[i,,])
    risk[i] <- mean((noised - sol[i,,])^2) + 2*s2*df[i]/(dim(sol)[2]*dim(sol)[3])
  }
  return(data.frame(lambda = lambda_seq, df = df, risk = risk))
}

##' estimates risk using CB estimator
##' 
##' @param noised noised image
##' @param lambda_seq sequence of lambda values to be used for fused lasso
##' @param alpha positive number, noise size
##' @param B number of noise replications
##' @param s2 provided variance of the noise to used in the risk estimation
##' @return data frame with columns lambda, df = unbiased estimate of df, risk = unbiased estimate of risk
CB_risk <- function(noised, lambda_seq, alpha, B, s2, sol = NULL){
  
  #progress bar
  iterations = B
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- makeCluster(detectCores())
  registerDoSNOW(cl)
  results <- foreach(i = 1:B, .packages=c("flsa"), .options.snow = opts) %dopar% {
    #cat(paste("i = ", i, "\n"))
    omega <- array(rnorm(nrow(noised)*ncol(noised), sd=sqrt(s2)), dim(noised))
    noised_alpha_plus <- noised + sqrt(alpha)*omega
    noised_alpha_minus <- noised - (1/sqrt(alpha))*omega
    #t0 <- Sys.time()
    sol <- flsa(noised_alpha_plus, lambda1 = 0, lambda2 = lambda_seq)
    #Sys.time()-t0
    out <- apply(sol, 1, function(denoised){
      mean((noised_alpha_minus - denoised)^2) - s2 - mean(omega^2)/alpha}) 
  }
  stopCluster(cl)
  rm(cl)
  
  return(reshape2::melt(colMeans(plyr::laply(results, function(el){el}))))
  
}