##############################
######## different mean estimators to be evaluated and true mean functions
##############################

#' fits linear regression model and predicts Y
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_lin_reg <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  mod <- lm(y ~ xtrain)
  if(is.vector(xtest)){
    g <- sum(c(1, xtest)*mod$coef)
  } else{
    g <- as.vector(cbind(1, xtest)%*%mod$coef)
  }
  return(g)
}

#' fits step AIC model selection process and predicts Y using selected model
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_stepAIC <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  df_train <- data.frame(xtrain)
  df_train$y <- y
  modFull <- lm(y ~ ., data = df_train)
  mod <- stepAIC(modFull, direction = "both", 
                 trace = FALSE)
  if(is.vector(xtest)){
    xtest <- rbind(xtest, rep(1, length(xtest)))
    df_test <- data.frame(xtest)
    g <- predict(mod, df_test)[1]
  } else{
    df_test <- data.frame(xtest)
    g <- predict(mod, df_test)
  }
  return(g)
}

#' fits ridge regression with fixed lambda = 6 and predicts Y
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_ridge06 <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  lambda = .6
  xtrain <- scale(xtrain)
  xtrain <- cbind(1, xtrain)
  beta <- solve(t(xtrain)%*%xtrain + diag(lambda*2*ncol(xtrain), ncol(xtrain)))%*%t(xtrain)%*%y
  if(is.vector(xtest)){
    return(as.numeric(c(1,xtest)%*%as.vector(beta)))
  } else{
    return(as.vector(xtrain%*%beta))
  }
}


#' LASSO model whose lambda is tuned using CV and predicts Y using selected model
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_lasso <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  mod <- cv.glmnet(xtrain, y, nfolds = 5)
  if(is.vector(xtest)){
    return(as.vector(predict(mod, matrix(xtest, nrow=1), "lambda.min")))
  } else{
    return(as.vector(predict(mod, xtest, "lambda.min")))
  }
}

#' forward stepwise model with fixed number of steps k = 90
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_forwardstep_90steps <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  mod <- fs(xtrain, y, maxsteps = 90)
  if(is.vector(xtest)){
    return(as.vector(predict(mod, xtest)[,90]))
  } else{
    return(as.vector(predict(mod, xtest)[,90]))
  }
}

#' forward stepwise model with fixed number of steps k = 10
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_forwardstep_10steps <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  mod <- fs(xtrain, y, maxsteps = 10)
  if(is.vector(xtest)){
    return(as.vector(predict(mod, xtest)[,10+1]))
  } else{
    return(as.vector(predict(mod, xtest)[,10+1]))
  }
}

#' forward stepwise model with fixed number of steps k = 3
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_forwardstep_3steps <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  mod <- fs(xtrain, y, maxsteps = 3)
  if(is.vector(xtest)){
    return(as.vector(predict(mod, xtest)[,3+1]))
  } else{
    return(as.vector(predict(mod, xtest)[,3+1]))
  }
}

#' forward stepwise model with fixed number of steps k = 2
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_forwardstep_2steps <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  mod <- fs(xtrain, y, maxsteps = 2)
  if(is.vector(xtest)){
    return(as.vector(predict(mod, xtest)[,2+1]))
  } else{
    return(as.vector(predict(mod, xtest)[,2+1]))
  }
}

#' relaxed LASSO model whose lambda is tuned using CV and predicts Y using selected model
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_relaxed_lasso <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  mod <- cv.glmnet(xtrain, y, nfolds = 5)
  lambmin <- mod$lambda.min
  mod_rel <- glmnet(xtrain, y, lambda = mod$lambda.min, relaxed = TRUE)
  if(is.vector(xtest)){
    return(as.vector(predict(mod_rel, matrix(xtest, nrow=1))))
  } else{
    return(as.vector(predict(mod_rel, xtest)))
  }
}

#' LASSO with fixed lambda = 0.31
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_lasso_fixedlambda <- function(y, xtrain, xtest = NULL, lambda = 0.31){
  if(is.null(xtest)){xtest = xtrain}
  mod <- glmnet(xtrain, y, lambda = lambda)
  if(is.vector(xtest)){
    return(as.vector(predict(mod, matrix(xtest, nrow=1), "lambda.min")))
  } else{
    return(as.vector(predict(mod, xtest, "lambda.min")))
  }
}

#' xgboost with fixed parameters
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_xgboost <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  dtrain <- xgb.DMatrix(xtrain, 
                        label = y)
  dtest <- xgb.DMatrix(xtest, 
                       label = y)
  
  params <- list(booster = "gbtree", 
                 eta = 0.05, 
                 gamma = 3, 
                 max_depth = 6,
                 objective = "reg:squarederror")
  mod <- xgb.train(params = params, 
                   nrounds = 50,
                   nthread = 1,
                   data = dtrain, 
                   colsample_bytree = 0.8,
                   subsample = 0.8,
                   maximize = FALSE)
  pred <- predict(mod, newdata = dtest)
  return(pred)
}

#' fits regression tree and predicts Y, default tree hyperparameters
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_tree <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  df <- data.frame(xtrain)
  df$y <- y
  mod <- rpart(y~., data = df, method='anova')
  if(is.vector(xtest)){
    return(as.vector(predict(mod, newdata = data.frame(matrix(xtest, nrow = 1, byrow = T))   )))
  } else{
    return(as.vector(predict(mod, data.frame(xtest))))
  }
}

#' fits regression tree and predicts Y, overfitting
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_tree_overfitting <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  df <- data.frame(xtrain)
  df$y <- y
  mod <- rpart(y~., data = df, method='anova', control = rpart.control(minsplit = 5, minbucket = 1, 
                                                                       cp = 0.0005, maxdepth = 15)) 
  if(is.vector(xtest)){
    return(as.vector(predict(mod, newdata = data.frame(matrix(xtest, nrow = 1, byrow = T))   )))
  } else{
    return(as.vector(predict(mod, data.frame(xtest))))
  }
}

#' fits regression tree and predicts Y, underfitting
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_tree_underfitting <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  df <- data.frame(xtrain)
  df$y <- y
  mod <- rpart(y~., data = df, method='anova', control = rpart.control(minsplit = 5, minbucket = 5, 
                                                                       cp = 0.001, maxdepth = 2)) 
  if(is.vector(xtest)){
    return(as.vector(predict(mod, newdata = data.frame(matrix(xtest, nrow = 1, byrow = T))   )))
  } else{
    return(as.vector(predict(mod, data.frame(xtest))))
  }
}


#' fits GAM and predicts Y
#' 
#' @param y response vector
#' @param xtrain matrix of features (fixed-X setting)
#' @param xtest matrix of features for testing. if not provided, uses xtrain
#' @return vector of \hat Y
g_gamsplines <- function(y, xtrain, xtest = NULL){
  if(is.null(xtest)){xtest = xtrain}
  
  p <- ncol(xtrain)
  df <- data.frame(y = y, xtrain)
  mod <- suppressWarnings(gam(as.formula(paste0("y~",paste0("s(X",1:p,")",collapse="+"))), data = df))
  
  Xtest <- data.frame(xtest)
  
  return(as.vector(predict(mod, Xtest)))
}
