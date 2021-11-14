#########
## methods comparison simulation
## multiple set-up combinations
## only ones in the paper are pushed to github
#########

root <- rprojroot::has_file(".git/index")
figsdir = root$find_file("figures")
resultsdir = root$find_file("savedfiles")
codedir = root$find_file("")

source(root$find_file("helpers_and_generators.R"))
source(root$find_file("estimators.R"))
source(root$find_file("gs_for_fitting.R"))
source(root$find_file("true_gs.R"))

## set-up
alpha_seq <- c(0.05, 0.1, 0.2, 0.5, 0.8, 1)
B <- 100
nsim <- 100
gen_y = gen_y_gaussian
nrep = 1000

## varying parameters
pars <- expand.grid(n = c(30, 100), p = c(10, 200), rho = c(0, 0.8), SNR = c(0.4, 2, 6),
                    true = c("linear", "linear5", "treenonlin"), g = c("underfittree", "fs2", "lasso31", "CVLasso", "ridge06"))
pars <- pars %>% dplyr::filter(!((n == 30) & (p == 200)))
pars$true <- as.character(pars$true)
pars$g <- as.character(pars$g)

for(i.sim in 1:nrow(pars)){
  if(pars$true[i.sim] == "linear") g_true = g_true_lin
  if(pars$true[i.sim] == "linear5") g_true = g_true_lin_sparse_5
  if(pars$true[i.sim] == "treenonlin") g_true = g_true_tree_notlin
  
  if(pars$g[i.sim] == "underfittree") g = g_tree_underfitting
  if(pars$g[i.sim] == "fs2") g = g_forwardstep_2steps
  if(pars$g[i.sim] == "lasso31") g = g_lasso_fixedlambda
  if(pars$g[i.sim] == "CVLasso") g = g_lasso
  if(pars$g[i.sim] == "ridge06") g = g_ridge06
  
  n = pars$n[i.sim]
  p = pars$p[i.sim]
  rho = pars$rho[i.sim]
  SNR = pars$SNR[i.sim]
  
  set.seed(202)
  X <- gen_x(n,p,rho,1)
  
  sig2 = var(g_true(X))/SNR
  
  filename = paste(gsub("[[:punct:]]", "", pars[i.sim,]), collapse = "_")
  
  cat(paste("\n\tIteration ", i.sim, "/", nrow(pars), ": ", filename, "\n", sep =))
  
  results <- NULL
  real <- NULL
  y <- replicate(nsim, gen_y_gaussian(X, sig2, g_true))
  if(n == 1) y <- matrix(y, 1, nsim)
  yhat <- apply(y, 2, g, X)
  if(n == 1) yhat <- matrix(yhat, 1, nsim)
  pb <- txtProgressBar(max = nsim, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- makeCluster(detectCores()/4)
  registerDoSNOW(cl)
  t0 <- Sys.time()
  results <- foreach(i = 1:nsim, .packages=c("mvtnorm", "glmnet", "dplyr", "rpart", "bestsubset"), .options.snow = opts) %dopar% {
    out <- data.frame()
    for(i.alpha in 1:length(alpha_seq)){
      w <- matrix(rnorm(n*B, 0, sqrt(sig2)), nrow = n, ncol = B)
      gyplus <- apply(y[,i] + sqrt(alpha_seq[i.alpha])*w, 2, g, X)
      temp <- onlycov_uncentered(y[,i], yhat[,i], w, gyplus, alpha_seq[i.alpha], sig2, nw_seq = B)
      res_i <- data.frame(nw = c(B,B), method = c("BY", "ParBoot"), alpha = rep(alpha_seq[i.alpha],2), 
                          training = c(temp$training[1], temp$training[1]), 
                          cov = temp$cov[1]*c(1/alpha_seq[i.alpha], 1), 
                          error_to_risk = sig2*c(1, 1))
      res_i$error <- res_i$training + 2*res_i$cov
      res_i$risk <- res_i$error - res_i$error_to_risk
      res_i$error_to_risk <- NULL
      out <- riskAlpha(y[,i], w, gyplus, alpha_seq[i.alpha], sig2, nw_seq = B) %>% 
        bind_rows(res_i) %>% bind_rows(out)
    }
    out
  }
  Sys.time() - t0
  stopCluster(cl)
  rm(cl)
  results <- plyr::ldply(results, function(el){el})
  real <- real_values_X(X, c(0,alpha_seq), sig2, g_true, g, nrep, gen_y_gaussian)
  saveRDS(list(results = results, real = real), paste(resultsdir, "/methods_comparison/results_", filename, sep = ""))
  
  p1 = results %>% dplyr::group_by(method, alpha) %>% dplyr::filter(method %in% c("BY", "CB")) %>% dplyr::summarise(avg = mean(risk), sd = sd(risk)) %>%
      ggplot(aes(x = factor(alpha), y = avg)) +
      geom_bar(aes(fill = method, color = method), alpha = 0.3, lwd = 1, stat="identity", position=position_dodge()) +
      geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd, color = method), width=.3, lwd = 1,
                    position=position_dodge()) +
      geom_hline(yintercept = real$riskalpha[1]) +
      geom_point(aes(x = factor(alpha), y = riskalpha), data = real %>% dplyr::filter(alpha != 0) %>% dplyr::select(alpha, riskalpha)) +
      theme_minimal(base_size = 20) +
      xlab(expression(alpha)) +
      ylab(expression(hat(Risk))) +
      theme(legend.position = "none", axis.text.x = element_text(angle=45))+
      facet_wrap(~method, ncol = 1) 
    
  ggsave(paste(figsdir, "/methods_comparison/plot_", filename, ".pdf", sep = ""), plot = p1, device = "pdf", dpi = "retina", width = 8, height = 10)
}



