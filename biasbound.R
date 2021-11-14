#########
## bias bound simulation
#########

root <- rprojroot::has_file(".git/index")
figsdir = root$find_file("figures")
resultsdir = root$find_file("savedfiles")
codedir = root$find_file("")

source(root$find_file("helpers_and_generators.R"))
source(root$find_file("estimators.R"))
source(root$find_file("gs_for_fitting.R"))
source(root$find_file("true_gs.R"))

#' computes bias and bias bound
#' 
#' @param X matrix of features
#' @param alphavec vector of values of alpha 
#' @param sig2 variance
#' @param g_true true mean function that takes X and outputs mu
#' @param g function that is being evaluated 
#' @param nrep number of MC approximations
#' @return data frame with real Risk_alpha, bound, and slope-bound for each value of alpha from alphavec 
bias_and_bound <- function(X, alphavec, sig2, g_true, g, nrep){
  n <- nrow(X)
  y <- replicate(nrep, gen_y_gaussian(X, sig2, g_true))
  if(n == 1) y <- matrix(y, 1, nrep)
  w <- matrix(rnorm(n*nrep, 0, sqrt(sig2)), nrow = n, ncol = nrep)
  muvec <- g_true(X)
  
  out <- data.frame(alpha = alphavec,
                    riskalpha = 0,
                    slope = 0,
                    bound = 0)
  for(i in 1:length(alphavec)){
    cat(paste(i, "/", length(alphavec), "alphas...\n", sep =))
    yplus <- y + sqrt(alphavec[i])*w
    gyplus <- apply(yplus, 2, g, X)
    if(n == 1) gyplus <- matrix(gyplus, 1, nrep)
    aux <- apply(gyplus, 2, function(gyplusi){mean((gyplusi - muvec)^2)})
    out$riskalpha[i] <- mean(aux)
    out$slope[i] <- sqrt(n)*sd(aux)/sqrt(2)
    out$bound[i] <- out$slope[i]*alphavec[i]
  }
  return(out)
}

## setup
n <- 100
p <- 200
g_true <- g_true_lin_sparse_5
gen_y <- gen_y_gaussian
rho <- 0
set.seed(202)
X <- gen_x(n,p,rho,5)
SNR <- 0.4
sig2 <- var(g_true(X))/SNR
alpha_seq <- c(0, 0.005, 0.01, 0.03, 0.05, 0.08, 0.1)

#### approximating bias and bias bound for FS with maximum of 90 steps
g <- g_forwardstep_90steps
filename <- "FS_90"
set.seed(202)
t0 <- Sys.time()
real <- bias_and_bound(X, alpha_seq, sig2, g_true, g, nrep = 1000)
Sys.time() - t0
real$bias <- real$riskalpha - real$riskalpha[1]
saveRDS(real, paste(resultsdir, "/bias_", filename, sep = ""))


#### approximating bias and bias bound for FS with maximum of 3 steps
g <- g_forwardstep_3steps
filename <- "FS_3"
set.seed(202)
t0 <- Sys.time()
real <- bias_and_bound(X, alpha_seq, sig2, g_true, g, nrep = 1000)
Sys.time() - t0
real$bias <- real$riskalpha - real$riskalpha[1]
saveRDS(real, paste(resultsdir, "/bias_", filename, sep = ""))

#### approximating bias and bias bound for FS with maximum of 10 steps
g <- g_forwardstep_10steps
filename <- "FS_10"
set.seed(202)
t0 <- Sys.time()
real <- bias_and_bound(X, alpha_seq, sig2, g_true, g, nrep = 1000)
Sys.time() - t0
real$bias <- real$riskalpha - real$riskalpha[1]
saveRDS(real, paste(resultsdir, "/bias_", filename, sep = ""))


# plots with all of them
filename <- "FS_90"
real_90 <- readRDS(paste(resultsdir, "/bias_", filename, sep = ""))
filename <- "FS_10"
real_10 <- readRDS(paste(resultsdir, "/bias_", filename, sep = ""))
filename <- "FS_3"
real_3 <- readRDS(paste(resultsdir, "/bias_", filename, sep = ""))

real_90$steps = "90"
real_10$steps = "10"
real_3$steps = "3"

df_plot <- bind_rows(real_90, bind_rows(real_10, real_3)) %>% 
  reshape2::melt(measure.vars = c("bound", "bias")) 
df_plot$steps <- factor(df_plot$steps, levels=c("3", "10", "90"))
df_plot$variable <- factor(df_plot$variable, levels=c("bias", "bound"))

df_plot$textsteps <- paste("k = ", df_plot$steps, sep = "")
df_plot$textvariable <- ifelse(df_plot$variable == "bound", "Bound", "Bias")
p1 <- df_plot %>% mutate(across(textsteps, factor, levels=c("k = 3","k = 10","k = 90"))) %>% 
  ggplot(aes(x = alpha, y = value, color = textsteps)) +
  geom_point() +
  theme_minimal(base_size = 16) +
  guides(color = "none") +
  geom_line(aes(linetype = textvariable)) +
  scale_alpha(guide = 'none') +
  xlab(expression(alpha)) + 
  facet_wrap(~textsteps, nrow = 1) +
  ylab("") + 
  labs(linetype="") + theme(panel.spacing = unit(2, "lines"))
p1
ggsave(paste(figsdir,"/FINALplot_04_bias.pdf", sep = ""), plot = p1, device = "pdf", dpi = "retina", width = 12, height = 4.5)

