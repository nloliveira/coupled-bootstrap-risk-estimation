#########
## DF estimates for lasso path and forward stepwise
#########

root <- rprojroot::has_file(".git/index")
figsdir = root$find_file("figures")
resultsdir = root$find_file("savedfiles")
codedir = root$find_file("")

source(root$find_file("helpers_and_generators.R"))
source(root$find_file("estimators.R"))
source(root$find_file("gs_for_fitting.R"))
source(root$find_file("true_gs.R"))

# set-up
alpha_seq <- c(0.05, 0.1, 0.2, 0.5, 0.8, 1)
B <- 100
nsim <- 50
gen_y = gen_y_gaussian

## setup
n <- 100
p <- 200
g_true <- g_true_lin_sparse_5
rho <- 0
set.seed(202)
X <- gen_x(n,p,rho,5)
sig2 <- 1
var(g_true(X))/sig2 ##SNR

#### Lasso
g <- g_lasso_alllambdas
lambda = c(0.984351830, 0.939611520, 0.896904726, 0.856139022, 0.817226182, 0.780081991, 0.744626061, 0.710781657, 0.678475534, 0.647637774, 0.618201637, 0.590103419, 0.563282308,
           0.537680258, 0.513241860, 0.489914225, 0.467646868, 0.446391596, 0.426102409, 0.406735397, 0.388248645, 0.370602145, 0.353757706, 0.337678872, 0.322330846, 0.307680412,
           0.293695862, 0.280346932, 0.267604731, 0.255441683, 0.243831465, 0.232748948, 0.222170149, 0.212072173, 0.202433166, 0.193232266, 0.184449561, 0.176066044, 0.168063571,
           0.160424822, 0.153133266, 0.146173123, 0.139529329, 0.133187505, 0.127133928, 0.121355495, 0.115839700, 0.110574607, 0.105548820, 0.100751464, 0.096172154, 0.091800981,
           0.087628485, 0.083645635, 0.079843812, 0.076214787, 0.072750708, 0.069444076, 0.066287735, 0.063274855, 0.060398916, 0.057653692, 0.055033242, 0.052531897, 0.050144241,
           0.047865108, 0.045689565, 0.043612904, 0.041630630, 0.039738453, 0.037932279, 0.036208199, 0.034562480, 0.032991562, 0.031492045, 0.030060683, 0.028694378, 0.027390174,
           0.026145249, 0.024956907, 0.023822577, 0.022739804, 0.021706245, 0.020719663, 0.019777922, 0.018878985, 0.018020906, 0.017201828, 0.016419978, 0.015673665, 0.014961273,
           0.014281260, 0.013632154, 0.013012552, 0.012421111, 0.011856553, 0.011317654, 0.010803249, 0.010312225, 0.009843518)
filename <- "Lasso"
set.seed(202)

  y <- replicate(nsim, gen_y_gaussian(X, sig2, g_true))
  if(n == 1) y <- matrix(y, 1, nsim)
  pb <- txtProgressBar(max = nsim, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- makeCluster(detectCores()/4)
  registerDoSNOW(cl)
  t0 <- Sys.time()
  results <- foreach(i = 1:nsim, .packages=c("mvtnorm", "glmnet", "dplyr"), .options.snow = opts) %dopar% {
    temp <- g(y[,i], X, lambda)
    yhat <- temp$yhat
    support <- temp$size_support
    w <- matrix(rnorm(n*B, 0, sqrt(sig2)), nrow = n, ncol = B)
    out <- data.frame()
    for(i.alpha in 1:length(alpha_seq)){
      gyplus <- plyr::aaply(y[,i] + sqrt(alpha_seq[i.alpha])*w, 2, function(yplusi){g(yplusi, X, lambda)$yhat})
      for(j in 1:dim(gyplus)[3]){
        temp <- onlycov_uncentered(y[,i], yhat[,j], w, gyplus[,,j] %>% t(), alpha_seq[i.alpha], sig2, nw_seq = B)
        res_i <- data.frame(nw = B, method = "BY", alpha = alpha_seq[i.alpha],
                            cov = temp$cov[1]*(1+alpha_seq[i.alpha])/alpha_seq[i.alpha])
        res_i$df <- res_i$cov*n/(sig2*(1+alpha_seq[i.alpha]))
        temp <- riskAlpha(y[,i], w, gyplus[,,j] %>% t(), alpha_seq[i.alpha], sig2, nw_seq = B)
        temp$df <- temp$cov*n/(sig2*(1+alpha_seq[i.alpha]))
        res_i <- bind_rows(res_i, temp %>% select(nw, method, alpha, cov, df))
        res_i$lambda <- lambda[j]
        out <- bind_rows(out, res_i)
      }
    }
    list(out = out, support = support)
  }
  Sys.time() - t0
  stopCluster(cl)
  rm(cl)
  support <- plyr::ldply(results, function(el){el$support})
  results <- plyr::ldply(results, function(el){el$out})
  set.seed(57)
  t0 <- Sys.time()
  real <- real_values_X_alllambdas(X, c(0,alpha_seq), lambda, sig2, g_true, g, 1000, gen_y_gaussian)
  Sys.time()-t0
  saveRDS(list(results = results, real = real, support = support), paste(resultsdir, "/DF_", filename, sep = ""))




#### FS
g <- g_subset_allsteps
steps = 98
filename <- "FS"
set.seed(57)
t0 <- Sys.time()
{
  y <- replicate(nsim, gen_y_gaussian(X, sig2, g_true))
  if(n == 1) y <- matrix(y, 1, nsim)
  pb <- txtProgressBar(max = nsim, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- makeCluster(detectCores()/4)
  registerDoSNOW(cl)
  results <- foreach(i = 1:nsim, .packages=c("mvtnorm", "glmnet", "dplyr", "bestsubset"), .options.snow = opts, .errorhandling = 'remove') %dopar% {
    temp <- g(y[,i], X, steps)
    yhat <- temp$yhat
    support <- temp$size_support
    w <- matrix(rnorm(n*B, 0, sqrt(sig2)), nrow = n, ncol = B)
    out <- data.frame()
    for(i.alpha in 1:length(alpha_seq)){
      gyplus <- plyr::aaply(y[,i] + sqrt(alpha_seq[i.alpha])*w, 2, function(yplusi){g(yplusi, X, steps)$yhat})
      for(j in 1:dim(gyplus)[3]){
        temp <- onlycov_uncentered(y[,i], yhat[,j], w, gyplus[,,j] %>% t(), alpha_seq[i.alpha], sig2, nw_seq = B)
        res_i <- data.frame(nw = B, method = "BY", alpha = alpha_seq[i.alpha],
                            cov = temp$cov[1]*(1+alpha_seq[i.alpha])/alpha_seq[i.alpha])
        res_i$df <- res_i$cov*n/(sig2*(1+alpha_seq[i.alpha]))
        temp <- riskAlpha(y[,i], w, gyplus[,,j] %>% t(), alpha_seq[i.alpha], sig2, nw_seq = B)
        temp$df <- temp$cov*n/(sig2*(1+alpha_seq[i.alpha]))
        res_i <- bind_rows(res_i, temp)
        res_i$steps <- j-1
        out <- bind_rows(out, res_i)
      }
    }
    list(out = out, support = support)
  }
  stopCluster(cl)
  rm(cl)
  support <- plyr::ldply(results, function(el){el$support})
  results <- plyr::ldply(results, function(el){el$out})
  t0 <- Sys.time()
  real <- real_values_X_allsteps(X, c(0,alpha_seq), steps, sig2, g_true, g, 1000, gen_y_gaussian)
  Sys.time()-t0
  saveRDS(list(results = results, real = real, support = support), paste(resultsdir, "/DF_", filename, sep = ""))
}
Sys.time() - t0








### plotting
filename <- "Lasso"
temp <- readRDS(paste(resultsdir, "/DF_", filename, sep = ""))
resultslasso <- temp$results
reallasso <- temp$real
reallasso$lambda = as.numeric(as.character(reallasso$lambda))

filename <- "FS"
temp <- readRDS(paste(resultsdir, "/DF_", filename, sep = ""))
results <- temp$results
real <- temp$real
real$steps = as.numeric(as.character(real$steps))

df_forwardstep <- results %>% dplyr::group_by(alpha, steps, method) %>% dplyr::summarise(sd = sd(df), df = mean(df)) %>% 
  dplyr::left_join(real %>% dplyr::filter(alpha > 0) %>% dplyr::select(-df), c("alpha", "steps")) %>% dplyr::filter(steps >0)
df_forwardstep$g = "Forward stepwise"

df_lasso <- resultslasso %>% dplyr::group_by(alpha, lambda, method) %>% dplyr::summarise(sd = sd(df), df = mean(df)) %>% 
  dplyr::left_join(reallasso %>% dplyr::filter(alpha > 0) %>% dplyr::select(-df), c("alpha", "lambda"))
df_lasso$g = "Lasso"

### plotting against noised support size
bind_rows(df_lasso, df_forwardstep)  %>% dplyr::filter(alpha == 0.1) %>% 
  ggplot() +
  geom_ribbon(aes(x = support, ymin = df - sd, ymax = df+sd, fill = g), alpha = 0.4) +
  geom_line(aes(x = support, y = df), data = (real %>% filter(alpha == 0.1, support < 99)), lty = 2) +
  theme_minimal(base_size = 16) +
  facet_wrap(~method) +
  theme(legend.position = "bottom") +
  xlab("Support size") +
  ylab("Degrees of freedom") +
  scale_fill_manual(values=c("gray13", "gray62")) +
  geom_abline(intercept = 0, slope = 1, lty = 2)


### plotting against noiseless support size
df_forwardstep_01 <- results %>% dplyr::group_by(alpha, steps, method) %>% dplyr::summarise(sd = sd(df), df = mean(df)) %>% dplyr::filter(alpha == 0.1) %>% 
  dplyr::left_join(real %>% dplyr::filter(alpha == 0) %>% dplyr::select(-df), c("steps")) %>% dplyr::filter(steps >0)
df_forwardstep_01$g = "Forward stepwise"

df_lasso_01 <- resultslasso %>% dplyr::group_by(alpha, lambda, method) %>% dplyr::summarise(sd = sd(df), df = mean(df)) %>% dplyr::filter(alpha == 0.1) %>% 
  dplyr::left_join(reallasso %>% dplyr::filter(alpha == 0) %>% dplyr::select(-df), c("lambda"))
df_lasso_01$g = "Lasso"


p1 <- bind_rows(df_lasso_01, df_forwardstep_01)  %>% 
  ggplot() +
  geom_ribbon(aes(x = support, ymin = df - sd, ymax = df+sd, fill = g), alpha = 0.4) +
  geom_line(aes(x = support, y = df), data = (real %>% filter(alpha == 0, support < 99)), lty = 2) +
  theme_minimal(base_size = 16) +
  facet_wrap(~method) +
  theme(legend.position = "bottom") +
  xlab("Support size") +
  ylab("Degrees of freedom") +
  scale_fill_manual(values=c("gray13", "gray62")) +
  geom_line(aes(x = x, y = x), data = data.frame(x = 1:98, y = 1:98),  lty = 2)

p1
ggsave(paste(figsdir, "/FINALplotDF.pdf", sep = ""), plot = p1, device = "pdf", dpi = "retina", width = 10, height = 6.5)
