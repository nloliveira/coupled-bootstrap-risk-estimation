######################################
#######  image denoising using fused LASSO
#######  compare unbiased error estimate and aux rand estimator
######################################
root <- rprojroot::has_file(".git/index")
figsdir = root$find_file("figures")
resultsdir = root$find_file("savedfiles")
codedir = root$find_file("")

source(root$find_file("helpers_and_generators.R"))
source(root$find_file("estimators.R"))
source(root$find_file("gs_for_fitting.R"))
source(root$find_file("true_gs.R"))
source(root$find_file("helpers_image_denoising.R"))


library(imager)
library(flsa)

set.seed(208)
fpath <- system.file('extdata/parrots.png',package='imager')
parrots <- load.image(fpath)
plot(grayscale(parrots))
img <- grayscale(parrots)[,,1,1] 
image(img[72:227, 172:311])
img <- img[72:227, 172:311]
img <- img[,ncol(img):1]
img <- normalize(img)
image(img, col = gray.colors(500))

## adding noise
set.seed(208)
s2 = .03
s2_est = s2
noised <- normalize(add_gauss_noise(img, sd = sqrt(s2)))
image(noised, col = gray.colors(500))

##############
### Many noise replications
##############
# nnoise <- 10
# lambda_seq <- c(exp(seq(log(0.007),log(2),length.out = 50)))
# 
# init_seed <- 202
# for(i in 1:nnoise){
#   set.seed((init_seed + i))
#   cat(paste("\n noise ", i, "...\n\n"))
#   cat(paste("\t seed: ", init_seed + i, "\n"))
#   print(Sys.time())
#   noised <- normalize(add_gauss_noise(img, sd = sqrt(s2)))
#   saveRDS(noised, paste(resultsdir, "noised_differentnoise_1it", (init_seed + i), "seed.RDS", sep = ""))
#   cat("Done with image noise generation\n")
#   t0 <- Sys.time()
#   sol <- flsa::flsa(noised, lambda1 = 0, lambda2 = lambda_seq)
#   print(Sys.time() - t0)
#   saveRDS(sol, paste(resultsdir, "sol_differentnoise_1it", (init_seed + i), "seed.RDS", sep = ""))
#   cat("Done with model fit\n")
#   ## unbiased DF estimator
#   t0 <- Sys.time()
#   fusedlasso_unbiasedDF_differentnoise <- unbiased_risk(noised, lambda_seq, s2 = s2_est, sol)
#   print(Sys.time() - t0)
#   saveRDS(fusedlasso_unbiasedDF_differentnoise, paste(resultsdir, "fusedlasso_unbiasedDF_differentnoise_1it", (init_seed + i), "seed.RDS", sep = ""))
#   cat("Done with SURE\n")
#   ## CB, alpha = 0.1
#   alpha = 0.1
#   B = 30
#   t0 <- Sys.time()
#   CB_est_risk_01_differentnoise <- CB_risk(noised, lambda_seq, alpha, B, s2 = s2_est, sol)
#   print(Sys.time() - t0)
#   saveRDS(CB_est_risk_01_differentnoise, paste(resultsdir, "CB_est_risk_fused_lasso_01_differentnoise_1it", (init_seed + i), "seed.RDS", sep = ""))
#   cat("Done with CB 0.1\n")
#   ## CB, alpha = 0.3
#   alpha = 0.3
#   B = 15
#   t0 <- Sys.time()
#   CB_est_risk_03_differentnoise <- CB_risk(noised, lambda_seq, alpha, B, s2 = s2_est, sol)
#   print(Sys.time() - t0)
#   saveRDS(CB_est_risk_03_differentnoise, paste(resultsdir, "CB_est_risk_fused_lasso_03_differentnoise_1it", (init_seed + i), "seed.RDS", sep = ""))
#   cat("Done with CB 0.3\n")
#   # CB, alpha = 0.5
#   alpha = 0.5
#   B = 5
#   t0 <- Sys.time()
#   CB_est_risk_05_differentnoise <- CB_risk(noised, lambda_seq, alpha, B, s2 = s2_est, sol)
#   print(Sys.time() - t0)
#   saveRDS(CB_est_risk_05_differentnoise, paste(resultsdir, "CB_est_risk_fused_lasso_05_differentnoise_1it", (init_seed + i), "seed.RDS", sep = ""))
#   cat("Done with CB 0.5\n")
#   cat(paste("Saved results for iteration ", i, "...\nStarting next iteration.\n"))
# }

## plotting
noised <- list()
sol <- list()
fusedlasso_unbiasedDF <- list()
CB_est_risk_fused_lasso_01 <- list()
CB_est_risk_fused_lasso_03 <- list()
CB_est_risk_fused_lasso_05 <- list()
for(i in 1:10){
  start_seed = 202
  cat(i, "\n")
  noised[[i]] <- readRDS(paste0(resultsdir, "/noised_differentnoise_1it", start_seed+i, "seed.RDS"))
  sol[[i]] <- readRDS(paste0(resultsdir, "/sol_differentnoise_1it", start_seed+i, "seed.RDS"))
  fusedlasso_unbiasedDF[[i]] <- readRDS(paste0(resultsdir, "/fusedlasso_unbiasedDF_differentnoise_1it", start_seed+i, "seed.RDS"))
  CB_est_risk_fused_lasso_01[[i]] <- readRDS(paste0(resultsdir, "/CB_est_risk_fused_lasso_01_differentnoise_1it", start_seed+i, "seed.RDS"))
  CB_est_risk_fused_lasso_03[[i]] <- readRDS(paste0(resultsdir, "/CB_est_risk_fused_lasso_03_differentnoise_1it", start_seed+i, "seed.RDS"))
  CB_est_risk_fused_lasso_05[[i]] <- readRDS(paste0(resultsdir, "/CB_est_risk_fused_lasso_05_differentnoise_1it", start_seed+i, "seed.RDS"))
}

get_optm_lambda_SURE <- function(dat){
  lambda <- dat[,1]
  risk  <- dat[,3]
  return(lambda[which.min(risk)])
}

get_optm_lambda_CB <- function(dat){
  lambda <- as.numeric(rownames(dat))
  risk  <- dat[,1]
  return(lambda[which.min(risk)])
}

optimal_lambdas <- data.frame(
  SURE = unlist(lapply(fusedlasso_unbiasedDF, get_optm_lambda_SURE)),
  CB01 = unlist(lapply(CB_est_risk_fused_lasso_01, get_optm_lambda_CB)),
  CB03 = unlist(lapply(CB_est_risk_fused_lasso_03, get_optm_lambda_CB)),
  CB05 = unlist(lapply(CB_est_risk_fused_lasso_05, get_optm_lambda_CB)))

optimal_lambdas %>%
  reshape2::melt() %>%
  mutate(value = as.factor(value)) %>%
  group_by(variable, value) %>%
  count(variable, as.factor(value), .drop = FALSE)

dat <- optimal_lambdas %>% reshape2::melt() %>%
  mutate(variable = as.character(variable),
         value = round(log(value),3))

p1 <- dat %>% ggplot() +
  geom_bar(aes(value, group = variable, fill = variable),
           stat = "count") +
  scale_fill_discrete(labels = c(expression(paste(alpha, " = 0.1")),
                                  expression(paste(alpha, " = 0.3")),
                                  expression(paste(alpha, " = 0.5")),
                                  "SURE")) +
  theme_minimal(base_size = 20) +
  theme(legend.title = element_blank()) +
  ylab("Count") +
  xlab(expression(log(lambda))) +
  xlim(range(log(fusedlasso_unbiasedDF[[1]][,1])))
p1

ggsave(paste(figsdir, "/selected_lambdas_multiple_noise.pdf", sep = ""), 
       plot = p1, device = "pdf", width = 10, height = 6)

####
pdf(file = paste(figsdir, "/imagedenoising_image_multiplelambdas.pdf", sep = ""),  
    width = 11, height = 7) 
par(mfrow = c(2,3))
image(noised[[1]], col = gray.colors(500), main = "Noise iteration 1")
image(sol[[1]][27,,], col = gray.colors(500), main = "SURE and CB0.1, lambda = 0.141")
image(sol[[1]][29,,], col = gray.colors(500), main = "CB0.5, lambda = 0.177")
image(noised[[4]], col = gray.colors(500), main = "Noise iteration 4")
image(sol[[4]][28,,], col = gray.colors(500), main = "SURE, lambda = 0.158")
image(sol[[4]][29,,], col = gray.colors(500), main = "CB0.5, lambda = 0.177")
dev.off()
par(mfrow=c(1,1))

p1 <- bind_rows(
  plyr::ldply(fusedlasso_unbiasedDF, function(el){
    data.frame(lambda = el$lambda, riskhat = el$risk, method = "SURE", id = rnorm(1))
  }),
  plyr::ldply(CB_est_risk_fused_lasso_03, function(el){
    data.frame(lambda = as.numeric(rownames(el)),
               riskhat = el[,1], method = "CB01", id = rnorm(1))
  }),
  plyr::ldply(CB_est_risk_fused_lasso_01, function(el){
    data.frame(lambda = as.numeric(rownames(el)),
               riskhat = el[,1], method = "CB03", id = rnorm(1))
  }),
  plyr::ldply(CB_est_risk_fused_lasso_05, function(el){
    data.frame(lambda = as.numeric(rownames(el)),
               riskhat = el[,1], method = "CB05", id = rnorm(1))
  })) %>% 
  ggplot(aes(x = log(lambda), y = riskhat, color = method, group = id)) +
  geom_line(alpha = 0.3) +
  theme_minimal(base_size = 20) +
  ylab(expression(hat(Risk))) +
  xlab(expression(log(lambda))) +
  theme(legend.title = element_blank())
p1
ggsave(paste(figsdir, "/risk_curves_multiple_noise.pdf", sep = ""), 
       plot = p1, device = "pdf", width = 10, height = 6)




# ---------
# results for one noise iteration simulation
df_plot <- bind_rows(
  plyr::ldply(list(fusedlasso_unbiasedDF[[1]]), function(el){
    data.frame(lambda = el$lambda, riskhat = el$risk, method = "SURE", id = rnorm(1))
  }),
  plyr::ldply(list(CB_est_risk_fused_lasso_01[[1]]), function(el){
    data.frame(lambda = as.numeric(rownames(el)),
               riskhat = el[,1], method = "CB01", id = rnorm(1))
  }),
  plyr::ldply(list(CB_est_risk_fused_lasso_03[[1]]), function(el){
    data.frame(lambda = as.numeric(rownames(el)),
               riskhat = el[,1], method = "CB03", id = rnorm(1))
  }),
  plyr::ldply(list(CB_est_risk_fused_lasso_05[[1]]), function(el){
    data.frame(lambda = as.numeric(rownames(el)),
               riskhat = el[,1], method = "CB05", id = rnorm(1))
  }))

p1 <- df_plot %>% ggplot() +
  geom_line(aes(x = log(lambda), y = riskhat, color = method, group = method),
            linewidth = 1) +
  scale_color_discrete(labels = c(expression(paste(alpha, " = 0.1")),
                                  expression(paste(alpha, " = 0.3")),
                                  expression(paste(alpha, " = 0.5")),
                                  "SURE")) +
  theme_minimal(base_size = 20) +
  ylab(expression(hat(Risk))) +
  xlab(expression(log(lambda))) +
  theme(legend.title = element_blank()) 

p1
ggsave(paste(figsdir, "/risk_parrots_CB_onereplicate.pdf", sep = ""), 
       plot = p1, device = "pdf", width = 10, height = 6)


pdf(file = paste(figsdir, "/imagedenoising_sidebyside_onereplicate.pdf", sep = ""),  
    width = 10, height = 3) 
par(mfrow = c(1,4))
image(img, col = gray.colors(500), main = "Original")
image(noised[[1]], col = gray.colors(500), main = "Noisy")
image(sol[[1]][29,,], col = gray.colors(500), main = "CB denoised")
image(sol[[1]][27,,], col = gray.colors(500), main = "SURE denoised")
dev.off()
par(mfrow=c(1,1))



