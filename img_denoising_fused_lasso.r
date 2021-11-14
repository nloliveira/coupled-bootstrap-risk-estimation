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


set.seed(208)
data(lena128)
img <- normalize(lena128)
image(img, col = gray.colors(500))

## adding noise
s2 = .01
s2_est = s2
noised <- normalize(add_gauss_noise(img, sd = sqrt(s2)))
image(noised, col = gray.colors(500))

## sequence of lambdas for fused lasso
lambda_seq <- c(0,exp(seq(log(0.001),log(1.5),length.out = 19)))

## fused lasso
t0 <- Sys.time()
sol <- flsa(noised, lambda1 = 0, lambda2 = lambda_seq)
Sys.time() - t0
saveRDS(sol, paste(resultsdir, "/sol", sep = ""))

## real risk
set.seed(57)
t0 <- Sys.time()
real_risk <- true_risk(img, lambda_seq, s2)
Sys.time() - t0
saveRDS(real_risk, paste(resultsdir, "/real_risk0_fused_lasso", sep = ""))

## unbiased DF estimator
t0 <- Sys.time()
fusedlasso_unbiasedDF <- unbiased_risk(noised, lambda_seq, s2 = s2_est, sol)
Sys.time() - t0
saveRDS(fusedlasso_unbiasedDF, paste(resultsdir, "/fusedlasso_unbiasedDF", sep = ""))

## CB, alpha = 0.1
set.seed(57)
alpha = 0.1
B = 30
t0 <- Sys.time()
CB_est_risk_01 <- CB_risk(noised, lambda_seq, alpha, B, s2 = s2_est, sol)
Sys.time() - t0
saveRDS(CB_est_risk_01, paste(resultsdir, "/CB_est_risk_fused_lasso_01", sep = ""))

## CB, alpha = 0.3
set.seed(57)
alpha = 0.3
B = 15
t0 <- Sys.time()
CB_est_risk_03 <- CB_risk(noised, lambda_seq, alpha, B, s2 = s2_est, sol)
Sys.time() - t0
saveRDS(CB_est_risk_03, paste(resultsdir, "/CB_est_risk_fused_lasso_03", sep = ""))

## CB, alpha = 0.5
set.seed(57)
alpha = 0.5
B = 5
t0 <- Sys.time()
CB_est_risk_05 <- CB_risk(noised, lambda_seq, alpha, B, s2 = s2_est, sol)
Sys.time() - t0
saveRDS(CB_est_risk_05, paste(resultsdir, "/CB_est_risk_fused_lasso_05", sep = ""))


### plotting
sol <- readRDS(paste(resultsdir, "/sol", sep = ""))
real_risk <- readRDS(paste(resultsdir, "/real_risk0_fused_lasso", sep = ""))
fusedlasso_unbiasedDF <- readRDS(paste(resultsdir, "/fusedlasso_unbiasedDF", sep = ""))
CB_est_risk_01 <- readRDS(paste(resultsdir, "/CB_est_risk_fused_lasso_01", sep = ""))
CB_est_risk_03 <- readRDS(paste(resultsdir, "/CB_est_risk_fused_lasso_03", sep = ""))
CB_est_risk_05 <- readRDS(paste(resultsdir, "/CB_est_risk_fused_lasso_05", sep = ""))


CB_est_risk_01$lambda <- as.numeric(rownames(CB_est_risk_01))
CB_est_risk_03$lambda <- as.numeric(rownames(CB_est_risk_03))
CB_est_risk_05$lambda <- as.numeric(rownames(CB_est_risk_05))
CB_est_risk_01$risk <- CB_est_risk_01$value
CB_est_risk_03$risk <- CB_est_risk_03$value
CB_est_risk_05$risk <- CB_est_risk_05$value

CB_est_risk_01$lambda[which.min(CB_est_risk_01$value)]
CB_est_risk_03$lambda[which.min(CB_est_risk_03$value)]
CB_est_risk_05$lambda[which.min(CB_est_risk_05$value)]
fusedlasso_unbiasedDF$lambda[which.min(fusedlasso_unbiasedDF$risk)]
real_risk$lambda[which.min(real_risk$risk)]

CB_est_risk_01$which <- "0.1"
CB_est_risk_03$which <- "0.3"
CB_est_risk_05$which <- "0.5"
fusedlasso_unbiasedDF$which <- "SURE"
real_risk$which <- "Risk"

real_risk$lambda <- as.numeric(as.vector(real_risk$lambda))

df_plot <- dplyr::bind_rows(CB_est_risk_01, CB_est_risk_03, CB_est_risk_05, 
                            fusedlasso_unbiasedDF, real_risk) %>% dplyr::select(lambda, risk, which)

p1 <- ggplot(df_plot, aes(x = log(lambda), y = risk, group = which, color = which)) + 
  geom_line(aes(linetype = which, size = which)) + 
  theme_minimal(base_size = 18) + 
  theme(legend.title = element_blank(), legend.position = "right") +
  scale_color_manual(values=c("gray78", "gray35", "black", "dodgerblue", "firebrick3"), 
                     labels = expression(paste(alpha, "= 0.1"), paste(alpha, "= 0.3"), paste(alpha, "= 0.5"), "Risk", "SURE")) + 
  scale_linetype_manual(values=c("solid", "solid", "solid", "dotted", "longdash"), 
                        labels = expression(paste(alpha, "= 0.1"), paste(alpha, "= 0.3"), paste(alpha, "= 0.5"), "Risk", "SURE")) + 
  scale_size_manual(values=c(0.6,0.6, 0.6, 1,0.8), 
                    labels = expression(paste(alpha, "= 0.1"), paste(alpha, "= 0.3"), paste(alpha, "= 0.5"), "Risk", "SURE")) + 
  xlab(expression(paste("log ", lambda))) + 
  ylab("Risk")
p1
ggsave(paste(figsdir, "/risk_lena_CB.pdf", sep = ""), 
       plot = p1, device = "pdf", width = 10, height = 6)

CB_est_risk <- CB_est_risk_01

pdf(file = paste(figsdir, "/imagedenoising_image.pdf", sep = ""),  
    width = 14, height = 4) 
par(mfrow = c(1,4))
image(img, col = gray.colors(500), main = "Original")
image(noised, col = gray.colors(500), main = "Noisy")
image(sol[which.min(CB_est_risk$value),,], col = gray.colors(500), main = "CB-denoised")
image(sol[which.min(fusedlasso_unbiasedDF$risk),,], col = gray.colors(500), main = "SURE-denoised")
dev.off()
par(mfrow=c(1,1))

