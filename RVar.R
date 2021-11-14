#########
## reducible variance simulation
#########

root <- rprojroot::has_file(".git/index")
figsdir = root$find_file("figures")
resultsdir = root$find_file("savedfiles")
codedir = root$find_file("")

source(root$find_file("helpers_and_generators.R"))
source(root$find_file("estimators.R"))
source(root$find_file("gs_for_fitting.R"))
source(root$find_file("true_gs.R"))

## setup 
n <- 100
p <- 200
g_true <- g_true_lin_sparse_5
rho <- 0
set.seed(202)
X <- gen_x(n,p,rho,5)
SNR <- 0.4
sig2 <- var(g_true(X))/SNR
g <- g_lasso_fixedlambda
gen_y <- gen_y_gaussian
filename <- "redvar_lassofixedlamb"

n_mc_var = n_mc_e = 500
alpha_seq <- seq(0.005, 0.1, length.out = 30)
B <- round(seq(1, 300, length.out = 100))

res <- data.frame(alpha = sort(rep(alpha_seq, length(B))), B = rep(B, length(alpha_seq)), redvar = NA)
set.seed(202)
for(i.a in 1:length(alpha_seq)){
  cat(paste("alpha = ", alpha_seq[i.a], "\n"))
  varw = rep(NA, n_mc_e)
  for(i.y in 1:n_mc_e){
    if(i.y %% 100 == 0) cat(paste("\t", i.y, "th y\n", sep = ""))
    y = gen_y_gaussian(X, sig2, g_true)
    varw[i.y] = var(replicate(n_mc_var, {
      w <- matrix(rnorm(n, 0, sqrt(sig2)), n, 1)
      gyplus <- matrix(g(y + sqrt(alpha_seq[i.a])*w, X), n, 1)
      riskAlpha(y, w, gyplus, alpha_seq[i.a], sig2, B_seq = 1)$risk
    }))
  }
  res[which(res["alpha"] == alpha_seq[i.a]),"redvar"] = mean(varw)
}
res$redvar = res$redvar/res$B
y <- replicate(nrep, gen_y(X, sig2, g_true))
gy <- apply(y, 2, g, X)
res$dominantterm = 4*sig2*mean(apply((y-gy), 2, function(yminusgy){mean(yminusgy^2)}))/(res$B*res$alpha)
saveRDS(res, paste(resultsdir, "/", filename, sep = ""))

## plotting
res <- readRDS(paste(resultsdir, "/", filename, sep = ""))
expr_terms = c("log true variance", "log dominant term")
names(expr_terms) <- c("redvar", "dominantterm")
p1 <- ggplot(res %>% reshape2::melt(measure.vars = c("redvar", "dominantterm")) %>% dplyr::filter(variable == "redvar")) + 
  geom_contour_fill(aes(x = alpha, y = B, z = log(value))) + 
  theme_minimal(base_size = 16)+ 
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  ylab("B") + 
  xlab(expression(alpha)) + 
  ggtitle("log RVar")
p1

p2 <- ggplot(res %>% reshape2::melt(measure.vars = c("redvar", "dominantterm")) %>% dplyr::filter(variable == "dominantterm")) + 
  geom_contour_fill(aes(x = alpha, y = B, z = log(value))) + 
  theme_minimal(base_size = 16)+ 
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  ylab("B") + 
  xlab(expression(alpha)) + 
  ggtitle("log dominant term")
p2
plt <- grid.arrange(p1, p2, ncol = 2)
ggsave(paste(figsdir,"/FINALplot_04_redvar.pdf", sep = ""), plot = plt, device = "pdf", dpi = "retina", width = 10, height = 6)
