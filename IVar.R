#################
### irreducible variance
### normal means problem set-up (no features)
#################

root <- rprojroot::has_file(".git/index")
figsdir = root$find_file("figures")
resultsdir = root$find_file("savedfiles")
codedir = root$find_file("")

source(root$find_file("helpers_and_generators.R"))
source(root$find_file("estimators.R"))
source(root$find_file("gs_for_fitting.R"))
source(root$find_file("true_gs.R"))

## set-up
n = 100
p = 200
rho = 0
SNR = 2
g_true = g_true_lin
n_mc_var = 300
n_mc_e = 300
g = g_lasso
set.seed(202)
X <- gen_x(n,p,rho,1)
sig2 = var(g_true(X))/SNR

alpha <- c(seq(0.01, 0.25, length.out = 15))
res_CB <- data.frame(method = "CB", alpha = alpha, irredvarT1 = NA, irredvarT2 = NA, irredvarCov = NA, redvar = NA, irredvar = NA, 
                     MCirredvarT1 = NA, MCirredvarT2 = NA, MCirredvarCov = NA, MCredvar = NA, MCirredvar = NA)
res_BY <- data.frame(method = "BY", alpha = alpha, irredvarT1 = NA, irredvarT2 = NA, irredvarCov = NA, redvar = NA, irredvar = NA, 
                     MCirredvarT1 = NA, MCirredvarT2 = NA, MCirredvarCov = NA, MCredvar = NA, MCirredvar = NA)
for(i.a in 1:length(alpha)){
  cat(paste("\nalpha = ", alpha[i.a], "\n"))
  varw_CB = rep(NA, n_mc_var)
  varw_BY = rep(NA, n_mc_var)
  expw_CB = rep(NA, n_mc_var)
  expw_BY = rep(NA, n_mc_var)
  expw_t1 = rep(NA, n_mc_var)
  expw_t2 = rep(NA, n_mc_var)
  expw_t1.2 = rep(NA, n_mc_var)
  for(i.y in 1:n_mc_var){
      cat(paste("\ti.y = ", i.y, "\n"))
    y = gen_y_gaussian(X, sig2, g_true)
    yhat = g(y, X)
    tmp <- data.frame(id = 1:n_mc_e, CB = NA, BY = NA, IVarT1 = NA, IVarT1.BY = NA, IVarT2 = NA)
    for(j in 1:n_mc_e){
      tmp[j,-1] <- riskAlpha_onew_ivarterms(y, yhat, g, sig2, alpha[i.a])
    }

    varw_CB[i.y] = var(tmp$CB)
    varw_BY[i.y] = var(tmp$BY)
    expw_CB[i.y] = mean(tmp$CB)
    expw_BY[i.y] = mean(tmp$BY)
    expw_t1[i.y] = mean(tmp$IVarT1)
    expw_t2[i.y] = mean(tmp$IVarT2)
    expw_t1.2[i.y] = mean(tmp$IVarT1.BY)
  }
  res_CB[i.a,"redvar"] = mean(varw_CB)
  res_CB[i.a,"MCredvar"] = var(varw_CB)/n_mc_var
  res_CB[i.a,"irredvar"] = var(expw_CB)
  res_CB[i.a,"MCirredvar"] = var((expw_CB - mean(expw_CB))^2)/n_mc_var
  res_CB[i.a,"irredvarT1"] = var(expw_t1)
  res_CB[i.a,"MCirredvarT1"] = var((expw_t1 - mean(expw_t1))^2)/n_mc_var
  res_CB[i.a,"irredvarT2"] = var(expw_t2)
  res_CB[i.a,"MCirredvarT2"] = var((expw_t2 - mean(expw_t2))^2)/n_mc_var
  res_CB[i.a,"irredvarCov"] = 2*cov(expw_t1, expw_t2)
  res_CB[i.a,"MCirredvarCov"] = 4*var((expw_t1 - mean(expw_t1))*(expw_t2 - mean(expw_t2)))/n_mc_var
  res_BY[i.a,"redvar"] = mean(varw_BY)
  res_BY[i.a,"MCredvar"] = var(varw_BY)/n_mc_var
  res_BY[i.a,"irredvar"] = var(expw_BY)
  res_BY[i.a,"MCirredvar"] = var((expw_BY - mean(expw_BY))^2)/n_mc_var
  res_BY[i.a,"irredvarT1"] = var(expw_t1.2)
  res_BY[i.a,"MCirredvarT1"] = var((expw_t1.2 - mean(expw_t1.2))^2)/n_mc_var
  res_BY[i.a,"irredvarT2"] = var(expw_t2)
  res_BY[i.a,"MCirredvarT2"] = var((expw_t2 - mean(expw_t2))^2)/n_mc_var
  res_BY[i.a,"irredvarCov"] = 2*cov(expw_t1.2, expw_t2)
  res_BY[i.a,"MCirredvarCov"] = 4*var((expw_t1.2 - mean(expw_t1.2))*(expw_t2 - mean(expw_t2)))/n_mc_var
}
res_CVLasso <- dplyr::bind_rows(res_CB, res_BY)
saveRDS(res_CVLasso, paste(resultsdir, "/IVarcomponents_CVLasso.RDS", sep = ""))

## plotting
res_CVLasso <- readRDS(paste(resultsdir, "/IVarcomponents_CVLasso.RDS", sep = ""))
df1 <- res_CVLasso %>% dplyr::select(method, alpha, irredvarT1, irredvarT2, irredvarCov, redvar) %>% 
  reshape2::melt(id.vars = c("alpha", "method"), value.name = "terms") 
df2 <- res_CVLasso %>% dplyr::select(method, alpha, MCirredvarT1, MCirredvarT2, MCirredvarCov, MCredvar) %>% 
  reshape2::melt(id.vars = c("alpha", "method"), value.name = "se")
df2$variable <- stringr::str_replace(df2$variable, "MC", "")
df2$se <- sqrt(df2$se)

df_plot <- dplyr::inner_join(df1, df2, by=c("method", "alpha", "variable"))
df_plot$variable <- ifelse(df_plot$variable == "irredvarT1", "IVar1", 
                           ifelse(df_plot$variable == "irredvarT2", "IVar2", 
                                  ifelse(df_plot$variable == "irredvarCov", "Cov1,2","RVar")))
p1 <- ggplot(df_plot,
       aes(x = alpha, y = terms, color = method)) + 
  geom_point(aes(shape = method), size = 2.5) + 
  geom_errorbar(aes(ymin = terms - se, ymax = terms + se, color = method, lty = method)) +
  theme_minimal(base_size = 16) + 
  #scale_color_manual(values = c("red", "black"))+
  facet_wrap(~variable, scales = "free", ncol = 4) + 
  theme(legend.position = "bottom", legend.title=element_blank()) +
  xlab(expression(alpha)) + 
  ylab("Variance")
ggsave(paste(figsdir, "/FINALplot_04_iredvar.pdf", sep = ""), plot = p1, device = "pdf", dpi = "retina", width = 12, height = 4.5)



