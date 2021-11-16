# coupled-bootstrap-risk-estimation

This repo contains the code used to produce the simulations in the paper "Unbiased Risk Estimation in the Normal Means Problem via Coupled Bootstrap Techniques", by Natalia L. Oliveira, Jing Lei, and Ryan J. Tibshirani.

## Abstract

We study a new method for estimating the risk of an arbitrary estimator of the mean vector in the classical normal means problem. The key idea is to generate two auxiliary data vectors, by adding two carefully constructed normal noise vectors to the original data vector; we then train the estimator of interest is on the first of these auxiliary data vectors and test it on the second. In order to stabilize the estimate of risk, we then average this procedure over multiple draws of the synthetic noise. A key aspect of this coupled bootstrap approach is that it delivers an unbiased estimate of risk under no assumptions on the estimator of the mean vector, albeit for a slightly “harder” version of the original normal means problem, where the error variance is inflated. We show that, under the assumptions required for Stein’s unbiased risk estimator (SURE), a limiting version of the coupled bootstrap estimator recovers SURE exactly (with an infinitesimal auxiliary noise variance and infinite bootstrap samples). We also analyze a bias-variance decomposition of the error of our risk estimator, to elucidate the effects of the variance of the auxiliary noise and the number of bootstrap samples on the accuracy of the risk estimator. Lastly, we demonstrate that our coupled bootstrap risk estimator performs quite favorably in simulated experiments.

## Code

### helper functions
- *estimators.R*: risk estimators
- *gs_for_fitting.R*: predictive modeling algorithms 
- *true_gs.R*: true mean functions
- *helpers_and_generators.R*: data generation, real risk approximation
- *helpers_dfsimulation.R.R*: functions that are used in the degrees of freedom simulation
- *helpers_image_denoising.R.R*: functions that are used in the image denoising application

### scripts
- *methods_comparison_comprehensive.R*: comprehensive set of simulations for comparing CB and BY estimators. Some of the result from this simulation are presented in Section 5.1. 
- *script_df_path.R*: simulation from Section 5.2
- *img_denoising_fused_lasso.r*: simulation from Section 5.3
- *biasbound.R*: simulation from Section F.1
- *RVar.R*: simulation from Section F.2
- *IVar.R*: simulation from Section F.3

### results and figures
- *savedfiles/* has the results in RDS format from the simulations
- *figures/* has the figures (in pdf) obtained from the simulations
