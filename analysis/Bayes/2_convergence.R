
#   Model description:
#   
#   * Hierarchical
# * Chains: 8
# * Burnin: 1000
# * Post-burnin: 1000
# * Thinning: 1
# * Total samples: 8000


library(rstan)
library(coda)
library(knitr)

library(shinystan)

launch_shinystan(fit)

models = list(hyperbolic = list(
  name = 'Hyperbolic model', 
  hypers=c('k_mean','k_sd','sigma_mean','sigma_sd')),
  
  exponential = list(
    name = 'Exponential model',
    hypers=c('k_mean','k_sd','sigma_mean','sigma_sd')),
  
  hyperbolic_gm = list( 
    name = 'Hyperbolic model (Green & Myerson version)',
    hypers=c('k_mean','k_sd','s_mean','s_sd','sigma_mean','sigma_sd')),
  
  prop_diff = list( 
    name = 'Proportional difference model',
    hypers=c('delta_mean','delta_sd','sigma_mean','sigma_sd')),  
  
  tradeoff = list(
    name = 'Tradeoff model',
    hypers=c("gamma_mean","gamma_sd","tau_mean","tau_sd","theta_mean","theta_sd",
             "kappa_mean","kappa_sd","alpha_mean","alpha_sd","eps_mean","eps_sd")),
  
  ITCH = list(
    name = 'ITCH model',
    hypers=c("B1_mean","B1_sd","BxA_mean","BxA_sd","BxR_mean","BxR_sd",
             "BtA_mean","BtA_sd","BtR_mean","BtR_sd")),
  
  const_sens = list(
    name = 'Constant Sensitivity Model (Ebert & Prelec, 2007)',
    hypers=c('alpha_mean','alpha_sd','beta_mean','beta_sd','sigma_mean','sigma_sd')),
  
  mazur1987 = list(
    name = 'Hyperboloid Model (Mazur, 1987)',
    hypers=c('k_mean','k_sd','s_mean','s_sd','sigma_mean','sigma_sd')),
  
  loewenstein1992 = list(
    name = "Loewenstein & Prelec's (1992) Model",
    hypers=c('alpha_mean','alpha_sd','beta_mean','beta_sd','sigma_mean','sigma_sd')),
  
  mcclure2007 = list(
    name = "Double Exponential Model (McClure et al., 2007)",
    hypers=c('omega_a','omega_b','beta_mean','beta_sd',"delta_mean","delta_sd","sigma_mean",'sigma_sd')),  
  
  killeen2009 = list(
    name = "Additive Utility Model (Killeen, 2009)",
    hypers=c('beta_mean','beta_sd',"lambda_mean","lambda_sd","sigma_mean",'sigma_sd'))
)

stan2coda <- function(fit) {
  mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}

for (i in 1:length(models)){

  print(models[[i]]$name)
  
  load(file=paste0("data/derived/fits/fit_Bayes_", names(models)[i],".RData"))
  
  smrystan = summary(fit)
  print(kable(smrystan[[1]][models[[i]]$hypers,]),format="markdown")
  
  print(kable(effectiveSize(stan2coda(fit))[models[[i]]$hypers],format="markdown",col.names="n_eff"))
  
  print(kable(gelman.diag(stan2coda(fit))[[1]][models[[i]]$hypers,],format="markdown"))
  
  traceplot = rstan::traceplot(fit,par=models[[i]]$hypers)
  ggsave(file=paste0("figures/convergence_trace_",names(models)[i],".pdf"),plot=traceplot)
  
  
  pdf(file=paste0("figures/convergence_pairs_",names(models)[i],".pdf"),width=10,height=10)
    pairs(fit,pars=models[[i]]$hypers)
  dev.off()
  
}


