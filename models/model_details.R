

models = list(hyperbolic = list(
  name = 'Hyperbolic model', 
  parms=c('k','sigma'),
  label=c('Hyperbolic'),
  hypers=c('k_mean','k_sd','sigma_mean','sigma_sd')),
  
  exponential = list(
    name = 'Exponential model',
    parms=c('k','sigma'),
    label=c('Exponential'),
    hypers=c('k_mean','k_sd','sigma_mean','sigma_sd')),
  
  hyperbolic_gm = list( 
    name = 'Hyperbolic model (Green & Myerson version)',
    parms=c('k','s','sigma'),
    label=c('Hyperboloid'),
    hypers=c('k_mean','k_sd','s_mean','s_sd','sigma_mean','sigma_sd')),
  
  prop_diff = list( 
    name = 'Proportional difference model',
    parms=c('delta','sigma'),
    label=c('Proportional Difference'),
    hypers=c('delta_mean','delta_sd','sigma_mean','sigma_sd')),  

  tradeoff = list(
    name = 'Tradeoff model',
    parms=c("gamma","tau","theta",
            "kappa","alpha","eps"),
    label=c('Tradeoff'),
    hypers=c("gamma_mean","gamma_sd","tau_mean","tau_sd","theta_mean","theta_sd",
             "kappa_mean","kappa_sd","alpha_mean","alpha_sd","eps_mean","eps_sd")),
      
  ITCH = list(
    name = 'ITCH model',
    parms=c("B1","BxA","BxR",
            "BtA","BtR"),
    label=c('ITCH'),
    hypers=c("B1_mean","B1_sd","BxA_mean","BxA_sd","BxR_mean","BxR_sd",
             "BtA_mean","BtA_sd","BtR_mean","BtR_sd")),
  
  const_sens = list(
    name = 'Constant Sensitivity Model (Ebert & Prelec, 2007)',
    parms=c('alpha','beta',"sigma"),
    label=c('Constant Sensitivity'),
    hypers=c('alpha_mean','alpha_sd','beta_mean','beta_sd','sigma_mean','sigma_sd')),
  
  mazur1987 = list(
    name = 'Hyperboloid Model (Mazur, 1987)',
    parms=c('k','s','sigma'),
    label=c('Generalized Hyperbolic'),
    hypers=c('k_mean','k_sd','s_mean','s_sd','sigma_mean','sigma_sd')),
  
  loewenstein1992 = list(
    name = "Loewenstein & Prelec's (1992) Model",
    parms=c('alpha','beta','sigma'),
    label=c('Generalized Hyperbola'),
    hypers=c('alpha_mean','alpha_sd','beta_mean','beta_sd','sigma_mean','sigma_sd')),
  
  mcclure2007 = list(
    name = "Double Exponential Model (McClure et al., 2007)",
    parms=c('omega','beta',"delta","sigma"),  
    label=c('Double Exponential'),
    hypers=c('omega_a','omega_b','beta_mean','beta_sd',"delta_mean","delta_sd","sigma_mean",'sigma_sd')), 
  
  killeen2009 = list(
    name = "Additive Utility Model (Killeen, 2009)",
    parms=c('alpha','beta',"lambda","sigma"),
    label=c('Additive Utility'),
    hypers=c('alpha_a','alpha_b','beta_a','beta_b',"lambda_mean","lambda_sd","sigma_mean",'sigma_sd'))
)


