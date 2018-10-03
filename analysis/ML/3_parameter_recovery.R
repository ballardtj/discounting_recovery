library(tidyverse)
library(rstan)

sim_dat = expand.grid(m_a = seq(120,500,by=20),
            d_a = seq(2,40,by=2),
            m_b = 100,
            d_b = 0)

#get parameters from Bayesian models
models = list(hyperbolic = list(
  name = 'Hyperbolic model', 
  parms=c('k','sigma')),
  
  exponential = list(
    name = 'Exponential model',
    parms=c('k','sigma')),
  
  hyperbolic_gm = list( 
    name = 'Hyperbolic model (Green & Myerson version)',
    parms=c('k','s','sigma')),
  
  prop_diff = list( 
    name = 'Proportional difference model',
    parms=c('delta','sigma')),  
  
  tradeoff = list(
    name = 'Tradeoff model',
    parms=c("gamma","tau","theta",
             "kappa","alpha","eps")),
  
  ITCH = list(
    name = 'ITCH model',
    parms=c("B1","BxA","BxR",
             "BtA","BtR")),
  
  const_sens = list(
    name = 'Constant Sensitivity Model (Ebert & Prelec, 2007)',
    parms=c('alpha','beta',"sigma")),
  
  mazur1987 = list(
    name = 'Hyperboloid Model (Mazur, 1987)',
    parms=c('k','s','sigma')),
  
  loewenstein1992 = list(
    name = "Loewenstein & Prelec's (1992) Model",
    parms=c('alpha','beta','sigma')),
  
  mcclure2007 = list(
    name = "Double Exponential Model (McClure et al., 2007)",
    parms=c('omega','beta',"delta","sigma")),  
  
  killeen2009 = list(
    name = "Additive Utility Model (Killeen, 2009)",
    parms=c('beta',"lambda","sigma"))
)

#generate random parameter values from the hyper distributions

i = 1
seed = 12345
nsamples = 100  #number of different parameter values tested
nreps = 100     #number of times recovery is run for each parameter value
iter_samples = sample(1:8000,size=nsamples,replace=TRUE)
subject_samples = sample(1:40,size=nsamples,replace=TRUE)

#load likelihood
source(paste0("models/",names(models)[i],".R"))

#load fit object
load(file=paste0("data/derived/fit_Bayes_", names(models)[i],".RData"))
  
#extract parameters
posts=rstan::extract(fit)
  
#parameters
parm_names = models[names(models)[i]][[1]]$parms

#preallocate matrix of samples
parms = matrix(NA,nrow=nsamples,ncol=length(parm_names))
colnames(parms) = parm_names

for(j in 1:nsample){
  #store sampled parameter values 
  for(k in 1:length(parm_names)){
  #for each sample re-centre parameters (will need to be updated for different models)
  parms[j,k] = posts[paste0(parm_names[k],'_raw')][[1]][iter_samples[j],subject_samples[j]]*
      posts[paste0(parm_names[k],'_sd')][[1]][iter_samples[j]] +
      posts[paste0(parm_names[k],'_mean')][[1]][iter_samples[j]]
  }
  
  #get choice probabilities under sampled parameter values
  sim_dat$prob_a = likelihood(parms[j,],sim_dat)
  
  #for each repetition, generate random choices from probabilities, and fit model to the randomly genereated choices
  for(l in 1:nreps){
    #randomly generate choices
    sim_dat$choose_a = as.numeric( runif(nrow(sim_dat)) < sim_dat$prob_a )
    
    #fit model to randomly generated choices
    
  }
  
  
}


 for(k in 1:nsamples){
 
   
   
   
 }  
}


#for each sample generate data for those parameter values
source(paste0("models/",names(models)[i],".R"))



names(models)[i]

run simulate data for each parameter
models[names(models)[i]][[1]]


  #set up matrix of recentered parameters
  if(names(models)[i] %in% c('hyperbolic','exponential')){
    par_mat = cbind(
      posts$k_raw[samples[i],]*posts$k_sd + posts$k_mean,
      posts$sigma_raw[samples[i],]*posts$sigma_sd[samples[i]] + posts$sigma_mean[samples[i]]
    )
  }
  if(names(models)[i] %in% c('hyperbolic_gm','mazur1987')){
    par_mat = cbind(
      posts$k_raw[samples[i],]*posts$k_sd[samples[i]] + posts$k_mean[samples[i]],
      posts$s_raw[samples[i],]*posts$s_sd[samples[i]] + posts$s_mean[samples[i]],
      posts$sigma_raw[samples[i],]*posts$sigma_sd[samples[i]] + posts$sigma_mean[samples[i]]
    )
  }
  if(names(models)[i] %in% c('prop_diff')){
    par_mat = cbind(
      posts$delta_raw[samples[i],]*posts$delta_sd[samples[i]] + posts$delta_mean[samples[i]],
      posts$sigma_raw[samples[i],]*posts$sigma_sd[samples[i]] + posts$sigma_mean[samples[i]]
    )
  }
  if(names(models)[i] %in% c('tradeoff')){
    par_mat = cbind(
      posts$gamma_raw[samples[i],]*posts$gamma_sd[samples[i]] + posts$gamma_mean[samples[i]],
      posts$tau_raw[samples[i],]*posts$tau_sd[samples[i]] + posts$tau_mean[samples[i]],
      posts$theta_raw[samples[i],]*posts$theta_sd[samples[i]] + posts$theta_mean[samples[i]],
      posts$kappa_raw[samples[i],]*posts$kappa_sd[samples[i]] + posts$kappa_mean[samples[i]],
      posts$alpha_raw[samples[i],]*posts$alpha_sd[samples[i]] + posts$alpha_mean[samples[i]],
      posts$eps_raw[samples[i],]*posts$eps_sd[samples[i]] + posts$eps_mean[samples[i]]
    )
  }
  if(names(models)[i] %in% c('ITCH')){
    par_mat = cbind(
      posts$B1_raw[samples[i],]*posts$B1_sd[samples[i]] + posts$B1_mean[samples[i]],
      posts$BxA_raw[samples[i],]*posts$BxA_sd[samples[i]] + posts$BxA_mean[samples[i]],
      posts$BxR_raw[samples[i],]*posts$BxR_sd[samples[i]] + posts$BxR_mean[samples[i]],
      posts$BtA_raw[samples[i],]*posts$BtA_sd[samples[i]] + posts$BtA_mean[samples[i]],
      posts$BtR_raw[samples[i],]*posts$BtR_sd[samples[i]] + posts$BtR_mean[samples[i]]
    )
  }
  if(names(models)[i] %in% c('const_sens','loewenstein1992')){
    par_mat = cbind(
      posts$alpha_raw[samples[i],]*posts$alpha_sd[samples[i]] + posts$alpha_mean[samples[i]],
      posts$beta_raw[samples[i],]*posts$beta_sd[samples[i]] + posts$beta_mean[samples[i]],
      posts$sigma_raw[samples[i],]*posts$sigma_sd[samples[i]] + posts$sigma_mean[samples[i]]
    )
  }
  if(names(models)[i] %in% c('mcclure2007')){
    par_mat = cbind(
      posts$omega[samples[i],],
      posts$beta_raw[samples[i],]*posts$beta_sd[samples[i]] + posts$beta_mean[samples[i]],
      posts$delta_raw[samples[i],]*posts$delta_sd[samples[i]] + posts$delta_mean[samples[i]],
      posts$sigma_raw[samples[i],]*posts$sigma_sd[samples[i]] + posts$sigma_mean[samples[i]]
    )
  }
  if(names(models)[i] %in% c('killeen2009')){
    par_mat = cbind(
      posts$beta_raw[samples[i],]*posts$beta_sd[samples[i]] + posts$beta_mean[samples[i]],
      posts$lambda_raw[samples[i],]*posts$lambda_sd[samples[i]] + posts$lambda_mean[samples[i]],
      posts$sigma_raw[samples[i],]*posts$sigma_sd[samples[i]] + posts$sigma_mean[samples[i]]
    )
  }
  
  
  
  
 
  
  par_inds_start = grep("choose_a",names(dat_tmp)) + 1
  par_inds_end = grep("p_a",names(dat_tmp)) - 1
  print(
    ggpairs(dat_tmp[dat_tmp$trial==1,par_inds_start:par_inds_end])
  )
  cat('\n') 
}

#load data
dat=readRDS(file="data/raw/data_delay_lba2018.rds")

#format data for model fitting
dat$m_a = dat$amount
dat$d_a = dat$delay
dat$m_b = 100
dat$d_b = 0    #note: b is smaller, sooner option
dat$choose_a = dat$Choice
Nsubj = length(unique(dat$subject))






#The conclusion is that all these models fit well. But theres a question
#as to how well their parameters can be interpretted. 
rm(list=ls())

#source libraries
library(tidyverse)

plot_pp = function(dat_tmp){
  dat_pl = dat_tmp %>%
    group_by(amount,delay) %>%
    summarise(obs = mean(choose_a),
              obs_se = sd(choose_a)/sqrt(length(choose_a)),
              pred = mean(p_a))
  
  ggplot(dat_pl,aes(x=delay)) +
    geom_point(aes(y=obs),colour="red") +
    geom_errorbar(aes(ymin=obs-obs_se,ymax=obs+obs_se),colour="red") +
    geom_line(aes(y=pred),colour="blue") +
    facet_wrap(~amount)
}

get_bic = function(dat_tmp){
  #make sure there are no cases where p_a = 0 or 1
  p_a = pmin(pmax(dat_tmp$p_a,0.001),0.999) 
  #calculate summed log likelihood
  lnl = sum(log(dat_tmp$choose_a*p_a + (1-dat_tmp$choose_a)*(1-p_a)))
  #compute number of parameters (based on how many columns are between observed choice and predicted choice)
  nparms = grep("p_a",names(dat_tmp))-grep("choose_a",names(dat_tmp)) - 1
  #compute bic
  bic = log(length(p_a))*nparms - 2*lnl 
  return(list(lnl=lnl,nparms=nparms,bic=bic))
}

#----------------------------------------------------
#hyperbolic
load(file="data/derived/fit_ML_hyperbolic.RData")
plot_pp(dat_tmp)
get_bic(dat_tmp)


#----------------------------------------------------
#exponential
load(file="data/derived/fit_ML_exponential.RData")
plot_pp(dat_tmp)
get_bic(dat_tmp)

#----------------------------------------------------
#green & myerson hyperbolic
load(file="data/derived/fit_ML_hyperbolic_gm.RData")
plot_pp(dat_tmp)
get_bic(dat_tmp)

#----------------------------------------------------
#fit proportional difference model
load(file="data/derived/fit_ML_prop_diff.RData")
plot_pp(dat_tmp)
get_bic(dat_tmp)

#----------------------------------------------------
#fit tradeoff model
load(file="data/derived/fit_ML_tradeoff.RData")
plot_pp(dat_tmp)
get_bic(dat_tmp)

#----------------------------------------------------
#fit ITCH model
load(file="data/derived/fit_ML_ITCH.RData")
plot_pp(dat_tmp)
get_bic(dat_tmp)

#----------------------------------------------------
#fit constant sensitivity model
load(file="data/derived/fit_ML_const_sens.RData")
plot_pp(dat_tmp)
get_bic(dat_tmp)

#----------------------------------------------------
#fit hyperbaloid model (mazur 1987)
load(file="data/derived/fit_ML_mazur1987.RData")
plot_pp(dat_tmp)
get_bic(dat_tmp)

#----------------------------------------------------
#fit Loewenstein 1992 model
load(file="data/derived/fit_ML_leowenstein1992.RData")
plot_pp(dat_tmp)
get_bic(dat_tmp)


#----------------------------------------------------
#fit McClure 2007 double exponential model
load(file="data/derived/fit_ML_mcclure2007.RData")
plot_pp(dat_tmp)
get_bic(dat_tmp)


#----------------------------------------------------
#fit Killeen 2009 model
load(file="data/derived/fit_ML_killeen2009.RData")
plot_pp(dat_tmp)
get_bic(dat_tmp)

# 
# 
# parms = dat %>% group_by(subject) %>%
#   summarise(k = mean(k),
#             sigma = mean(sigma))
# 
# ggplot(parms,aes(x=k,y=sigma)) +
#   geom_density_2d()  +
#   geom_point()

