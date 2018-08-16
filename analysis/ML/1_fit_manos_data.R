rm(list=ls())

#source libraries
library(tidyverse)
library(DEoptim)

#load data
dat=readRDS(file="data/raw/data_delay_lba2018.rds")

#format data for model fitting
dat$m_a = dat$amount
dat$d_a = dat$delay
dat$m_b = 100
dat$d_b = 0    #note: b is smaller, sooner option
dat$choose_a = dat$Choice
Nsubj = length(unique(dat$subject))

#wrapper function 
wrapper=function(pars,dat){
  p_a = likelihood(pars,dat)
  p_a = 0.001 + p_a*0.998  
  neglnLs = -log(dat$choose_a*p_a + (1-dat$choose_a)*(1-p_a)) 
  return(sum(neglnLs))
}

#par list
model_pars = list(
  hyperbolic = list(pars=c('k','sigma'),
                    upper=c(10,100),
                    lower=c(0,0)),
  exponential = list(pars=c('k','sigma'),
                     upper=c(10,100),
                     lower=c(0,0)),
  hyperbolic_gm = list(pars=c('k','s','sigma'),  
                       upper=c(10,100,100),
                       lower=c(0,0,0)),
  prop_diff = list(pars=c('delta','sigma'),
                   upper=c(1,100),
                   lower=c(-1,0)),
  tradeoff = list(pars=c('gamma','tau','theta','kappa','alpha','eps'),
                  upper=c(1,1,10,1,1,1),
                  lower=c(1e-10,1e-10,1,1e-10,1e-10,1e-2)),
  ITCH = list(pars=c("B1","BxA","BxR","BtA","BtR"),
              upper=c(50,50,50,0,0), 
              lower=c(-50,0,0,-50,-50)),
  const_sens = list(pars=c("alpha","beta","sigma"),
                    upper=c(100,100,100),lower=c(0,0,0)),
  mazur1987 = list(pars=c("k","s","sigma"),upper=c(100,100,100),
                   lower=c(0,0,0)),
  loewenstein1992 = list(pars=c("alpha","beta","sigma"),  
                         upper=c(1000,1000,1000), 
                         lower=c(0,0,0)),
  mcclure2007 = list(pars=c("omega","beta","delta","sigma"),
                     upper=c(1,10,10,100),
                     lower=c(0,0,0,0)),
  killeen2009 = list(pars=c('beta','lambda','sigma'),
                     upper=c(10,1000,1000),
                     lower=c(0,0,0))
)

for(i in 1:length(model_pars)){
  
  #source model
  model_name = names(model_pars)[i]
  source(paste0("models/",model_name,".R"))
  
  #add empty columns in data frame for free parameters
  dat_tmp = dat
  for(j in 1:length(model_pars[[i]]$pars)){
    dat_tmp[,model_pars[[i]]$pars[j]] = NA
  }
  dat_tmp$p_a = NA
  
  for(k in 1:Nsubj){
    output=DEoptim(fn=wrapper,
                  dat=dat_tmp[dat_tmp$subject==k,],
                  lower=model_pars[[i]]$lower,
                  upper=model_pars[[i]]$upper,
                  control=list(trace=0))
  
    output2 = optim(par=output$optim$bestmem,
                    fn=wrapper,
                    dat=dat_tmp[dat_tmp$subject==k,],
                    method="L-BFGS-B",
                    lower=model_pars[[i]]$lower,
                    upper=model_pars[[i]]$upper)
  
    for(j in 1:length(model_pars[[i]]$pars)){
      dat_tmp[dat_tmp$subject==k,model_pars[[i]]$pars[j]] = output2$par[j]
    }
  
  dat_tmp$p_a[dat_tmp$subject==k] = likelihood(output2$par,dat_tmp[dat_tmp$subject==k,])
  }
  
  save(dat_tmp,file=paste0("data/derived/fit_ML_", model_name,".RData"))
}
                  



