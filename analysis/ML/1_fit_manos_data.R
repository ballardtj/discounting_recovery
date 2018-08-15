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
model_pars = list(hyperbolic = c('k','sigma'),
                  exponential = c('k','sigma'),
                  hyperbolic_gm = c('k','s','sigma'),
                  prop_diff = c('delta','sigma'),
                  tradeoff = c('gamma','tau','theta','kappa','alpha','eps'),
                  ITCH = c("B1","BxA","BxR","BtA","BtR"),
                  const_sens = c("alpha","beta","sigma"),
                  mazur1987 = c("k","s","sigma"),
                  leowenstein1992 = c("alpha","beta","sigma"),
                  mcclure2007 = c("omega","beta","delta","sigma"),
                  killeen2009 = c('beta','lambda','sigma'))

for(i in 1:length(model_pars)){
  
  #source model
  model_name = names(model_pars)[i]
  source(paste0("models/",model_name,".R"))
  
  #add empty columns in data frame for free parameters
  dat_tmp = dat
  for(j in 1:length(model_pars[[i]])){
    dat_tmp[,model_pars[[i]][j]] = NA
  }
  
}
                  


#----------------------------------------------------
#fit hyperbolic
source("models/hyperbolic.R")
dat_tmp=dat
dat_tmp$k = NA
dat_tmp$sigma = NA

for(i in 1:Nsubj){
  lower=c(0,0)
  upper=c(10,100)
  
  output=DEoptim(fn=hyperbolic_wrapper,
                 dat=dat_tmp[dat_tmp$subject==i,],
                 lower=lower,
                 upper=upper,
                 control=list(trace=0))
  
  output2 = optim(par=output$optim$bestmem,
        fn=hyperbolic_wrapper,
        dat=dat_tmp[dat_tmp$subject==i,],
        method="L-BFGS-B",
        lower=lower,
        upper=upper)

  dat_tmp$k[dat_tmp$subject==i] = output2$par[1]
  dat_tmp$sigma[dat_tmp$subject==i] = output2$par[2]
  dat_tmp$p_a[dat_tmp$subject==i] = hyperbolic_likelihood(dat_tmp$m_a[dat_tmp$subject==i],dat_tmp$d_a[dat_tmp$subject==i],dat_tmp$m_b[dat_tmp$subject==i],dat_tmp$d_b[dat_tmp$subject==i], 
                                                  output2$par[1],output2$par[2])
}

save(dat_tmp,file="data/derived/fit_ML_hyperbolic.RData")


#----------------------------------------------------
#fit exponential
source("models/exponential.R")
dat_tmp=dat
dat_tmp$k = NA
dat_tmp$sigma = NA

for(i in 1:Nsubj){
  lower=c(0,0)
  upper=c(10,100)
  
  output=DEoptim(fn=exponential_wrapper,
                 dat=dat_tmp[dat_tmp$subject==i,],
                 lower=lower,
                 upper=upper,
                 control=list(trace=0))
  
  output2 = optim(par=output$optim$bestmem,
                  fn=exponential_wrapper,
                  dat=dat_tmp[dat_tmp$subject==i,],
                  method="L-BFGS-B",
                  lower=lower,
                  upper=upper)
  
  
  dat_tmp$k[dat_tmp$subject==i] = output2$par[1]
  dat_tmp$sigma[dat_tmp$subject==i] = output2$par[2]
  dat_tmp$p_a[dat_tmp$subject==i] = exponential_likelihood(dat_tmp$m_a[dat_tmp$subject==i],dat_tmp$d_a[dat_tmp$subject==i],dat_tmp$m_b[dat_tmp$subject==i],dat_tmp$d_b[dat_tmp$subject==i], 
                                                  output2$par[1],output2$par[2])
}

save(dat_tmp,file="data/derived/fit_ML_exponential.RData")

#----------------------------------------------------
#fit nyerson & hyperbolic
source("models/hyperbolic_gm.R")
dat_tmp=dat
dat_tmp$k = NA
dat_tmp$s = NA
dat_tmp$sigma = NA

for(i in 1:Nsubj){
  lower=c(0,0,0)
  upper=c(10,100,100)
  
  output=DEoptim(fn=hyperbolic_gm_wrapper,
                 dat=dat_tmp[dat_tmp$subject==i,],
                 lower=lower,
                 upper=upper,
                 control=list(trace=0))
  
  output2 = optim(par=output$optim$bestmem,
                  fn=hyperbolic_gm_wrapper,
                  dat=dat_tmp[dat_tmp$subject==i,],
                  method="L-BFGS-B",
                  lower=lower,
                  upper=upper)
  
  
  dat_tmp$k[dat_tmp$subject==i] = output2$par[1]
  dat_tmp$s[dat_tmp$subject==i] = output2$par[2]
  dat_tmp$sigma[dat_tmp$subject==i] = output2$par[3]
  dat_tmp$p_a[dat_tmp$subject==i] = hyperbolic_gm_likelihood(dat_tmp$m_a[dat_tmp$subject==i],dat_tmp$d_a[dat_tmp$subject==i],dat_tmp$m_b[dat_tmp$subject==i],dat_tmp$d_b[dat_tmp$subject==i], 
                                                  output2$par[1],output2$par[2],output2$par[3])
}

save(dat_tmp,file="data/derived/fit_ML_hyperbolic_gm.RData")

#----------------------------------------------------
#fit proportional difference model
source("models/prop_diff.R")
dat_tmp=dat
dat_tmp$delta = NA
dat_tmp$sigma = NA

for(i in 1:Nsubj){
  lower=c(-1,0)
  upper=c(1,100)
  
  output=DEoptim(fn=prop_diff_wrapper,
                 dat=dat_tmp[dat_tmp$subject==i,],
                 lower=lower,
                 upper=upper,
                 control=list(trace=0))
  
  output2 = optim(par=output$optim$bestmem,
                  fn=prop_diff_wrapper,
                  dat=dat_tmp[dat_tmp$subject==i,],
                  method="L-BFGS-B",
                  lower=lower,
                  upper=upper)
  
  dat_tmp$delta[dat_tmp$subject==i] = output2$par[1]
  dat_tmp$sigma[dat_tmp$subject==i] = output2$par[2]
  dat_tmp$p_a[dat_tmp$subject==i] = prop_diff_likelihood(dat_tmp$m_a[dat_tmp$subject==i],dat_tmp$d_a[dat_tmp$subject==i],dat_tmp$m_b[dat_tmp$subject==i],dat_tmp$d_b[dat_tmp$subject==i], 
                                                     output2$par[1],output2$par[2])
  
}

save(dat_tmp,file="data/derived/fit_ML_prop_diff.RData")

#----------------------------------------------------
#fit tradeoff model
source("models/tradeoff.R")
dat_tmp=dat
dat_tmp$gamma = NA
dat_tmp$tau = NA
dat_tmp$theta = NA
dat_tmp$kappa = NA
dat_tmp$alpha = NA
dat_tmp$eps = NA

for(i in 1:Nsubj){
  lower=c(1e-10,1e-10,1,1e-10,1e-10,1e-2)
  upper=c(1,1,10,1,1,1)
  
  output=DEoptim(fn=tradeoff_wrapper,
                 dat=dat_tmp[dat_tmp$subject==i,],
                 lower=lower,
                 upper=upper,
                 control=list(trace=F))
  
  output2 = optim(par=output$optim$bestmem,
                  fn=tradeoff_wrapper,
                  dat=dat_tmp[dat_tmp$subject==i,],
                  method="L-BFGS-B",
                  lower=lower,
                  upper=upper)
  
  dat_tmp$gamma[dat_tmp$subject==i] = output2$par[1]
  dat_tmp$tau[dat_tmp$subject==i] = output2$par[2]
  dat_tmp$theta[dat_tmp$subject==i] = output2$par[3]
  dat_tmp$kappa[dat_tmp$subject==i] = output2$par[4]
  dat_tmp$alpha[dat_tmp$subject==i] = output2$par[5]
  dat_tmp$eps[dat_tmp$subject==i] = output2$par[6]

  dat_tmp$p_a[dat_tmp$subject==i] = tradeoff_likelihood(dat_tmp$m_a[dat_tmp$subject==i],dat_tmp$d_a[dat_tmp$subject==i],dat_tmp$m_b[dat_tmp$subject==i],dat_tmp$d_b[dat_tmp$subject==i], 
                                                 output2$par[1],output2$par[2],output2$par[3],
                                                 output2$par[4],output2$par[5],output2$par[6])
  
}

save(dat_tmp,file="data/derived/fit_ML_tradeoff.RData")

#----------------------------------------------------
#fit ITCH model
source("models/ITCH.R")
dat_tmp=dat
dat_tmp$B1 = NA
dat_tmp$BxA = NA
dat_tmp$BxR = NA
dat_tmp$BtA = NA
dat_tmp$BtR = NA

for(i in 1:Nsubj){
  lower=c(-50,0,0,-50,-50)
  upper=c(50,50,50,0,0)
  
  output=DEoptim(fn=ITCH_wrapper,
                 dat=dat_tmp[dat_tmp$subject==i,],
                 lower=lower,
                 upper=upper,
                 control=list(trace=0))
  
  output2 = optim(par=output$optim$bestmem,
                  fn=ITCH_wrapper,
                  dat=dat_tmp[dat_tmp$subject==i,],
                  method="L-BFGS-B",
                  lower=lower,
                  upper=upper)
  
  
  dat_tmp$B1[dat_tmp$subject==i] = output2$par[1]
  dat_tmp$BxA[dat_tmp$subject==i] = output2$par[2]
  dat_tmp$BxR[dat_tmp$subject==i] = output2$par[3]
  dat_tmp$BtA[dat_tmp$subject==i] = output2$par[4]
  dat_tmp$BtR[dat_tmp$subject==i] = output2$par[5]

  dat_tmp$p_a[dat_tmp$subject==i] = ITCH_likelihood(dat_tmp$m_a[dat_tmp$subject==i],dat_tmp$d_a[dat_tmp$subject==i],dat_tmp$m_b[dat_tmp$subject==i],dat_tmp$d_b[dat_tmp$subject==i],  
                                                output2$par[1],output2$par[2],output2$par[3],
                                                output2$par[4],output2$par[5])
  
}

save(dat_tmp,file="data/derived/fit_ML_ITCH.RData")


#----------------------------------------------------
#fit constant sensitivity model

source("models/const_sens.R")
dat_tmp=dat
dat_tmp$alpha = NA
dat_tmp$beta = NA
dat_tmp$sigma = NA

for(i in 1:Nsubj){
  lower=c(0,0,0)
  upper=c(100,100,100)
  
  output=DEoptim(fn=const_sens_wrapper,
                 dat=dat_tmp[dat_tmp$subject==i,],
                 lower=lower,
                 upper=upper,
                 control=list(trace=0))
  
  output2 = optim(par=output$optim$bestmem,
                  fn=const_sens_wrapper,
                  dat=dat_tmp[dat_tmp$subject==i,],
                  method="L-BFGS-B",
                  lower=lower,
                  upper=upper)
  
  
  dat_tmp$alpha[dat_tmp$subject==i] = output2$par[1]
  dat_tmp$beta[dat_tmp$subject==i] = output2$par[2]
  dat_tmp$sigma[dat_tmp$subject==i] = output2$par[3]

  
  dat_tmp$p_a[dat_tmp$subject==i] = const_sens_likelihood(dat_tmp$m_a[dat_tmp$subject==i],dat_tmp$d_a[dat_tmp$subject==i],dat_tmp$m_b[dat_tmp$subject==i],dat_tmp$d_b[dat_tmp$subject==i],  
                                                    output2$par[1],output2$par[2],output2$par[3])
  
}

save(dat_tmp,file="data/derived/fit_ML_const_sens.RData")


#----------------------------------------------------
#fit hyperbaloid model (mazur 1987)

source("models/mazur1987.R")
dat_tmp=dat
dat_tmp$k = NA
dat_tmp$s = NA
dat_tmp$sigma = NA

for(i in 1:Nsubj){
  lower=c(0,0,0)
  upper=c(10,10,100)
  
  output=DEoptim(fn=mazur_wrapper,
                 dat=dat_tmp[dat_tmp$subject==i,],
                 lower=lower,
                 upper=upper,
                 control=list(trace=0))
  
  output2 = optim(par=output$optim$bestmem,
                  fn=mazur_wrapper,
                  dat=dat_tmp[dat_tmp$subject==i,],
                  method="L-BFGS-B",
                  lower=lower,
                  upper=upper)
  
  dat_tmp$k[dat_tmp$subject==i] = output2$par[1]
  dat_tmp$s[dat_tmp$subject==i] = output2$par[2]
  dat_tmp$sigma[dat_tmp$subject==i] = output2$par[3]
  
  
  dat_tmp$p_a[dat_tmp$subject==i] = mazur_likelihood(dat_tmp$m_a[dat_tmp$subject==i],dat_tmp$d_a[dat_tmp$subject==i],dat_tmp$m_b[dat_tmp$subject==i],dat_tmp$d_b[dat_tmp$subject==i],  
                                                          output2$par[1],output2$par[2],output2$par[3])
  
}

save(dat_tmp,file="data/derived/fit_ML_mazur1987.RData")



#----------------------------------------------------
#fit Loewenstein 1992 model

source("models/loewenstein1992.R")
dat_tmp=dat
dat_tmp$alpha = NA
dat_tmp$beta = NA
dat_tmp$sigma = NA

for(i in 1:Nsubj){
  lower=c(0,0,0)
  upper=c(1000,1000,1000)
  
  output=DEoptim(fn=leow_wrapper,
                 dat=dat_tmp[dat_tmp$subject==i,],
                 lower=lower,
                 upper=upper,
                 control=list(trace=0))
  
  output2 = optim(par=output$optim$bestmem,
                  fn=leow_wrapper,
                  dat=dat_tmp[dat_tmp$subject==i,],
                  method="L-BFGS-B",
                  lower=lower,
                  upper=upper)
  
  
  dat_tmp$alpha[dat_tmp$subject==i] = output2$par[1]
  dat_tmp$beta[dat_tmp$subject==i] = output2$par[2]
  dat_tmp$sigma[dat_tmp$subject==i] = output2$par[3]
  
  
  dat_tmp$p_a[dat_tmp$subject==i] = leow_likelihood(dat_tmp$m_a[dat_tmp$subject==i],dat_tmp$d_a[dat_tmp$subject==i],dat_tmp$m_b[dat_tmp$subject==i],dat_tmp$d_b[dat_tmp$subject==i],  
                                                          output2$par[1],output2$par[2],output2$par[3])
  
}

save(dat_tmp,file="data/derived/fit_ML_leowenstein1992.RData")



#----------------------------------------------------
#fit McClure 2007 double exponential model

source("models/mcclure2007.R")
dat_tmp=dat
dat_tmp$omega = NA
dat_tmp$delta = NA
dat_tmp$beta = NA
dat_tmp$sigma = NA

for(i in 1:Nsubj){
  lower=c(0,0,0,0)
  upper=c(1,10,10,100)
  
  output=DEoptim(fn=double_exp_wrapper,
                 dat=dat_tmp[dat_tmp$subject==i,],
                 lower=lower,
                 upper=upper,
                 control=list(trace=0))
  
  output2 = optim(par=output$optim$bestmem,
                  fn=double_exp_wrapper,
                  dat=dat_tmp[dat_tmp$subject==i,],
                  method="L-BFGS-B",
                  lower=lower,
                  upper=upper)
  
  dat_tmp$omega[dat_tmp$subject==i] = output2$par[1]
  dat_tmp$delta[dat_tmp$subject==i] = output2$par[2]
  dat_tmp$beta[dat_tmp$subject==i] = output2$par[3]
  dat_tmp$sigma[dat_tmp$subject==i] = output2$par[4]
  
  
  dat_tmp$p_a[dat_tmp$subject==i] = double_exp_likelihood(dat_tmp$m_a[dat_tmp$subject==i],dat_tmp$d_a[dat_tmp$subject==i],dat_tmp$m_b[dat_tmp$subject==i],dat_tmp$d_b[dat_tmp$subject==i],  
                                                    output2$par[1],output2$par[2],output2$par[3],output2$par[4])
  
}

save(dat_tmp,file="data/derived/fit_ML_mcclure2007.RData")


#----------------------------------------------------
#fit Killeen 2009 model

source("models/killeen2009.R")
dat_tmp=dat
dat_tmp$beta = NA
dat_tmp$lambda = NA
dat_tmp$sigma = NA

for(i in 1:Nsubj){
  lower=c(0,0,0)
  upper=c(10,1000,1000)
  
  output=DEoptim(fn=killeen_wrapper,
                 dat=dat_tmp[dat_tmp$subject==i,],
                 lower=lower,
                 upper=upper,
                 control=list(trace=0))
  
  output2 = optim(par=output$optim$bestmem,
                  fn=killeen_wrapper,
                  dat=dat_tmp[dat_tmp$subject==i,],
                  method="L-BFGS-B",
                  lower=lower,
                  upper=upper)
  
  dat_tmp$beta[dat_tmp$subject==i] = output2$par[1]
  dat_tmp$lambda[dat_tmp$subject==i] = output2$par[2]
  dat_tmp$sigma[dat_tmp$subject==i] = output2$par[3]
  
  
  dat_tmp$p_a[dat_tmp$subject==i] = killeen_likelihood(dat_tmp$m_a[dat_tmp$subject==i],dat_tmp$d_a[dat_tmp$subject==i],dat_tmp$m_b[dat_tmp$subject==i],dat_tmp$d_b[dat_tmp$subject==i],  
                                                    output2$par[1],output2$par[2],output2$par[3])
  
}

save(dat_tmp,file="data/derived/fit_ML_killeen2009.RData")



# 
# 
# parms = dat %>% group_by(subject) %>%
#   summarise(k = mean(k),
#             sigma = mean(sigma))
# 
# ggplot(parms,aes(x=k,y=sigma)) +
#   geom_density_2d()  +
#   geom_point()

