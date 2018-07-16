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
  lower=c(0,0)
  upper=c(100,1000)
  
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
  lower=c(0,0,0,0,0,0)
  upper=c(50,50,50,50,50,50)
  
  output=DEoptim(fn=tradeoff_wrapper,
                 dat=dat_tmp[dat_tmp$subject==i,],
                 lower=lower,
                 upper=upper,
                 control=list(trace=0))
  
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
  lower=c(-50,-50,-50,-50,-50,-50)
  upper=c(50,50,50,50,50,50)
  
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


# 
# 
# parms = dat %>% group_by(subject) %>%
#   summarise(k = mean(k),
#             sigma = mean(sigma))
# 
# ggplot(parms,aes(x=k,y=sigma)) +
#   geom_density_2d()  +
#   geom_point()

