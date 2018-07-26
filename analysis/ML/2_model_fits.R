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

