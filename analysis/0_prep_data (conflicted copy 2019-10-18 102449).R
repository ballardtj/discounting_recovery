rm(list=ls())

#load packages
library(plyr)
library(tidyverse)
library(R.matlab)
library(rstan)
library(shinystan)

####################
for(exp in 1:2){
  
  if(exp == 1){
    #load data
    dat=readRDS(file="data/raw/delay_data_lba_2019.rds") %>%
      mutate(subject = as.numeric(as.factor(subject)))
    
    #format data for model fitting  
    dat$m_a = dat$amount
    dat$d_a = dat$delay
    dat$m_b = 100
    dat$d_b = 0    #note: b is smaller, sooner option
    dat$choose_a = dat$Choice
  }
  
  if(exp == 2){
    #load data and transform into R data frame (a = LL, b = SS)
    tmp = readMat("data/raw/datS3.mat")
    tmp = tmp[[1]]
    dat = adply(tmp,c(2,1))
    names(dat) = c("trial","subject","m_b","d_b","m_a","d_a","choose_a")
  }
  
  
  Nsubj = length(unique(dat$subject))
  
  #Get number of observations per subject
  Obs_each_subj = dat %>%
    group_by(subject) %>%
    dplyr::summarise(count = n()) %>% .$count
  
  #Get max number of observations per subject
  Max_obs_per_subj = max(Obs_each_subj)
  
  #Initialise objects to store data for each subject
  real_data = matrix(0,nrow = Nsubj,ncol= Max_obs_per_subj*4)
  int_data = matrix(0,nrow  = Nsubj,ncol= Max_obs_per_subj+2)
  
  #Loop through subjects filling in real and int data matricies
  for(subj in 1:Nsubj){
    
    #select predictor variables
    tmp = filter(dat,subject==subj) %>%
      select(m_a,m_b,d_a,d_b)
    
    #create matrix of 999s to pad rows
    tmp2 = matrix(999,nrow=Max_obs_per_subj-nrow(tmp),4)
    
    #create matrix of data and combine with padding
    tmp3 = rbind(data.matrix(tmp,rownames.force = NA),tmp2)
    
    #reshape matrix to vector and store
    real_data[subj,] = as.vector(tmp3) #stacks variables on top of each other.
    
    #get choices for relevant subject
    choice = filter(dat,subject==subj) %>% .$choose_a
    
    #get vector of 999s to pad choices
    choice_padded = c(choice,rep(999,Max_obs_per_subj - Obs_each_subj[subj]))
    
    #store integer
    int_data[subj,1] = Obs_each_subj[subj]
    int_data[subj,2] = Max_obs_per_subj
    int_data[subj,3:(Max_obs_per_subj+2)] = choice_padded
    
  }
  
  #save relevant variables into rdump file
  stan_rdump(c('Nsubj','Max_obs_per_subj','real_data','int_data'),
             file=paste0('data/clean/data_list_e',exp,'_rdump.R'))
  
  stan_list = list(Nsubj = Nsubj,
                   Max_obs_per_subj = Max_obs_per_subj,
                   real_data = real_data,
                   int_data = int_data)
  
  save(stan_list, file=paste0('data/clean/data_list_e',exp,'.RData'))
  
}

#########################################################
# Check bounds on reasonable parameters for k and sigma #
#########################################################

# dat=readRDS(file="data/raw/delay_data_lba_2019.rds") %>%
#   mutate(subject = as.numeric(as.factor(subject)))
# 
# #format data for model fitting
# dat$m_a = dat$amount
# dat$d_a = dat$delay
# dat$m_b = 100
# dat$d_b = 0    #note: b is smaller, sooner option
# dat$choose_a = dat$Choice
# 
# k = 0.1
# sigma = 0.01
# 
# filter(dat,subject==1) %>%
#   mutate(u_a = m_a*exp(-k*d_a),
#          u_b = m_b*exp(-k*d_b),
#          p_a = 1/(1+exp(-sigma*(u_a-u_b)))) %>%
#   ggplot() +
#   geom_line(aes(x=d_a,y=p_a,group=m_a,colour=m_a))




#DEBUG MODELS 
data_list = list(Nsubj = Nsubj,
                 Max_obs_per_subj = Max_obs_per_subj,
                 real_data = real_data,
                 int_data = int_data)

tic = Sys.time()
Sys.setenv(STAN_NUM_THREADS=7)
fit_leow=stan(file=paste0("models/loewenstein1992_gamma_mpi.stan"),
              data=data_list,
              iter=4000,
              cores=1,
              chains=4,
              control=list(max_treedepth=10,adapt_delta=0.90),
              seed=12345)
toc = Sys.time()-tic
print(toc)

fit_leow=stan(file=paste0("models/loewenstein1992_gamma_reparm_mpi.stan"),
              data=data_list,
              iter=200,
              cores=4,
              chains=4,
              control=list(max_treedepth=10,adapt_delta=0.80),
              seed=12345)


fit_nt=stan(file=paste0("models/hyperbolic_nontruncated_mpi.stan"),
         data=data_list,
         iter=2000,
         cores=4,
         chains=4,
         seed=12345)

t = 1:100
alpha = 1/10
beta = 1

y = (1+alpha*t)^-(beta/alpha)
plot(t,y)

#model fails miserably when no truncation is used. Chains don't move
#at all. How does it go when we use truncation?

fit=stan(file=paste0("models/hyperbolic_mpi.stan"),
         data=data_list,
         iter=2000,
         cores=4,
         chains=4,
         seed=12345)

# Failed miserably - but tried again also truncating k_raw to be positive.

fit_gamma=stan(file=paste0("models/hyperbolic_gamma_mpi.stan"),
         data=data_list,
         iter=2000,
         cores=4,
         chains=4,
         seed=12345)

#gamma model does well. But some of the chains aren't finding the typical set until later on.

fit_lnorm=stan(file=paste0("models/hyperbolic_lognormal_mpi.stan"),
               data=data_list,
               iter=2000,
               cores=4,
               chains=4,
               seed=12345)

#lognormal model doesn't do well because it allows seriously high values of l (e.g., 10^5).

#gamme model with doubled iterations
fit_gamma=stan(file=paste0("models/hyperbolic_gamma_mpi.stan"),
               data=data_list,
               iter=4000,
               cores=4,
               chains=4,
               seed=12345)
#excellent convergence, but max_depth exceeded heaps.

fit_gamma=stan(file=paste0("models/hyperbolic_gamma_mpi.stan"),
               data=data_list,
               iter=2000,
               cores=4,
               chains=4,
               control=list(max_treedepth=20,adapt_delta=0.9),
               seed=12345)

#actually got more divergences! Must be because of the shorter burnin.

#delta down to 0.85
fit_gamma=stan(file=paste0("models/hyperbolic_gamma_mpi.stan"),
               data=data_list,
               iter=2000,
               cores=4,
               chains=4,
               control=list(max_treedepth=20,adapt_delta=0.85),
               seed=12345)



#it appears model converges well when delta = 0.85

#5 Oct 2019 - run the harder models with delta = 0.85

fit_gamma80=stan(file=paste0("models/hyperbolic_gamma_mpi.stan"),
               data=data_list,
               iter=4000,
               cores=4,
               chains=4,
               control=list(max_treedepth=20,adapt_delta=0.8),
               seed=12345)

save(fit_gamma80,file="data/derived/fits/test_e2_hyperbolic_gamma_80.RData")
#Not great, chains reach slightly different maximimums. One chain not pulled in until later on.

fit_gamma85=stan(file=paste0("models/hyperbolic_gamma_mpi.stan"),
               data=data_list,
               iter=4000,
               cores=4,
               chains=4,
               control=list(max_treedepth=20,adapt_delta=0.85),
               seed=12345)

save(fit_gamma85,file="data/derived/fits/test_e2_hyperbolic_gamma_85.RData")
#This one converges well.

fit_exp=stan(file=paste0("models/exponential_gamma_mpi.stan"),
               data=data_list,
               iter=4000,
               cores=4,
               chains=4,
               control=list(max_treedepth=20,adapt_delta=0.85),
               seed=12345)

save(fit_exp,file="data/derived/fits/test_e2_exponential_gamma_85.RData")
#This one fits quite poorly. Lots of divergences. Resembles hyper_gamma80. Probably need to increase delta.

fit_hyp_gm=stan(file=paste0("models/hyperbolic_gm_gamma_mpi.stan"),
               data=data_list,
               iter=4000,
               cores=4,
               chains=4,
               control=list(max_treedepth=20,adapt_delta=0.85),
               seed=12345)

save(fit_hyp_gm,file="data/derived/fits/test_e2_hyperbolic_gm_gamma_85.RData")
#Pretty poor convergence. Chains reach different lp's. Fewer divergences, but still bad.

fit_cs=stan(file=paste0("models/const_sens_mpi.stan"),
               data=data_list,
               iter=4000,
               cores=4,
               chains=4,
               control=list(max_treedepth=20,adapt_delta=0.85),
               seed=12345)

save(fit_cs,file="data/derived/fits/test_e2_const_sens_85.RData")
#Chains take awhile to converge. Reasonable number of divergences. No real 
#correlations in parameters. But alpha_mean is piling up around lower bound.

fit_leow=stan(file=paste0("models/loewenstein1992_mpi.stan"),
               data=data_list,
               iter=4000,
               cores=4,
               chains=4,
               control=list(max_treedepth=20,adapt_delta=0.85),
               seed=12345)

save(fit_leow,file="data/derived/fits/test_e2_loewenstein1992_85.RData")
#pretty much the same as cs (above), but no real divergences. 
#alpha_mean and beta_mean (to a lesser extent) pile up around zero.
#I think this is a scaling issue build into the Loewenstein model.
#Alpha and beta are contingent on one another. This might be one of
#those issues that we can't get around. Unless, perhaps we could estimate
#the ratio of alpha to beta?

##################################################################################
### 7 October, 2019. Run the models above but with delta of 0.9 instead of
### 0.85. A lot of these issues might be solved by tuning the sampler. If not,
### we might have to look into separate tunings for each model.



fit_gamma85=stan(file=paste0("models/hyperbolic_gamma_mpi.stan"),
                 data=data_list,
                 iter=4000,
                 cores=4,
                 chains=4,
                 control=list(max_treedepth=20,adapt_delta=0.9),
                 seed=12345)

save(fit_gamma,file="data/derived/fits/test_e2_hyperbolic_gamma_90.RData")
print(fit_gamma)

fit_exp=stan(file=paste0("models/exponential_gamma_mpi.stan"),
             data=data_list,
             iter=4000,
             cores=4,
             chains=4,
             control=list(max_treedepth=20,adapt_delta=0.90),
             seed=12345)

save(fit_exp,file="data/derived/fits/test_e2_exponential_gamma_90.RData")
print(fit_exp)


fit_hyp_gm=stan(file=paste0("models/hyperbolic_gm_gamma_mpi.stan"),
                data=data_list,
                iter=4000,
                cores=4,
                chains=4,
                control=list(max_treedepth=20,adapt_delta=0.90),
                seed=12345)

save(fit_hyp_gm,file="data/derived/fits/test_e2_hyperbolic_gm_gamma_90.RData")
print(fit_hyp_gm)

fit_cs=stan(file=paste0("models/const_sens_mpi.stan"),
            data=data_list,
            iter=4000,
            cores=4,
            chains=4,
            control=list(max_treedepth=20,adapt_delta=0.90),
            seed=12345)

save(fit_cs,file="data/derived/fits/test_e2_const_sens_90.RData")
print(fit_cs)

fit_leow=stan(file=paste0("models/loewenstein1992_mpi.stan"),
              data=data_list,
              iter=4000,
              cores=4,
              chains=4,
              control=list(max_treedepth=20,adapt_delta=0.90),
              seed=12345)

save(fit_leow,file="data/derived/fits/test_e2_loewenstein1992_90.RData")
print(fit_leow)




#test to make sure no coding errors in models
fit_gamma85=stan(file=paste0("models/killeen2009_gamma_mpi.stan"),
                 data=data_list,
                 iter=1,
                 cores=1,
                 chains=1)


hist(rinvgamma(10000,shape=1,scale=.1))

shape = 1
scale = 1
hist(rgamma(10000,shape,scale))

shape = 1
scale = 5
hist(rgamma(10000,shape,1/scale))


a = 5
b = 9

mean = a / (a+b)
prec = a+b

#a = prec - b

#mean = a/(prec)




a <- mean * prec 
b <- (1-mean) * prec


a = 5
b = 3
eps = 0.001
a^eps / (a^eps + b^eps)

# save(data_list_mpi, file = "data/derived/data_list_mpi.R")

#prep data
data_list = list(Ntotal = dim(dat)[1],
                 Nsubj = Nsubj,
                 subj = as.numeric(dat$subject),
                 y=dat$choose_a,
                 d_a = dat$d_a,
                 d_b = dat$d_b,
                 m_a = dat$m_a,
                 m_b = dat$m_b)

source("models/model_details.R")

i = 10

set.seed(12345)

fit=stan(file=paste0("models/",names(models)[i],"_PL.stan"),
         data=data_list,
         iter=10,
         cores=1,
         chains=1,
         seed=12345)


fit_mpi=stan(file=paste0("models/",names(models)[i],"_PL_mpi.stan"),
         data=data_list,
         iter=10,
         cores=1,
         chains=1,
         seed=12345)
# 
# 
# 
# fit=stan(file=paste0("models/",names(models)[i],".stan"),
#          data=data_list,
#          iter=2000,
#          warmup=1000,
#          cores=8,
#          chains=8,
#          control=list(max_treedepth=20,adapt_delta=0.99),
#          seed=12345)
# 
# save(fit,file=paste0("data/derived/fits/manos_fit_Bayes_",models[i],".RData"))


