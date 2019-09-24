rm(list=ls())

#load packages
library(plyr)
library(tidyverse)
library(R.matlab)
library(rstan)

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
  
}

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
