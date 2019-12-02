#clear workspace
rm(list=ls())

#load libraries
library(tidyverse)
library(R.matlab)
library(plyr)
library(rstan)

#load model details
source("models/model_details.R")

#set analysis parameters
set.seed(12345)
nsamples = 100  #number of different parameter values tested
nreps = 100    #number of times recovery is run for each parameter value
iter_samples = sample(1:8000,size=nsamples,replace=TRUE) #vector of sampled iterations

#loop through experiments
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
    sim_dat = filter(dat,subject==1) #all p's had same design, so just need one p's worth of data
  }
  
  if(exp == 2){
    #load data and transform into R data frame (a = LL, b = SS)
    tmp = readMat("data/raw/datS3.mat")
    tmp = tmp[[1]]
    dat = adply(tmp,c(2,1))
    names(dat) = c("trial","subject","m_b","d_b","m_a","d_a","choose_a")
  }
  
  #vector of sampled subjects
  subject_samples = sample(1:length(unique(dat$subject)),size=nsamples,replace=TRUE)
  
  #loop through models
  for(model in 6){#1:length(models)){
    print(model)
    
    #load likelihood
    source(paste0("models/",names(models)[model],".R"))
    
    #load fit object
    load(file=paste0("data/derived/fits99/fit_exp",exp,"_",names(models)[model],".RData"))
    
    #extract parameters
    posts=rstan::extract(fit)
    rm(fit)
    
    #parameters
    parm_names = models[names(models)[model]][[1]]$parms
    parm_names
    names(posts)
    
    #preallocate matrix of samples
    parms = matrix(NA,nrow=nsamples,ncol=length(parm_names))
    colnames(parms) = parm_names
    
    for(sample in 1:nsamples){
      
      #in the second dataset, participants had different items. We extract the items
      #relevant to the participant whose parameter was sampled
      if(exp == 2){
        sim_dat = filter(dat,subject==subject_samples[sample])
      }
      
      #transform (if necessary) and store sampled parameter
      if(names(models)[model] == "hyperbolic"){
        #no transforms necessary
        parms[sample,1] = posts$k[iter_samples[sample],subject_samples[sample]]
        parms[sample,2] = posts$sigma[iter_samples[sample],subject_samples[sample]]
      }
      
      if(names(models)[model] == "exponential"){
        #no transforms necessary
        parms[sample,1] = posts$k[iter_samples[sample],subject_samples[sample]]
        parms[sample,2] = posts$sigma[iter_samples[sample],subject_samples[sample]]
      }
      
      if(names(models)[model] == "hyperbolic_gm"){
        #no transforms necessary
        parms[sample,1] = posts$k[iter_samples[sample],subject_samples[sample]]
        parms[sample,2] = posts$s[iter_samples[sample],subject_samples[sample]]
        parms[sample,3] = posts$sigma[iter_samples[sample],subject_samples[sample]]
      }
      
      if(names(models)[model] == "prop_diff"){
        #delta needs to be unstandardised
        parms[sample,1] = posts$delta_raw[iter_samples[sample],subject_samples[sample]]*
          posts$delta_sd[iter_samples[sample]] + posts$delta_mean[iter_samples[sample]]
        parms[sample,2] = posts$sigma[iter_samples[sample],subject_samples[sample]]
      }
      
      if(names(models)[model] == "tradeoff"){
        #1 is added to theta
        parms[sample,1] = posts$gamma[iter_samples[sample],subject_samples[sample]]
        parms[sample,2] = posts$tau[iter_samples[sample],subject_samples[sample]]
        parms[sample,3] = 1+posts$theta_raw[iter_samples[sample],subject_samples[sample]]
        parms[sample,4] = posts$kappa[iter_samples[sample],subject_samples[sample]]
        parms[sample,5] = posts$alpha[iter_samples[sample],subject_samples[sample]]
        parms[sample,6] = posts$eps[iter_samples[sample],subject_samples[sample]]
      }
      
      if(names(models)[model] == "ITCH"){
        #all parameters need to be unstandardised
        parms[sample,1] = posts$B1_raw[iter_samples[sample],subject_samples[sample]]*
          posts$B1_sd[iter_samples[sample]] + posts$B1_mean[iter_samples[sample]]
        parms[sample,2] = posts$BxA_raw[iter_samples[sample],subject_samples[sample]]*
          posts$BxA_sd[iter_samples[sample]] + posts$BxA_mean[iter_samples[sample]]
        parms[sample,3] = posts$BxR_raw[iter_samples[sample],subject_samples[sample]]*
          posts$BxR_sd[iter_samples[sample]] + posts$BxR_mean[iter_samples[sample]]
        parms[sample,4] = posts$BtA_raw[iter_samples[sample],subject_samples[sample]]*
          posts$BtA_sd[iter_samples[sample]] + posts$BtA_mean[iter_samples[sample]]
        parms[sample,5] = posts$BtR_raw[iter_samples[sample],subject_samples[sample]]*
          posts$BtR_sd[iter_samples[sample]] + posts$BtR_mean[iter_samples[sample]]
      }
      
      if(names(models)[model] == "const_sens"){
        #no transforms necessary
        parms[sample,1] = posts$alpha[iter_samples[sample],subject_samples[sample]]
        parms[sample,2] = posts$beta[iter_samples[sample],subject_samples[sample]]
        parms[sample,3] = posts$sigma[iter_samples[sample],subject_samples[sample]]
      }
      
      if(names(models)[model] == "mazur1987"){
        #no transforms necessary
        parms[sample,1] = posts$k[iter_samples[sample],subject_samples[sample]]
        parms[sample,2] = posts$s[iter_samples[sample],subject_samples[sample]]
        parms[sample,3] = posts$sigma[iter_samples[sample],subject_samples[sample]]
      }
      
      if(names(models)[model] == "loewenstein1992"){
        #beta_on_alpha needs to be multipled by alpha to calculate beta
        parms[sample,1] = posts$alpha[iter_samples[sample],subject_samples[sample]]
        parms[sample,2] = posts$beta_on_alpha[iter_samples[sample],subject_samples[sample]]*
          posts$alpha[iter_samples[sample],subject_samples[sample]]
        parms[sample,3] = posts$sigma[iter_samples[sample],subject_samples[sample]]
      }
      
      if(names(models)[model] == "mcclure2007"){
        #no transform necessary
        parms[sample,1] = posts$omega[iter_samples[sample],subject_samples[sample]]
        parms[sample,2] = posts$beta[iter_samples[sample],subject_samples[sample]]
        parms[sample,3] = posts$delta[iter_samples[sample],subject_samples[sample]]
        parms[sample,4] = posts$sigma[iter_samples[sample],subject_samples[sample]]
      }
      
      if(names(models)[model] == "killeen2009"){
        #no transform necessary
        parms[sample,1] = posts$alpha[iter_samples[sample],subject_samples[sample]]
        parms[sample,2] = posts$beta[iter_samples[sample],subject_samples[sample]]
        parms[sample,3] = posts$lambda[iter_samples[sample],subject_samples[sample]]
        parms[sample,4] = posts$sigma[iter_samples[sample],subject_samples[sample]]
      }
        
      #get choice probabilities under sampled parameter values
      sim_dat$prob_a = likelihood(parms[sample,],sim_dat)
      
      Nsubj = nreps
      nobs = nrow(sim_dat) 
      Max_obs_per_subj = nobs
      
      #Initialise objects to store data for each repetition
      real_data = matrix(0,nrow = nreps,ncol = Max_obs_per_subj*4)
      int_data = matrix(0,nrow  = nreps,ncol = Max_obs_per_subj+2)
      
      #for each repetition, generate random choices from probabilities
      sim_list = list()
      for(rep in 1:nreps){
        
        #select predictor variables
        tmp = select(sim_dat,m_a,m_b,d_a,d_b)
        
        #create matrix of 999s to pad rows (not needed for recovery since all repetitions are the same size)
        tmp2 = matrix(999,nrow=Max_obs_per_subj-nrow(tmp),4)
        
        #create matrix of data and combine with padding
        tmp3 = rbind(data.matrix(tmp,rownames.force = NA),tmp2)
        
        #reshape matrix to vector and store
        real_data[rep,] = as.vector(tmp3) #stacks variables on top of each other.
        
        #simulate choices using choice probabiilty generated by the model
        choice =  as.numeric( runif(nrow(sim_dat)) < sim_dat$prob_a )
        
        #get vector of 999s to pad choices
        choice_padded = c(choice,rep(999,Max_obs_per_subj - nobs))
        
        #store integer
        int_data[rep,1] = nobs
        int_data[rep,2] = Max_obs_per_subj
        int_data[rep,3:(Max_obs_per_subj+2)] = choice_padded
        
      } #close replication loop
      
      #save relevant variables into rdump file
      #stan_rdump(c('Nsubj','Max_obs_per_subj','real_data','int_data'),
      #           file=paste0('data/derived/simulated2/data_list_e',exp,'_',names(models)[model],'_',sample,'_rdump.R'))
      
      stan_list = list(Nsubj = Nsubj,
                       Max_obs_per_subj = Max_obs_per_subj,
                       real_data = real_data,
                       int_data = int_data)
      
      save(stan_list, file=paste0('data/derived/simulated/data_list_exp',exp,'_',names(models)[model],'_',sample,'.RData'))
      
    } #close sample loop
    
    #save generating parameters
    save(parms,file=paste0("data/derived/simulated/generating_parms_exp",exp,"_",names(models)[model],".RData"))
    
    
    
  } #close model loop
}#close experiment loop

