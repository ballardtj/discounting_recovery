rm(list=ls())

#source libraries
library(tidyverse)
library(DEoptim)
library(GGally)

#load data
dat=readRDS(file="../../data/data_delay_lba2018.rds")

########################
##### Experiment 1 #####
########################

#Set up dataset for Experiment 1

dat=expand.grid(ll_m=600,ll_d=c(1,10,30,50,100,150),
                ss_m=seq(20,600,20),ss_d=0)

#--------------------
#HYPERBOLIC MODEL

hyperbolic.likelihood=function(ll_m,ll_d,ss_m,ss_d,
                             k,sigma){
  #Equation 5
  uLL = ll_m/(1+k*ll_d)
  uSS = ss_m/(1+k*ss_d)
  pLL = pnorm((uLL-uSS)/sigma)
  return(pLL)
}

hyperbolic.wrapper=function(pars,dat){
  #extract pars
  k=pars[1]
  sigma=pars[2]

  pLLs = hyperbolic.likelihood(dat$ll_m,dat$ll_d,dat$ss_m,dat$ss_d,
                             k,sigma)
  pLLs = pmin(pmax(pLLs,0.001),0.999)  
  neglnLs = -log(dat$chooseLL*pLLs + (1-dat$chooseLL)*(1-pLLs)) 
  return(sum(neglnLs))
}

#data generating values
Nparms=100
Nsims=100
ks=runif(Nparms)          #bounded between 0 and 1    #0.025
sigmas=runif(Nparms)*100   #bounded between 0 and 100

##initialise lists to store results online
recov_list1=list() #list the stores data from outer loop
recov_list2=list() #stores data from inner loop

#initialise data frame for simulations
recov_tmp2=data.frame(parname=c('k','sigma'),
                     generating=NA,
                     estimated=NA,
                     sim=NA)

pb <- txtProgressBar(min = 0, max = Nparms*Nsims, style = 3)
ctr=0
for(p in 1:Nparms){
 
  #set generating parms
  k=ks[p]
  sigma=sigmas[p]
  
  #calculate choice probabilities under data generating parms
  dat_parm = dat %>% 
    mutate(pLL=hyperbolic.likelihood(ll_m,ll_d,ss_m,ss_d,k,sigma))

for(i in 1:Nsims){
  #Generate one choice for each stimulus
  dat_tmp = dat_parm %>% mutate(chooseLL = (pLL > runif(length(pLL)))*1)
  
  #Maximum likelihood optimisation
  output=DEoptim(fn=hyperbolic.wrapper,
                 dat=dat_tmp,
                 lower=c(0,0),
                 upper=c(1000,1000),
                 control=list(trace=0))
  
  recov_tmp2$generating = c(k,sigma)
  recov_tmp2$estimated = output$optim$bestmem
  recov_tmp2$sim=i
  recov_list2[[i]]=recov_tmp2
  
  ctr=ctr+1
  setTxtProgressBar(pb, ctr)
}


recov_tmp1=bind_rows(recov_list2) %>%
  group_by(parname) %>%
  summarise(generating = mean(generating),
            estimated_mean = mean(estimated),
            estimated_upper = estimated_mean + sd(estimated)/sqrt(length(estimated)),
            estimated_lower = estimated_mean - sd(estimated)/sqrt(length(estimated)))

recov_list1[[p]]=recov_tmp1

}

recov = bind_rows(recov_list1)

save(recov,file="../output/Cheng_e1_hyperbolic_recovery.Rdata")


#Plots
load(file="../output/Cheng_e1_hyperbolic_recovery.Rdata")

ggplot(data=recov) +
  geom_point(aes(x=generating,y=estimated_mean)) +
  geom_errorbar(aes(x=generating,ymin=estimated_lower,ymax=estimated_upper)) +
  facet_wrap(~parname,scale="free") +
  geom_abline(aes(intercept=0,slope=1))
ggsave(file="../output/Cheng_e1_hyperbolic_recovery_error.pdf",height=6,width=8)




pd=recov %>%
 select(parname,estimated_mean) %>%
  mutate(rep = rep(1:Nparms,each=length(unique(parname)))) %>%   
  spread(parname,estimated_mean) %>%
  select(-rep)

ggpairs(data = pd)
ggsave(file="../output/Cheng_e1_hyperbolic_recovery_pairs.pdf",height=8,width=8)



#---------------------------------
#PROPORTIONAL DIFFERENCE MODEL

propdiff.likelihood=function(ll_m,ll_d,ss_m,ss_d,
                               delta,sigma){
  #Equation 1
pi_m = (pmax(abs(ll_m),abs(ss_m)) - pmin(abs(ll_m),abs(ss_m))) / pmax(abs(ll_m),abs(ss_m)) #max minus min of each pair divided by max of all options
pi_d = (pmax(abs(ll_d),abs(ss_d)) - pmin(abs(ll_d),abs(ss_d))) / pmax(abs(ll_d),abs(ss_d))
d = pi_m-pi_d
pLL = pnorm((d-delta) / sigma)
  return(pLL)
}

propdiff.wrapper=function(pars,dat){
  #extract pars
  delta=pars[1]
  sigma=pars[2]
  
  pLLs = hyperbolic.likelihood(dat$ll_m,dat$ll_d,dat$ss_m,dat$ss_d,
                               delta,sigma)
  pLLs = pmin(pmax(pLLs,0.001),0.999)  
  neglnLs = -log(dat$chooseLL*pLLs + (1-dat$chooseLL)*(1-pLLs)) 
  return(sum(neglnLs))
}


#data generating values
Nparms=100
Nsims=100
deltas=runif(Nparms,min=-10,max=10)          #bounded between 0 and 1    #0.025
sigmas=runif(Nparms)*100   #bounded between 0 and 100

##initialise lists to store results online
recov_list1=list() #list the stores data from outer loop
recov_list2=list() #stores data from inner loop

#initialise data frame for simulations
recov_tmp2=data.frame(parname=c('delta','sigma'),
                      generating=NA,
                      estimated=NA,
                      sim=NA)

pb <- txtProgressBar(min = 0, max = Nparms*Nsims, style = 3)
ctr=0
for(p in 1:Nparms){
  
  #set generating parms
  delta=deltas[p]
  sigma=sigmas[p]
  
  #calculate choice probabilities under data generating parms
  dat_parm = dat %>% 
    mutate(pLL=propdiff.likelihood(ll_m,ll_d,ss_m,ss_d,delta,sigma))
  
  for(i in 1:Nsims){
    #Generate one choice for each stimulus
    dat_tmp = dat_parm %>% mutate(chooseLL = (pLL > runif(length(pLL)))*1)
    
    #Maximum likelihood optimisation
    output=DEoptim(fn=propdiff.wrapper,
                   dat=dat_tmp,
                   lower=c(-1000,0),
                   upper=c(1000,1000),
                   control=list(trace=0))
    
    recov_tmp2$generating = c(delta,sigma)
    recov_tmp2$estimated = output$optim$bestmem
    recov_tmp2$sim=i
    recov_list2[[i]]=recov_tmp2
    
    ctr=ctr+1
    setTxtProgressBar(pb, ctr)
  }
  
  
  recov_tmp1=bind_rows(recov_list2) %>%
    group_by(parname) %>%
    summarise(generating = mean(generating),
              estimated_mean = mean(estimated),
              estimated_upper = estimated_mean + sd(estimated)/sqrt(length(estimated)),
              estimated_lower = estimated_mean - sd(estimated)/sqrt(length(estimated)))
  
  recov_list1[[p]]=recov_tmp1
  
}

recov = bind_rows(recov_list1)

save(recov,file="../output/Cheng_e1_propdiff_recovery.Rdata")
load(file="../output/Cheng_e1_propdiff_recovery.Rdata")

ggplot(data=recov) +
  geom_point(aes(x=generating,y=estimated_mean)) +
  geom_errorbar(aes(x=generating,ymin=estimated_lower,ymax=estimated_upper)) +
  facet_wrap(~parname,scale="free") +
  geom_abline(aes(intercept=0,slope=1))
ggsave(file="../output/Cheng_e1_propdiff_recovery_error.pdf",height=6,width=8)



pd=recov %>%
  select(parname,estimated_mean) %>%
  mutate(rep = rep(1:Nparms,each=length(unique(parname)))) %>%   
  spread(parname,estimated_mean) %>%
  select(-rep)

ggpairs(data = pd)
ggsave(file="../output/Cheng_e1_propdiff_recovery_pairs.pdf",height=8,width=8)



#--------------------
#TRADEOFF MODEL


tradeoff.likelihood=function(ll_m,ll_d,ss_m,ss_d,
                                 gamma,eps,tau,theta,kappa,alpha){
  #Equation 5
  term1 = ( (1/gamma)*log(1+gamma*ll_m) - (1/gamma)*log(1+gamma*ss_m) )^(1/eps) 
  term2 = ( ( (1/tau)*log(1+tau+ll_d) - (1/tau)*log(1+tau*ss_d) )/theta )^theta
  term3 = (kappa/alpha)*log(1+alpha*term2)^(1/eps)
  pLL = term1/(term1+term3)
  
  #This happens at some parameter values and causes NaN
  pLL[term1==0 & term3==0] <- 0.5
  
  return(pLL)
}

tradeoff.wrapper=function(pars,dat){

  #extract pars
  gamma=pars[1]
  eps=pars[2]
  tau=pars[3]
  theta=pars[4]
  kappa=pars[5]
  alpha=pars[6]
 
  pLLs = tradeoff.likelihood(dat$ll_m,dat$ll_d,dat$ss_m,dat$ss_d,
            gamma,eps,tau,theta,kappa,alpha)
  pLLs = pmin(pmax(pLLs,0.001),0.999)  
  neglnLs = -log(dat$chooseLL*pLLs + (1-dat$chooseLL)*(1-pLLs)) 
 
  return(sum(neglnLs))
}
  

#data generating values
Nparms=100
Nsims=100
gammas=runif(Nparms)*10          #bounded between 0 and 10    
epss=runif(Nparms)*10          #bounded between 0 and 10    
taus=runif(Nparms)*10          #bounded between 0 and 10    
thetas=runif(Nparms)*10          #bounded between 0 and 10    
kappas=runif(Nparms)*10          #bounded between 0 and 10    
alphas=runif(Nparms)          #bounded between 0 and 1    


##initialise lists to store results online
recov_list1=list() #list the stores data from outer loop
recov_list2=list() #stores data from inner loop

#initialise data frame for simulations
recov_tmp2=data.frame(parname=c('gamma','eps','tau','theta','kappa','alpha'),
                      generating=NA,
                      estimated=NA,
                      sim=NA)

pb <- txtProgressBar(min = 0, max = Nparms*Nsims, style = 3)
ctr=0
for(p in 1:Nparms){
  
  #set generating parms
  gamma=gammas[p]
  eps=epss[p]
  tau=taus[p]
  theta=thetas[p]
  kappa=kappas[p]
  alpha=alphas[p]
  
  #calculate choice probabilities under data generating parms
  dat_parm = dat %>% 
    mutate(pLL=tradeoff.likelihood(ll_m,ll_d,ss_m,ss_d,
                                   gamma,eps,tau,theta,kappa,alpha))
  
  for(i in 1:Nsims){
    #Generate one choice for each stimulus
    dat_tmp = dat_parm %>% mutate(chooseLL = (pLL > runif(length(pLL)))*1)
    
    #Maximum likelihood optimisation
    output=DEoptim(fn=tradeoff.wrapper,
                   dat=dat_tmp,
                   lower=c(0,0,0,0,0,0),
                   upper=c(1000,1000,1000,1000,1000,1000),
                   control=list(trace=0))
    
    recov_tmp2$generating = c(gamma,eps,tau,theta,kappa,alpha)
    recov_tmp2$estimated = output$optim$bestmem
    recov_tmp2$sim=i
    recov_list2[[i]]=recov_tmp2
    
    ctr=ctr+1
    setTxtProgressBar(pb, ctr)
  }
  
  
  recov_tmp1=bind_rows(recov_list2) %>%
    group_by(parname) %>%
    summarise(generating = mean(generating),
              estimated_mean = mean(estimated),
              estimated_upper = estimated_mean + sd(estimated)/sqrt(length(estimated)),
              estimated_lower = estimated_mean - sd(estimated)/sqrt(length(estimated)))
  
  recov_list1[[p]]=recov_tmp1
  
}

recov = bind_rows(recov_list1)

save(recov,file="../output/Cheng_e1_tradeoff_recovery.Rdata")

load(file="../output/Cheng_e1_tradeoff_recovery.Rdata")


ggplot(data=recov) +
  geom_point(aes(x=generating,y=estimated_mean)) +
  geom_errorbar(aes(x=generating,ymin=estimated_lower,ymax=estimated_upper)) +
  facet_wrap(~parname,scale="free") +
  geom_abline(aes(intercept=0,slope=1))
ggsave(file="../output/Cheng_e1_tradoff_recovery_error.pdf",height=6,width=8)

pd=recov %>%
  select(parname,estimated_mean) %>%
  mutate(rep = rep(1:Nparms,each=length(unique(parname)))) %>%   
  spread(parname,estimated_mean) %>%
  select(-rep)
ggpairs(data = pd)
ggsave(file="../output/Cheng_e1_tradoff_recovery_pairs.pdf",height=8,width=8)

#5) May need to do ITCH model at some stage

