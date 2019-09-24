rm(list=ls())

###read in the arguments listed at the command line
args <- commandArgs(trailingOnly = F)

print(args)

i <- args[length(args)]
i <- strsplit(i,"--")[[1]][2]
i <- as.numeric(i)

print(i)

j = i %% 100
i = ceiling( i / 100 )

if(j==0){j=100}

#set working directory
setwd("/QRISdata/Q0992")

library(R.matlab)
library(plyr)
library(rstan)

#load data and transform into R data frame (a = LL, b = SS)
tmp = readMat("data/raw/datS3.mat")
tmp = tmp[[1]]
dat = adply(tmp,c(2,1))
names(dat) = c("trial","subject","m_b","d_b","m_a","d_a","choose_a")

dat1 = filter(dat,subject==1)

parms = expand.grid(rate = seq(0,5,by=0.1),
                    sigma = seq(0,5,by=0.1),
                    mean_prob_a = NA)

for(i in 1:nrow(parms)){
  util_a = dat1$m_a / 1+parms$rate[i]*dat1$d_a
  util_b = dat1$m_b / 1+parms$rate[i]*dat1$d_b
  prob_a = 1/(1+exp(-parms$sigma[i]*(util_a-util_b)))
  parms$mean_prob_a[i] = mean(prob_a)
}

ggplot(parms,aes(x=rate,y=mean_prob_a)) +
  geom_line() + facet_wrap(~sigma)

x  = rlnorm(10000,0,1)
plot(density(x))
quantile(x)

x = rcauchy(10000,0,1)
plot(density(x))
quantile(x)

x = rt(10000,1)



dat %>% filter(subject == 1) %>%
  summarise(prop_LL = mean(choose_a)) %>%
  ggplot(aes(x=y=prop_LL)) + geom_point()


dat %>% mutate(m_diff = m_a - m_b,
               d_diff = d_a - d_b )  %>%
  group_by(m_diff,d_diff) %>%
  summarise(prob_a = mean(choose_a)) %>%
  ggplot() +
  geom_line(aes(x=m_diff,y=prob_a,colour=as.factor(d_diff)))

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

set.seed(12345)
nsamples = 100  #number of different parameter values tested
nreps = 100     #number of times recovery is run for each parameter value
iter_samples = sample(1:8000,size=nsamples,replace=TRUE)
subject_samples = sample(1:length(unique(dat$subject)),size=nsamples,replace=TRUE)

#load likelihood
source(paste0("models/",names(models)[i],".R"))

#load fit object
load(file=paste0("data/derived/fits/junyi_fit_Bayes_", names(models)[i],".RData"))
  
#extract parameters
posts=rstan::extract(fit)
rm(fit)
  
#parameters
parm_names = models[names(models)[i]][[1]]$parms
parm_names
names(posts)

#preallocate matrix of samples
parms = matrix(NA,nrow=nsamples,ncol=length(parm_names))
colnames(parms) = parm_names

#for(j in 1:nsamples){
  
  #store sampled parameter values 
  for(k in 1:length(parm_names)){
  #for each sample re-centre parameters (will need to be updated for different models)
  if( (names(models)[i] == "mcclure2007" & parm_names[k] == "omega") | 
      (names(models)[i] == "kileen2009" & parm_names[k] == "alpha") | 
      names(models)[i] == "kileen2009" & parm_names[k] == "beta")  {
    #some parameters don't need to be transformed. So the value can be read directly from the posts matrix
    parms[j,k] = posts[parm_names[k]][[1]][iter_samples[j],subject_samples[j]]
  } else {
    #all other parameters need to be transformed from their raw values to the uncentered version
    parms[j,k] = posts[paste0(parm_names[k],'_raw')][[1]][iter_samples[j],subject_samples[j]]*
      posts[paste0(parm_names[k],'_sd')][[1]][iter_samples[j]] +
      posts[paste0(parm_names[k],'_mean')][[1]][iter_samples[j]]
  }
}

  #get data for relevant participant
  sim_dat = filter(dat,subject==subject_samples[i])

  #get choice probabilities under sampled parameter values
  sim_dat$prob_a = likelihood(parms[j,],sim_dat)
  
  #for each repetition, generate random choices from probabilities
  sim_list = list()
  for(l in 1:nreps){
    #randomly generate choices
    sim_dat$y = as.numeric( runif(nrow(sim_dat)) < sim_dat$prob_a )
    sim_dat$subj = l #here, each replication is a different 'subject'. That way, we can fit all replications at once within the same model
    sim_list[[l]] = sim_dat
  }
  stan_list = as.list(bind_rows(sim_list))
  stan_list$Ntotal = length(stan_list$y)
  stan_list$Nsubj = nreps
  
  #run person level model in stan with each value of 'subj' representing a seperate repetition.
  fit = stan(file=paste0('models/',names(models)[i],'_PL.stan'),data=stan_list,cores=4,seed=12345)
  #save results of person level model
  #save(fit,file=paste0("data/derived/recovery_fit_",names(models)[i],"_",j,".RData"))
  save(fit,file=paste0("data/derived/recovery/junyi_recovery_fit_",names(models)[i],"_",j,".RData"))

#}

#save generating parameters
#save(parms,file=paste0("data/derived/recovery_generating_parms_",names(models)[i],".RData"))
save(parms,file=paste0("data/derived/recovery/junyi_recovery_generating_parms_",names(models)[i],"_",j,".RData"))


#tag trials according to the effects they are used to generate
# dat = dat %>%
#   mutate(trial = as.numeric(trial),
#          effect = case_when(
#            #delay duration effect (delay for LL is twice as long as delay for SS)
#            #magnitude of b (SS) varied between participants. Magnitude of a (LL) always 36.
#            trial <= 100 ~ 1, 
#            #common difference effect (magnitude for LL was always 30 more than SS).
#            #magnitude of a (LL) always 32. magnitude of b (SS) varied across participants. Delay for a (SS) and (LL) varied,
#            #but were constant across participants.
#            trial > 100 & trial <=200 ~ 2,
#            #magnitude effect (magnitude of LL always twice magnitude of SS). Magnitude of b (SS) varied from
#            #2 - 40 in steps of 2. Delay of b (SS) fixed at 12. Delay of a (LL) varied across subjects.
#            trial > 200 ~ 3
#          )
#   )
# 


       
       
       

