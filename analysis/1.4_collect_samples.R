library(rstan)

args <- commandArgs(trailingOnly = T)

print(args)

run <- as.numeric(  gsub("--","",args[1]) )
chain <- as.numeric(  gsub("--","",args[2]) )

print(paste0("Run: ",run,"; Chain: ",chain))

#set directory
#dir = "~/Nextcloud/DISCRECOV-Q0992"
dir = "/30days/uqtballa/Q0992"

#source model details
source(paste0(dir,"/models/model_details.R"))

if(run < 12){
  exp = 1
  model = names(models)[run]
  if(run != 4 & run != 6){model = paste0(model,"_gamma")}
}

if(run > 11){
  exp = 2
  model = names(models)[run-11]
  if(run != 15 & run != 16){model = paste0(model,"_gamma")}
}

#load data list
load(paste0(dir,"/data/clean/data_list_e",exp,'.RData'))

#call stan
tic = Sys.time()

Sys.setenv(STAN_NUM_THREADS=11)

fit=stan(file=paste0(dir,"/models/",model,"_mpi.stan"),
              data=stan_list,
              iter=4000,
              cores=1,
              chains=1,
              control=list(max_treedepth=20,adapt_delta=0.90),
              seed=12345,
              chain_id = chain)

toc = Sys.time()-tic

print(toc)

#save fit object
save(fit,file=paste0(dir,"/analysis/output/fit_exp",exp,"_",model,"_",chain,".RData"))
