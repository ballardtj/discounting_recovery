library(rstan)
options(mc.cores = parallel::detectCores())

args <- commandArgs(trailingOnly = T)

print(args)

run <- as.numeric(  gsub("--","",args[1]) )

#identify chain, replication, and model from run number 
n_chains = 4
n_reps = 100
n_models = 11

design = expand.grid(chain = 1:n_chains,
                     rep = 1:n_reps,
                     model = 1:n_models,
                     exp = 1:2)

print(paste0("Exp: ",design$exp[run],
             "  Model: ",design$model[run],
             "  Replication: ",design$rep[run],
             "  Exp: ",design$exp[run]))

#set directory
dir = "/30days/uqtballa/Q0992"

#source model details
source(paste0(dir,"/models/model_details.R"))
model = names(models)[design$model[run]]

#load data list
load(paste0(dir,"/data/derived/simulated/data_list_exp",design$exp[run],'_',model,'_',design$rep[run],'.RData'))

#set number of threads
Sys.setenv(STAN_NUM_THREADS=11)

#set target acceptance rate
delta = 0.99
if( design$exp[run] == 2 & design$model[run] == 10  ) {delta = 0.85} #Mcclure2007 wouldn't finish for E2 with delta = 0.9

#set r seed

#call stan
tic = Sys.time()

set.seed(12345)

fit=stan(file=paste0(dir,"models/",model,"_PL_mpi.stan"),
              data=stan_list,
              iter=4,
              cores=1,
              chains=1,
              control=list(max_treedepth=20,adapt_delta=delta),
              seed=12345,
              chain_id = design$chain[run])

toc = Sys.time()-tic

print(toc)

#save fit object
save(fit,file=paste0(dir,"/analysis/output/fit_recovered_exp",design$exp[run],"_",model,"_",design$rep[run],"_",design$chain[run],".RData"))
