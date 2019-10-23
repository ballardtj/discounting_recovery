library(rstan)

args <- commandArgs(trailingOnly = T)

print(args)

run <- as.numeric(  gsub("--","",args[1]) )

#identify chain, replication, and model from run number 
n_reps = 100
n_models = 11

design = expand.grid(rep = 1:n_reps,
                     model = 1:n_models,
                     exp = 1:2)

print(paste0("Exp: ",design$exp[run],
             "  Model: ",design$model[run],
             "  Replication: ",design$rep[run],
             "  Exp: ",design$exp[run]))

#set directory
#dir = "~/Nextcloud/DISCRECOV-Q0992"
dir = "/30days/uqtballa/Q0992"

#source model details
source(paste0(dir,"/models/model_details.R"))
model = names(models)[design$model[run]]

#load individual stanfit objects from each chain
sflist = list()
for(chain in 1:4){
  load(paste0(dir,"/analysis/output/fit_recovered_exp",design$exp[run],"_",model,"_",design$rep[run],"_",chain,".RData"))
  sflist[[chain]] = fit
  rm(fit)
}

fit = sflist2stanfit(sflist) 

#save fit object
save(fit,file=paste0(dir,"/data/derived/recovery/fit_recovery_exp",design$exp[run],"_",model,"_",design$rep[run],".RData"))
