library(rstan)

args <- commandArgs(trailingOnly = T)

print(args)

run <- as.numeric(  gsub("--","",args[1]) )

print(paste0("Run: ",run))

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

#load individual stanfit objects from each chain
sflist = list()
for(chain in 1:4){
  load(paste0(dir,"/analysis/output/fit_exp",exp,"_",model,"_",chain,".RData"))
  sflist[[chain]] = fit
  rm(fit)
}

fit = sflist2stanfit(sflist) 

#save fit object
save(fit,file=paste0(dir,"/data/derived/fits/fit_exp",exp,"_",model,".RData"))
