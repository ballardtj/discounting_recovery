library(rstan)

args <- commandArgs(trailingOnly = T)

print(args)

i <- args[1]

run <- as.numeric( strsplit(i[1],"--")[[1]][2] )

source("/QRISdata/Q0992/models/model_details.R")

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


regexstr = paste0("samples_e",exp,"_",model,"_.?.csv")

path_to_model = "/QRISdata/Q0992/analysis/Bayes/"

#get relevant csv files
csvfiles=dir(path=path_to_model,pattern=regexstr)

print(csvfiles)

#create fit object from those csv files
library(rstan)
fit=read_stan_csv(paste0(path_to_model,csvfiles))

#save fit object
save(fit,file=paste0("/QRISdata/Q0992/data/derived/fits/fit_exp",exp,"_",model,".RData"))
