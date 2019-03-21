rm(list=ls())

###read in the arguments from bash script
args <- commandArgs(trailingOnly = F)

i <- args[length(args)]
i <- strsplit(i,"--")[[1]][2]
i <- as.numeric(i)

#set working directory
setwd("/QRISdata/Q0992")

#set seed
set.seed(12345)

#load packages
library(R.matlab)
library(plyr)
library(rstan)

#load data and transform into R data frame
tmp <- readMat("data/raw/datS3.mat")
tmp <- tmp[[1]]
dat <- adply(tmp,c(2,1))
names(dat) <- c("trial","subject","m_b","d_b","m_a","d_a","choose_a")

Nsubj = length(unique(dat$subject))

#prep data
data_list = list(Ntotal = dim(dat)[1],
                 Nsubj = Nsubj,
                 subj = as.numeric(dat$subject),
                 y=dat$choose_a,
                 d_a = dat$d_a,
                 d_b = dat$d_b,
                 m_a = dat$m_a,
                 m_b = dat$m_b)

models = c('hyperbolic',
           'exponential',
           'hyperbolic_gm',
           'prop_diff',
           'tradeoff',
           'ITCH',
           'const_sens',
           'mazur1987',
           'loewenstein1992',
           'mcclure2007',
           'killeen2009')


fit=stan(file=paste0("models/",models[i],".stan"),
         data=data_list,
         iter=2000,
         warmup=1000,
         cores=8,
         chains=8,
         control=list(max_treedepth=20,adapt_delta=0.99),
         seed=12345)

save(fit,file=paste0("data/derived/fits/fit_Bayes_",models[i],".RData"))
