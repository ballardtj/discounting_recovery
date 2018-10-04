rm(list=ls())

# ###read in the arguments listed at the command line
# args <- commandArgs(trailingOnly = F)
# 
# print(args)
# 
# i <- args[length(args)]
# i <- strsplit(i,"--")[[1]][2]
# i <- as.numeric(i)
# 
# print(i)
# 
# #set working directory
# setwd("~/discounting_recovery")

i = 5 #refit tradeoff model

#load packages
library(tidyverse)
library(rstan)

#load data
dat=readRDS(file="data/raw/data_delay_lba2018.rds")

#format data for model fitting
dat$m_a = dat$amount
dat$d_a = dat$delay
dat$m_b = 100
dat$d_b = 0    #note: b is smaller, sooner option
dat$choose_a = dat$Choice
Nsubj = length(unique(dat$subject))

#prep data
data_list = list(Ntotal = dim(dat)[1],
                 Nsubj = Nsubj,
                 subj = dat$subject,
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
         cores=4,
         chains=4,
         seed=12345)

save(fit,file=paste0("data/derived/fit_Bayes_",models[i],".RData"))
