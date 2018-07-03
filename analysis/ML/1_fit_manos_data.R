rm(list=ls())

#source libraries
library(tidyverse)
library(DEoptim)

#load data
dat=readRDS(file="../../data/data_delay_lba2018.rds")

#format data for model fitting
dat$m_a = dat$amount
dat$d_a = dat$delay
dat$m_b = 100
dat$d_b = 0
dat$choose_a = dat$Choice
dat$k = NA
dat$sigma = NA

#fit hyperbolic
source("../../models/hyperbolic.R")

Nsubj = length(unique(dat$subject))

for(i in 1:Nsubj){
  output=DEoptim(fn=hyperbolic_wrapper,
                 dat=dat[dat$subject==i,],
                 lower=c(0,0),
                 upper=c(10,100),
                 control=list(trace=0))

  dat$k[dat$subject==i] = output$optim$bestmem[1]
  dat$sigma[dat$subject==i] = output$optim$bestmem[2]
}

parms = dat %>% group_by(subject) %>%
  summarise(k = mean(k),
            sigma = mean(sigma))

ggplot(parms,aes(x=k,y=sigma)) +
  geom_density_2d()  +
  geom_point()

