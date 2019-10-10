
#   Model description:
#   
#   * Hierarchical
# * Chains: 4
# * Burnin: 1000
# * Post-burnin: 1000
# * Thinning: 1
# * Total samples: 8000


library(rstan)
library(coda)
library(knitr)

# library(shinystan)
# 
# launch_shinystan(fit)

source("models/model_details.R")

stan2coda <- function(fit) {
  mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}

for( exp in 1){
for (i in 1:length(models)){

  if(i != 4 & i != 6){names(models)[i] = paste0(names(models)[i],"_gamma")}
  
  print(models[[i]]$name)
  
  load(file=paste0("data/derived/fits/fit_exp",exp,"_", names(models)[i],".RData"))
  
  smrystan = summary(fit)
                                                                                         
  print(kable(round(smrystan[[1]],3)),format="markdown")

  write.csv(smrystan[[1]], file=paste0("figures/summary_exp",exp,"_",names(models)[i],".csv"))
  
  #grep(models[[i]]$parms,rownames(smrystan[[1]]))
  
  smrystan[[1]][grep(models[[i]]$parms,rownames(smrystan[[1]])),]

  print(kable(effectiveSize(stan2coda(fit))[models[[i]]$hypers],format="markdown",col.names="n_eff"))

  try(print(kable(gelman.diag(stan2coda(fit))[[1]][models[[i]]$hypers,],format="markdown")))

  traceplot = rstan::traceplot(fit,par=models[[i]]$hypers)
  ggsave(file=paste0("figures/convergence_trace_exp",exp,"_",names(models)[i],".pdf"),plot=traceplot)


  pdf(file=paste0("figures/convergence_pairs_exp",exp,"_",names(models)[i],".pdf"),width=10,height=10)
     pairs(fit,pars=models[[i]]$hypers)
  dev.off()
  
}
}

#For experiment 1, all models except double exponential converge with 10000 burnin, 5000 samples, and delta of 0.8

