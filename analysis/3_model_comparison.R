
library(rstan)
library(tidyverse)
library(R.matlab)
library(ggthemes)
library(loo)
library(xtable)


#Specify model details
source("models/model_details.R")

# #Extract random selection of iterations
# nsamples = 8000
# iter_samples = 1:8000
# 
# #Calculate Loo and WAIC for Experiment
# loo_list_e1 = list()
# waic_list_e1 = list()
# loo_list_e2 = list()
# waic_list_e2 = list()
# 
# for( exp in 1:2){
#   
#   ctr1=0
#   ctr2=0
#   
#   if(exp == 1){
#     #load data
#     dat=readRDS(file="data/raw/delay_data_lba_2019.rds") %>%
#       mutate(subject = as.numeric(as.factor(subject)))
#     
#     #format data for model fitting  
#     dat$m_a = dat$amount
#     dat$d_a = dat$delay
#     dat$m_b = 100
#     dat$d_b = 0    #note: b is smaller, sooner option
#     dat$choose_a = dat$Choice
#     
#   }
#   
#   if(exp == 2){
#     #load data and transform into R data frame (a = LL, b = SS)
#     tmp = readMat("data/raw/datS3.mat")
#     tmp = tmp[[1]]
#     dat = plyr::adply(tmp,c(2,1))
#     names(dat) = c("trial","subject","m_b","d_b","m_a","d_a","choose_a")
#   }
#   
#   Nsubj = length(unique(dat$subject))
#   
#   for (model in 1:length(models)){
#     
#     print(model)
#     
#     ctr2=ctr2+1
#     
#     #Preallocate matrix for storing log likelihoods
#     log_lik_mat = matrix(NA,nsamples,nrow(dat))
#     
#     #load likelihood
#     source(paste0("models/",names(models)[model],".R"))
#     
#     #load fit object
#     load(file=paste0("data/derived/fits99/fit_exp",exp,"_",names(models)[model],".RData"))
#     
#     #extract parameters
#     posts=rstan::extract(fit)
#     rm(fit)
#     
#     #parameters
#     parm_names = models[names(models)[model]][[1]]$parms
#     #parm_names
#     names(posts)
#     
#     for(sample in 1:nsamples){
#       
#       #index for inserting log_liks into the matrix
#       start = 1
#       #print(sample)
#       for(subj in 1:Nsubj){
#         
#         ctr1=ctr1+1
#         
#         #extract decisions for each subject
#         sim_dat = filter(dat,subject==subj)
#         
#         #transform (if necessary) and store sampled parameter
#         if(names(models)[model] == "hyperbolic"){
#           #no transforms necessary
#           parms = c(posts$k[iter_samples[sample],subj],
#                     posts$sigma[iter_samples[sample],subj])
#         }
#         
#         if(names(models)[model] == "exponential"){
#           #no transforms necessary
#           parms = c(posts$k[iter_samples[sample],subj],
#                     posts$sigma[iter_samples[sample],subj])
#         }
#         
#         if(names(models)[model] == "hyperbolic_gm"){
#           #no transforms necessary
#           parms = c(posts$k[iter_samples[sample],subj],
#                     posts$s[iter_samples[sample],subj],
#                     posts$sigma[iter_samples[sample],subj])
#         }
#         
#         if(names(models)[model] == "prop_diff"){
#           #delta needs to be unstandardised
#           parms = c(posts$delta_raw[iter_samples[sample],subj]*
#                       posts$delta_sd[iter_samples[sample]] + posts$delta_mean[iter_samples[sample]],
#                     posts$sigma[iter_samples[sample],subj])
#         }
#         
#         if(names(models)[model] == "tradeoff"){
#           #1 is added to theta
#           parms = c(posts$gamma[iter_samples[sample],subj],
#                     posts$tau[iter_samples[sample],subj],
#                     1+posts$theta_raw[iter_samples[sample],subj],
#                     posts$kappa[iter_samples[sample],subj],
#                     posts$alpha[iter_samples[sample],subj],
#                     posts$eps[iter_samples[sample],subj])
#         }
#         
#         if(names(models)[model] == "ITCH"){
#           #all parameters need to be unstandardised
#           parms = c(posts$B1_raw[iter_samples[sample],subj]*
#                       posts$B1_sd[iter_samples[sample]] + posts$B1_mean[iter_samples[sample]],
#                     posts$BxA_raw[iter_samples[sample],subj]*
#                       posts$BxA_sd[iter_samples[sample]] + posts$BxA_mean[iter_samples[sample]],
#                     posts$BxR_raw[iter_samples[sample],subj]*
#                       posts$BxR_sd[iter_samples[sample]] + posts$BxR_mean[iter_samples[sample]],
#                     posts$BtA_raw[iter_samples[sample],subj]*
#                       posts$BtA_sd[iter_samples[sample]] + posts$BtA_mean[iter_samples[sample]],
#                     posts$BtR_raw[iter_samples[sample],subj]*
#                       posts$BtR_sd[iter_samples[sample]] + posts$BtR_mean[iter_samples[sample]])
#         }
#         
#         if(names(models)[model] == "const_sens"){
#           #no transforms necessary
#           parms = c(posts$alpha[iter_samples[sample],subj],
#                     posts$beta[iter_samples[sample],subj],
#                     posts$sigma[iter_samples[sample],subj])
#         }
#         
#         if(names(models)[model] == "mazur1987"){
#           #no transforms necessary
#           parms = c(posts$k[iter_samples[sample],subj],
#                     posts$s[iter_samples[sample],subj],
#                     posts$sigma[iter_samples[sample],subj])
#         }
#         
#         if(names(models)[model] == "loewenstein1992"){
#           #beta_on_alpha needs to be multipled by alpha to calculate beta
#           parms = c(posts$alpha[iter_samples[sample],subj],
#                     posts$beta_on_alpha[iter_samples[sample],subj]*
#                       posts$alpha[iter_samples[sample],subj],
#                     posts$sigma[iter_samples[sample],subj])
#         }
#         
#         if(names(models)[model] == "mcclure2007"){
#           #no transform necessary
#           parms = c(posts$omega[iter_samples[sample],subj],
#                     posts$beta[iter_samples[sample],subj],
#                     posts$delta[iter_samples[sample],subj],
#                     posts$sigma[iter_samples[sample],subj])
#         }
#         
#         if(names(models)[model] == "killeen2009"){
#           #no transform necessary
#           parms = c(posts$alpha[iter_samples[sample],subj],
#                     posts$beta[iter_samples[sample],subj],
#                     posts$lambda[iter_samples[sample],subj],
#                     posts$sigma[iter_samples[sample],subj])
#         }
#         
#         #get choice probabilities under sampled parameter values
#         prob_a = likelihood(parms,sim_dat)
#         # sim_dat$choose_a_pp = as.numeric(runif(n=nrow(sim_dat)) < sim_dat$prob_a)
#         # sim_dat$sample = iter_samples[sample]
#         # sim_dat$model = names(models)[model] 
#         # sim_dat$exp = exp
#         
#         #specify end of range in log lik matrix
#         end = start + (nrow(sim_dat)-1)
#         
#         #fill log lik matrix for that participant
#         log_lik_mat[sample,start:end] = ifelse(sim_dat$choose_a==1,prob_a,1-prob_a)
#         
#         #set start of next participant range
#         start = end+1
#         
#         # if(exp==1){
#         #   pp_list_e1[[ctr1]] = sim_dat
#         # }
#         # if(exp==2){
#         #   pp_list_e2[[ctr1]] = sim_dat
#         # }
#         
#       } #end subject loop
#     } #end samples loop
#     
#     #save(log_lik_mat,file=paste0("data/derived/predictives/log_lik_",names(models)[model],"_exp_",exp,".RData"))
#     
#     #enter log_lik into ll list
#     if(exp==1){
#       log_lik_mat[log_lik_mat==0] <- 0.001
#       log_lik_mat[log_lik_mat==1] <- 0.999
#       loo_list_e1[[model]] = loo(log(log_lik_mat))
#       waic_list_e1[[model]] = waic(log(log_lik_mat))
#     }
#     if(exp==2){
#       log_lik_mat[log_lik_mat==0] <- 0.001
#       log_lik_mat[log_lik_mat==1] <- 0.999
#       loo_list_e2[[model]] = loo(log(log_lik_mat))
#       waic_list_e2[[model]] = waic(log(log_lik_mat))
#     }
#     
#   } #end model loop
# }
# 
# save(loo_list_e1,file="data/derived/loo_list_e1.RData")
# save(waic_list_e1,file="data/derived/waic_list_e1.RData")
# save(loo_list_e2,file="data/derived/loo_list_e2.RData")
# save(waic_list_e2,file="data/derived/waic_list_e2.RData")

load(file="data/derived/loo_list_e1.RData")
load(file="data/derived/waic_list_e1.RData")
load(file="data/derived/loo_list_e2.RData")
load(file="data/derived/waic_list_e2.RData")


names(loo_list_e1) = sapply(models,function(x) x[["label"]])   
names(waic_list_e1) = sapply(models,function(x) x[["label"]])  
names(loo_list_e2) = sapply(models,function(x) x[["label"]])  
names(waic_list_e2) = sapply(models,function(x) x[["label"]])  

loo_compare(loo_list_e1)
loo_compare(waic_list_e1)
loo_compare(loo_list_e2)
loo_compare(waic_list_e2)

x=compare(waic_list_e1[[1]],waic_list_e1[[2]])
#Create latex formatted table from data frame

names = c("Exponential","Hyperbolic","Hyperboloid",
          "Generalized Hyperbolic","Generalized Hyperbola",
          "Proportional Difference","Constant Sensitivity",
          "Additive Utility","Double Exponential","ITCH",
          "Tradeoff")

e1_waic = round(loo_compare(waic_list_e1),2)[names,c("waic","se_waic")]
e2_waic = round(loo_compare(waic_list_e2),2)[names,c("waic","se_waic")]

waic_table = cbind(e1_waic,rank(e1_waic[,"waic"]),
                 rep(" ",length(models)),
                 e2_waic,rank(e2_waic[,"waic"]))

colnames(waic_table) = c("WAIC","SE","Rank"," ","WAIC","SE","Rank")


latex_table=xtable(waic_table,
                   digits=c(0,2,2,0,0,2,2,0),
                   align=c("l","l",rep("r",ncol(waic_table)-1)),
                   caption="Results of Model Comparisons for Experiments 1 and 2",
                   label = "tab:waic_table")
# addtorow <- list()
# addtorow$pos <- list(nrow(waic_table))
# addtorow$command <- "\\hline  \\multicolumn{4}{p\\textwidth}{Note: 
# The `Participants Best Fit' column reports the number of participants for whom the relevant model yielded the lowest WAIC value.} \\\\ "
# 
print(latex_table,
      #add.to.row = addtorow,
      include.rownames=T,
      #include.colnames=F,
      hline.after=c(-1,0),
      caption.placement = "top")

#Need to add top row of table which indicates experiment manually in overleaf.