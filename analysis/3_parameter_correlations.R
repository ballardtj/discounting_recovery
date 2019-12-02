
rm(list=ls())

library(rstan)
#library(plyr)
library(tidyverse)
library(R.matlab)
library(reshape2)
library(ggthemes)
library(loo)
library(xtable)

#Specify model details
source("models/model_details.R")

#Extract random selection of iterations
set.seed(71342)
nsamples = 100  #number of different parameter values tested
iter_samples = sample(1:nsamples,size=nsamples,replace=TRUE) #vector of sampled iterations

for (model in 1:length(models)){
    
  #load data and transform into R data frame (a = LL, b = SS)
  tmp = readMat("data/raw/datS3.mat")
  tmp = tmp[[1]]
  dat = plyr::adply(tmp,c(2,1))
  names(dat) = c("trial","subject","m_b","d_b","m_a","d_a","choose_a")
  Nsubj = length(unique(dat$subject))
  
  parms_frame = expand.grid(subject = 1:Nsubj,
                            sample = 1:nsamples)
  
  for(p in 1:length(models[[model]]$parms)){
    parms_frame[models[[model]]$parms[p]] = NA
  }
  
  #List to store parameter estimates for each experiment
  parm_list_e345 = list()
  
  print(model)

  for( exp in 3:5){
    
    #load fit object
    load(file=paste0("data/derived/fits99/fit_exp",exp,"_",names(models)[model],".RData"))
    
    #extract parameters
    posts=rstan::extract(fit)
    rm(fit)
    
    #parameters
    parm_names = models[names(models)[model]][[1]]$parms
    parm_names
    names(posts)
    
    #for(sample in 1:nsamples){
      
      #for(subj in 1:Nsubj){
        
        #transform (if necessary) and store sampled parameter
        if(names(models)[model] == "hyperbolic"){
          #no transforms necessary
          parms = cbind( melt(posts$k), 
                         melt(posts$sigma)$value ) 
          names(parms) = c("sample","subject","k","sigma")
        }
        
        if(names(models)[model] == "exponential"){
          #no transforms necessary
          parms = cbind( melt(posts$k), 
                         melt(posts$sigma)$value ) 
          names(parms) = c("sample","subject","k","sigma")
        }
        
        if(names(models)[model] == "hyperbolic_gm"){
          #no transforms necessary
          parms = cbind( melt(posts$k), 
                         melt(posts$s)$value, 
                         melt(posts$sigma)$value ) 
          names(parms) = c("sample","subject","k","s","sigma")
        }
        
        if(names(models)[model] == "prop_diff"){
          #delta needs to be unstandardised
          delta_tmp = sweep(posts$delta_raw,1,posts$delta_sd,`*`)
          delta_tmp2 = sweep(delta_tmp,1,posts$delta_mean,`+`)  
          parms = cbind(melt(delta_tmp2),
                        melt(posts$sigma)$value)
          names(parms) = c("sample","subject","delta","sigma")
        }
        
        if(names(models)[model] == "tradeoff"){
          #1 is added to theta
          parms = cbind( melt(posts$gamma), 
                         melt(posts$tau)$value, 
                         melt(posts$theta_raw)$value+1,
                         melt(posts$kappa)$value, 
                         melt(posts$alpha)$value, 
                         melt(posts$eps)$value)
          names(parms) = c("sample","subject","gamma","tau","theta","kappa","alpha","eps")
        }
        
        if(names(models)[model] == "ITCH"){
          
          #all parameters need to be unstandardised
          B1_tmp = sweep(posts$B1_raw,1,posts$B1_sd,`*`)
          B1_tmp2 = sweep(B1_tmp,1,posts$B1_mean,`+`)  
          BxA_tmp = sweep(posts$BxA_raw,1,posts$BxA_sd,`*`)
          BxA_tmp2 = sweep(BxA_tmp,1,posts$BxA_mean,`+`)  
          BxR_tmp = sweep(posts$BxR_raw,1,posts$BxR_sd,`*`)
          BxR_tmp2 = sweep(BxR_tmp,1,posts$BxR_mean,`+`) 
          BtA_tmp = sweep(posts$BtA_raw,1,posts$BtA_sd,`*`)
          BtA_tmp2 = sweep(BtA_tmp,1,posts$BtA_mean,`+`) 
          BtR_tmp = sweep(posts$BtR_raw,1,posts$BtR_sd,`*`)
          BtR_tmp2 = sweep(BtR_tmp,1,posts$BtR_mean,`+`) 
          
          parms = cbind(melt(B1_tmp2),
                        melt(BxA_tmp2)$value,
                        melt(BxR_tmp2)$value,
                        melt(BtA_tmp2)$value,
                        melt(BtR_tmp2)$value)
          
          names(parms) = c("sample","subject","B1","BxA","BxR","BtA","BtR")
          
        }
        
        if(names(models)[model] == "const_sens"){
          #no transforms necessary
          parms = cbind( melt(posts$alpha), 
                         melt(posts$beta)$value, 
                         melt(posts$sigma)$value ) 
          names(parms) = c("sample","subject","alpha","beta","sigma")
        }
        
        if(names(models)[model] == "mazur1987"){
          #no transforms necessary
          parms = cbind( melt(posts$k), 
                         melt(posts$s)$value, 
                         melt(posts$sigma)$value ) 
          names(parms) = c("sample","subject","k","s","sigma")
        }
        
        if(names(models)[model] == "loewenstein1992"){
          #beta_on_alpha needs to be multipled by alpha to calculate beta
          parms = cbind( melt(posts$alpha), 
                         melt(posts$beta_on_alpha)$value * melt(posts$alpha)$value , 
                         melt(posts$sigma)$value ) 
          names(parms) = c("sample","subject","alpha","beta","sigma")
        }
        
        if(names(models)[model] == "mcclure2007"){
          #no transform necessary
          parms = cbind( melt(posts$omega), 
                         melt(posts$beta)$value, 
                         melt(posts$delta)$value,
                         melt(posts$sigma)$value) 
          names(parms) = c("sample","subject","omega","beta","delta","sigma")
        }
        
        if(names(models)[model] == "killeen2009"){
          #no transform necessary
          parms = cbind( melt(posts$alpha), 
                         melt(posts$beta)$value, 
                         melt(posts$lambda)$value,
                         melt(posts$sigma)$value) 
          names(parms) = c("sample","subject","alpha","beta","lambda","sigma")
        }
        
        parms$exp=exp
        parm_list_e345[[exp-2]] = parms 
        
      #} #end subject loop
    #} #end samples loop
  } #end experiment loop
  save(parm_list_e345,file=paste0("data/derived/predictives/parm_consistency_",names(models)[model],".RData"))
} #end model loop


#extract correlations and put them in table
corr_list = list()
for(model in 1:length(models)){
  
  load(paste0("data/derived/predictives/parm_consistency_",names(models)[model],".RData")) 
  
  corrs_tmp = bind_rows(parm_list_e345) %>%
      mutate(exp = factor(exp,3:5,labels=c('e3','e4','e5'))) %>%
      gather(key=parm,value=value,models[[model]]$parms) %>%
      spread(key=exp,value=value) %>%
      group_by(sample,parm) %>%
      summarise(corr34 = cor(e3,e4),
             corr45 = cor(e4,e5),
             corr35 = cor(e3,e5)) %>%
      gather(key=corr,value=value,corr34:corr35) %>%
      mutate(corr = factor(corr,levels=c("corr34","corr35","corr45"),
                           labels=c("DD/CD","DD/M","CD/M"))) %>%
      group_by(parm,corr) %>%
      summarise(lower = quantile(value,0.025),
                upper = quantile(value,0.975)) %>%
      mutate(model = models[[model]]$label)
  
    corr_list[[model]] = corrs_tmp
  
}

order = c("Exponential","Hyperbolic","Hyperboloid",
          "Generalized Hyperbolic","Generalized Hyperbola",
          "Proportional Difference","Constant Sensitivity",
          "Additive Utility","Double Exponential","ITCH",
          "Tradeoff")

table = bind_rows(corr_list) %>%
  ungroup() %>%
  mutate(model = factor(model,levels=order),
         lower = sprintf("%.2f", round(lower,2)),
         upper = sprintf("%.2f", round(upper,2))) %>%
  unite(col=ci,lower,upper,sep=" ") %>%
  spread(key=corr,value=ci) %>%
  arrange(model,parm) %>% 
  select(model,everything())

colnames(table) = c("Model","Parameter","DD / CD","DD / M","CD / M")

#replace parameter names with latex math format
table$Parameter[table$Parameter=="k"] <- "$k$"
table$Parameter[table$Parameter=="s"] <- "$s$"
table$Parameter[table$Parameter=="sigma"] <- "$\\sigma$"
table$Parameter[table$Parameter=="alpha"] <- "$\\alpha$"
table$Parameter[table$Parameter=="beta"] <- "$\\beta$"
table$Parameter[table$Parameter=="delta"] <- "$\\delta$"
table$Parameter[table$Parameter=="lambda"] <- "$\\lambda$"
table$Parameter[table$Parameter=="theta"] <- "$\\theta$"
table$Parameter[table$Parameter=="gamma"] <- "$\\gamma$"
table$Parameter[table$Parameter=="tau"] <- "$\\tau$"
table$Parameter[table$Parameter=="eps"] <- "$\\epsilon$"
table$Parameter[table$Parameter=="kappa"] <- "$\\kappa$"
table$Parameter[table$Parameter=="omega"] <- "$\\omega$"
table$Parameter[table$Parameter=="B1"] <- "$\\beta_1$"
table$Parameter[table$Parameter=="BxA"] <- "$\\beta_{xA}$"
table$Parameter[table$Parameter=="BxR"] <- "$\\beta_{xR}$"
table$Parameter[table$Parameter=="BtA"] <- "$\\beta_{tA}$"
table$Parameter[table$Parameter=="BtR"] <- "$\\beta_{tR}$"

latex_table=xtable(table,
                   digits=c(0,0,0,0,0,0),
                   align=c("l","l","l",rep("r",ncol(table)-2)),
                   caption="Correlations between Parameters Estimated from Three Subsets of the Dai et al. (2014) Data.",
                   label = "tab:corr_table")
# addtorow <- list()
# addtorow$pos <- list(nrow(waic_table))
# addtorow$command <- "\\hline  \\multicolumn{4}{p\\textwidth}{Note: 
# The `Participants Best Fit' column reports the number of participants for whom the relevant model yielded the lowest WAIC value.} \\\\ "
# 
print(latex_table,
      #add.to.row = addtorow,
      include.rownames=F,
      #include.colnames=F,
      hline.after=c(-1,0,2,4,7,10,13,15,18,22,26,31,37),
      caption.placement = "top",
      sanitize.text.function=function(x){x})

#parse parameter symbol






# #correlations
# x=bind_rows(parm_list_e345) %>%
#   mutate(exp = factor(exp,3:5,labels=c('e3','e4','e5'))) %>%
#   gather(key=parm,value=value,k:sigma) %>%
#   spread(key=exp,value=value) %>%
#   group_by(sample,parm) %>%
#   summarise(corr34 = cor(e3,e4),
#          corr45 = cor(e4,e5),
#          corr35 = cor(e3,e5)) %>%
#   gather(key=corr,value=value,corr34:corr35) %>%
#   group_by(parm,corr) %>%
#   summarise(lower = quantile(value,0.025),
#             upper = quantile(value,0.975))


bind_rows(parm_list_e345) %>%
  mutate(exp = factor(exp,3:5,labels=c('e3','e4','e5'))) %>%
  gather(key=parm,value=value,models[[model]]$parms) %>%
  spread(key=exp,value=value) %>%
  group_by(subject,parm) %>%
  summarise(median_e3 = median(e3),
            lower_e3 = quantile(e3,0.025),
            upper_e3 = quantile(e3,0.975),
            median_e4 = median(e4),
            lower_e4 = quantile(e4,0.025),
            upper_e4 = quantile(e4,0.975),
            median_e5 = median(e5),
            lower_e5 = quantile(e5,0.025),
            upper_e5 = quantile(e5,0.975)) %>%
  ggplot(aes(x=median_e3,y=median_e5)) +
  geom_point() +
  geom_errorbar(aes(ymin=lower_e5,ymax=upper_e5),alpha=0.2) +
  geom_errorbarh(aes(xmin=lower_e3,xmax=upper_e3),alpha=0.2) +
  facet_wrap(~parm,scale="free")
