
library(rstan)
library(tidyverse)

dat=readRDS(file="data/raw/data_delay_lba2018.rds")

#format data for model fitting
dat$m_a = dat$amount
dat$d_a = dat$delay
dat$m_b = 100
dat$d_b = 0    #note: b is smaller, sooner option
dat$choose_a = dat$Choice
Nsubj = length(unique(dat$subject))

#Specify model details

models = list(hyperbolic = list(
  name = 'Hyperbolic model', 
  hypers=c('k_mean','k_sd','sigma_mean','sigma_sd')),
  
  exponential = list(
    name = 'Exponential model',
    hypers=c('k_mean','k_sd','sigma_mean','sigma_sd')),
  
  hyperbolic_gm = list( 
    name = 'Hyperbolic model (Green & Myerson version)',
    hypers=c('k_mean','k_sd','s_mean','s_sd','sigma_mean','sigma_sd')),
  
  prop_diff = list( 
    name = 'Proportional difference model',
    hypers=c('delta_mean','delta_sd','sigma_mean','sigma_sd')),  
  
  tradeoff = list(
    name = 'Tradeoff model',
    hypers=c("gamma_mean","gamma_sd","tau_mean","tau_sd","theta_mean","theta_sd",
             "kappa_mean","kappa_sd","alpha_mean","alpha_sd","eps_mean","eps_sd")),
  
  ITCH = list(
    name = 'ITCH model',
    hypers=c("B1_mean","B1_sd","BxA_mean","BxA_sd","BxR_mean","BxR_sd",
             "BtA_mean","BtA_sd","BtR_mean","BtR_sd")),
  
  const_sens = list(
    name = 'Constant Sensitivity Model (Ebert & Prelec, 2007)',
    hypers=c('alpha_mean','alpha_sd','beta_mean','beta_sd','sigma_mean','sigma_sd')),
  
  mazur1987 = list(
    name = 'Hyperboloid Model (Mazur, 1987)',
    hypers=c('k_mean','k_sd','s_mean','s_sd','sigma_mean','sigma_sd')),
  
  loewenstein1992 = list(
    name = "Loewenstein & Prelec's (1992) Model",
    hypers=c('alpha_mean','alpha_sd','beta_mean','beta_sd','sigma_mean','sigma_sd')),
  
  mcclure2007 = list(
    name = "Double Exponential Model (McClure et al., 2007)",
    hypers=c('omega_a','omega_b','beta_mean','beta_sd',"delta_mean","delta_sd","sigma_mean",'sigma_sd')),  
  
  killeen2009 = list(
    name = "Additive Utility Model (Killeen, 2009)",
    hypers=c('beta_mean','beta_sd',"lambda_mean","lambda_sd","sigma_mean",'sigma_sd'))
)


pp_list=list()
ctr=0
for (i in 1:length(models)){
 
  print(models[[i]]$name)
  #cat('\n')
  #load fit object
  load(file=paste0("data/derived/fits/fit_Bayes_", names(models)[i],".RData"))
  
  #load likelihood
  source(paste0("models/",names(models)[i],".R"))
  
  Nsamp=100
  posts=rstan::extract(fit)
  Npost=dim(posts[[1]])[1]
  samples=sample(x=1:Npost,size=Nsamp)
  
  for(j in 1:Nsamp){
    
    #set up matrix of recentered parameters
    if(names(models)[i] %in% c('hyperbolic','exponential')){
      par_mat = cbind(
        posts$k_raw[samples[i],]*posts$k_sd[samples[i]] + posts$k_mean[samples[i]],
        posts$sigma_raw[samples[i],]*posts$sigma_sd[samples[i]] + posts$sigma_mean[samples[i]]
      )
    }
    if(names(models)[i] %in% c('hyperbolic_gm','mazur1987')){
      par_mat = cbind(
        posts$k_raw[samples[i],]*posts$k_sd[samples[i]] + posts$k_mean[samples[i]],
        posts$s_raw[samples[i],]*posts$s_sd[samples[i]] + posts$s_mean[samples[i]],
        posts$sigma_raw[samples[i],]*posts$sigma_sd[samples[i]] + posts$sigma_mean[samples[i]]
      )
    }
    if(names(models)[i] %in% c('prop_diff')){
      par_mat = cbind(
        posts$delta_raw[samples[i],]*posts$delta_sd[samples[i]] + posts$delta_mean[samples[i]],
        posts$sigma_raw[samples[i],]*posts$sigma_sd[samples[i]] + posts$sigma_mean[samples[i]]
      )
    }
    if(names(models)[i] %in% c('tradeoff')){
      par_mat = cbind(
        posts$gamma_raw[samples[i],]*posts$gamma_sd[samples[i]] + posts$gamma_mean[samples[i]],
        posts$tau_raw[samples[i],]*posts$tau_sd[samples[i]] + posts$tau_mean[samples[i]],
        posts$theta_raw[samples[i],]*posts$theta_sd[samples[i]] + posts$theta_mean[samples[i]],
        posts$kappa_raw[samples[i],]*posts$kappa_sd[samples[i]] + posts$kappa_mean[samples[i]],
        posts$alpha_raw[samples[i],]*posts$alpha_sd[samples[i]] + posts$alpha_mean[samples[i]],
        posts$eps_raw[samples[i],]*posts$eps_sd[samples[i]] + posts$eps_mean[samples[i]]
      )
    }
    if(names(models)[i] %in% c('ITCH')){
      par_mat = cbind(
        posts$B1_raw[samples[i],]*posts$B1_sd[samples[i]] + posts$B1_mean[samples[i]],
        posts$BxA_raw[samples[i],]*posts$BxA_sd[samples[i]] + posts$BxA_mean[samples[i]],
        posts$BxR_raw[samples[i],]*posts$BxR_sd[samples[i]] + posts$BxR_mean[samples[i]],
        posts$BtA_raw[samples[i],]*posts$BtA_sd[samples[i]] + posts$BtA_mean[samples[i]],
        posts$BtR_raw[samples[i],]*posts$BtR_sd[samples[i]] + posts$BtR_mean[samples[i]]
      )
    }
    if(names(models)[i] %in% c('const_sens','loewenstein1992')){
      par_mat = cbind(
        posts$alpha_raw[samples[i],]*posts$alpha_sd[samples[i]] + posts$alpha_mean[samples[i]],
        posts$beta_raw[samples[i],]*posts$beta_sd[samples[i]] + posts$beta_mean[samples[i]],
        posts$sigma_raw[samples[i],]*posts$sigma_sd[samples[i]] + posts$sigma_mean[samples[i]]
      )
    }
    if(names(models)[i] %in% c('mcclure2007')){
      par_mat = cbind(
        posts$omega[samples[i],],
        posts$beta_raw[samples[i],]*posts$beta_sd[samples[i]] + posts$beta_mean[samples[i]],
        posts$delta_raw[samples[i],]*posts$delta_sd[samples[i]] + posts$delta_mean[samples[i]],
        posts$sigma_raw[samples[i],]*posts$sigma_sd[samples[i]] + posts$sigma_mean[samples[i]]
      )
    }
    if(names(models)[i] %in% c('killeen2009')){
      par_mat = cbind(
        posts$beta_raw[samples[i],]*posts$beta_sd[samples[i]] + posts$beta_mean[samples[i]],
        posts$lambda_raw[samples[i],]*posts$lambda_sd[samples[i]] + posts$lambda_mean[samples[i]],
        posts$sigma_raw[samples[i],]*posts$sigma_sd[samples[i]] + posts$sigma_mean[samples[i]]
      )
    }
    
    for(k in 1:Nsubj){
      ctr=ctr+1
      pp_dat_tmp = dat[dat$subject==k,]
      pp_dat_tmp$p_a= likelihood(par_mat[k,],pp_dat_tmp)
      pp_dat_tmp$sample = samples[j]
      pp_dat_tmp$model = names(models)[i] 
      pp_list[[ctr]] = pp_dat_tmp
    }
  }
}


pp_data=bind_rows(pp_list) %>%
  mutate(amount = factor(amount,levels=seq(120,500,by=20),labels=paste0("$",sort(unique(amount)))),
         model = factor(model,levels=names(models),labels=c('Hyperbolic','Exponential','Hyperboloid',
                                                             'Proportional Difference','Tradeoff','ITCH',
                                                             'Constant Sensitivity','Generalized Hyperbolic',
                                                             'Generalized Hyperbola','Double Exponential',
                                                             'Additive Utility')))
                                                      
  
  data = pp_data %>%
    filter(model == "Exponential")  %>%
    filter(sample == min(sample)) %>%
    group_by(amount,delay) %>%
    summarise(prob_m = mean(choose_a),
              prob_u = prob_m + sd(choose_a)/sqrt(n()),
              prob_l = prob_m - sd(choose_a)/sqrt(n()))
  
  sim = pp_data %>%
    group_by(amount,delay,model,sample) %>%
    summarise(p_a = mean(p_a)) %>%
    group_by(amount,delay,model) %>%
    summarise(prob_m = mean(p_a),
              prob_u = quantile(p_a,0.975), 
              prob_l = quantile(p_a,0.025))
  
  library(RColorBrewer)
  
  display.brewer.pal(6,"Dark2")
  
  bcols = brewer.pal(6,"Dark2")
  
  m_colours = c(1:6,1:5)
  
  m_colours = c(bcols,bcols[1:5])
  
  #m_colours = c('#e41a1c',"#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33",
  #              '#e41a1c',"#377eb8","#4daf4a","#984ea3","#ff7f00")
  m_linetypes = c("solid","solid","solid","solid","solid","solid",
                  "dashed","dashed","dashed","dashed","dashed")
  
#https://community.rstudio.com/t/ggplot2-change-legend-title-while-controlling-line-types-and-colors/14966/2  
  
 pp_fig = ggplot(sim,aes(x=delay)) +
      #geom_ribbon(aes(ymin=prob_l,ymax=prob_u),colour="skyblue") +
      geom_line(aes(y=prob_m,colour=model,linetype=model)) +
      geom_point(data=data,aes(y=prob_m),colour="black",alpha=0.5) +
     # geom_errorbar(data=data,aes(ymin=prob_l,ymax=prob_u),colour="red") +
      facet_wrap(~amount) + labs(x="Delay",y="Probability of Choosing Later/Larger Option",colour="Model") +
   theme(legend.position = "bottom") +
   scale_colour_manual(name="Model",values = m_colours) +
   scale_linetype_manual(name="Model",values = m_linetypes)
  
pp_fig

ggsave(file="figures/predictives.pdf",plot=pp_fig,width=9,height=7)
