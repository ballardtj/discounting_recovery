
library(rstan)
library(tidyverse)
library(R.matlab)
library(ggthemes)
library(loo)
library(xtable)

#Specify model details
source("models/model_details.R")

#Extract random selection of iterations
set.seed(1234)
nsamples = 100  #number of different parameter values tested
iter_samples = sample(1:8000,size=nsamples,replace=TRUE) #vector of sampled iterations

#Generate predictives for each sample for each participant, model, and experiment
pp_list_e1=list()
pp_list_e2=list()

for( exp in 1:2){
  
  ctr1=0
  ctr2=0
  
  if(exp == 1){
    #load data
    dat=readRDS(file="data/raw/delay_data_lba_2019.rds") %>%
      mutate(subject = as.numeric(as.factor(subject)))
    
    #format data for model fitting  
    dat$m_a = dat$amount
    dat$d_a = dat$delay
    dat$m_b = 100
    dat$d_b = 0    #note: b is smaller, sooner option
    dat$choose_a = dat$Choice
    
  }
  
  if(exp == 2){
    #load data and transform into R data frame (a = LL, b = SS)
    tmp = readMat("data/raw/datS3.mat")
    tmp = tmp[[1]]
    dat = plyr::adply(tmp,c(2,1))
    names(dat) = c("trial","subject","m_b","d_b","m_a","d_a","choose_a")
  }
  
  Nsubj = length(unique(dat$subject))
  
  for (model in 1:length(models)){
    
    print(model)
    
    ctr2=ctr2+1
    
    #load likelihood
    source(paste0("models/",names(models)[model],".R"))
    
    #load fit object
    load(file=paste0("data/derived/fits99/fit_exp",exp,"_",names(models)[model],".RData"))
    
    #extract parameters
    posts=rstan::extract(fit)
    rm(fit)
    
    #parameters
    parm_names = models[names(models)[model]][[1]]$parms
    parm_names
    names(posts)
    
    for(sample in 1:nsamples){
      
      #print(sample)
      for(subj in 1:Nsubj){
        
        ctr1=ctr1+1
        
        #extract decisions for each subject
        sim_dat = filter(dat,subject==subj)
        
        #transform (if necessary) and store sampled parameter
        if(names(models)[model] == "hyperbolic"){
          #no transforms necessary
          parms = c(posts$k[iter_samples[sample],subj],
                    posts$sigma[iter_samples[sample],subj])
        }
        
        if(names(models)[model] == "exponential"){
          #no transforms necessary
          parms = c(posts$k[iter_samples[sample],subj],
                    posts$sigma[iter_samples[sample],subj])
        }
        
        if(names(models)[model] == "hyperbolic_gm"){
          #no transforms necessary
          parms = c(posts$k[iter_samples[sample],subj],
                    posts$s[iter_samples[sample],subj],
                    posts$sigma[iter_samples[sample],subj])
        }
        
        if(names(models)[model] == "prop_diff"){
          #delta needs to be unstandardised
          parms = c(posts$delta_raw[iter_samples[sample],subj]*
                      posts$delta_sd[iter_samples[sample]] + posts$delta_mean[iter_samples[sample]],
                    posts$sigma[iter_samples[sample],subj])
        }
        
        if(names(models)[model] == "tradeoff"){
          #1 is added to theta
          parms = c(posts$gamma[iter_samples[sample],subj],
                    posts$tau[iter_samples[sample],subj],
                    1+posts$theta_raw[iter_samples[sample],subj],
                    posts$kappa[iter_samples[sample],subj],
                    posts$alpha[iter_samples[sample],subj],
                    posts$eps[iter_samples[sample],subj])
        }
        
        if(names(models)[model] == "ITCH"){
          #all parameters need to be unstandardised
          parms = c(posts$B1_raw[iter_samples[sample],subj]*
                      posts$B1_sd[iter_samples[sample]] + posts$B1_mean[iter_samples[sample]],
                    posts$BxA_raw[iter_samples[sample],subj]*
                      posts$BxA_sd[iter_samples[sample]] + posts$BxA_mean[iter_samples[sample]],
                    posts$BxR_raw[iter_samples[sample],subj]*
                      posts$BxR_sd[iter_samples[sample]] + posts$BxR_mean[iter_samples[sample]],
                    posts$BtA_raw[iter_samples[sample],subj]*
                      posts$BtA_sd[iter_samples[sample]] + posts$BtA_mean[iter_samples[sample]],
                    posts$BtR_raw[iter_samples[sample],subj]*
                      posts$BtR_sd[iter_samples[sample]] + posts$BtR_mean[iter_samples[sample]])
        }
        
        if(names(models)[model] == "const_sens"){
          #no transforms necessary
          parms = c(posts$alpha[iter_samples[sample],subj],
                    posts$beta[iter_samples[sample],subj],
                    posts$sigma[iter_samples[sample],subj])
        }
        
        if(names(models)[model] == "mazur1987"){
          #no transforms necessary
          parms = c(posts$k[iter_samples[sample],subj],
                    posts$s[iter_samples[sample],subj],
                    posts$sigma[iter_samples[sample],subj])
        }
        
        if(names(models)[model] == "loewenstein1992"){
          #beta_on_alpha needs to be multipled by alpha to calculate beta
          parms = c(posts$alpha[iter_samples[sample],subj],
                    posts$beta_on_alpha[iter_samples[sample],subj]*
                      posts$alpha[iter_samples[sample],subj],
                    posts$sigma[iter_samples[sample],subj])
        }
        
        if(names(models)[model] == "mcclure2007"){
          #no transform necessary
          parms = c(posts$omega[iter_samples[sample],subj],
                    posts$beta[iter_samples[sample],subj],
                    posts$delta[iter_samples[sample],subj],
                    posts$sigma[iter_samples[sample],subj])
        }
        
        if(names(models)[model] == "killeen2009"){
          #no transform necessary
          parms = c(posts$alpha[iter_samples[sample],subj],
                    posts$beta[iter_samples[sample],subj],
                    posts$lambda[iter_samples[sample],subj],
                    posts$sigma[iter_samples[sample],subj])
        }
        
        #get choice probabilities under sampled parameter values
        sim_dat$prob_a = likelihood(parms,sim_dat)
        sim_dat$choose_a_pp = as.numeric(runif(n=nrow(sim_dat)) < sim_dat$prob_a)
        sim_dat$sample = iter_samples[sample]
        sim_dat$model = names(models)[model] 
        sim_dat$exp = exp
        
        if(exp==1){
          pp_list_e1[[ctr1]] = sim_dat
        }
        if(exp==2){
          pp_list_e2[[ctr1]] = sim_dat
        }
        
      } #end subject loop
    } #end samples loop
  } #end model loop
}

save(pp_list_e1,file="data/derived/predictives/pp_list_e1.RData")
save(pp_list_e1,file="data/derived/predictives/pp_list_e2.RData")

#Generate plot for Experiment 1

pp_data_e1 = bind_rows(pp_list_e1) %>%
  mutate(amount = factor(m_a,levels=seq(120,500,by=20),labels=paste0("$",sort(unique(amount)))),
         delay = d_a,
         model = factor(model,levels=names(models),labels=c('Hyperbolic','Exponential','Hyperboloid',
                                                             'Proportional Difference','Tradeoff','ITCH',
                                                             'Constant Sensitivity','Generalized Hyperbolic',
                                                             'Generalized Hyperbola','Double Exponential',
                                                             'Additive Utility')),
         model = factor(model,levels=c("Exponential","Hyperbolic","Hyperboloid",
                                       "Generalized Hyperbolic","Generalized Hyperbola",
                                       "Proportional Difference","Constant Sensitivity",
                                       "Additive Utility","Double Exponential","ITCH",
                                       "Tradeoff")))

data = pp_data_e1 %>%
  filter(model == "Exponential")  %>%
  filter(sample == min(sample)) %>%
  group_by(amount,delay) %>%
  summarise(prob_m = mean(choose_a),
            prob_u = NA,
            prob_l = NA)

sim = pp_data_e1 %>%
  group_by(amount,delay,model,sample) %>%
  summarise(choose_a = mean(choose_a_pp)) %>%
  group_by(amount,delay,model) %>%
  summarise(prob_m = mean(choose_a),
            prob_u = quantile(choose_a,0.975), 
            prob_l = quantile(choose_a,0.025))

# library(RColorBrewer)

bcols = brewer.pal(6,"Dark2")
bcols = colorblind_pal()(6)

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
  facet_wrap(~amount) + labs(x="Larger-Later Delay",y="Proportion Choosing Larger-Later Option",colour="Model") +
  theme(legend.position = "bottom") +
  scale_colour_manual(name="Model",values = m_colours) +
  scale_linetype_manual(name="Model",values = m_linetypes)

pp_fig

ggsave(file="figures/predictives_e1.pdf",plot=pp_fig,width=9,height=7)


#Generate plot for Experiment 2

pp_data_e2 = bind_rows(pp_list_e2) %>%
  mutate(trial = as.numeric(trial),
         effect = case_when(
              trial < 101 ~ "Delay Duration Effect",
              trial >= 101 & trial < 201 ~ "Common Difference Effect",
              trial >= 201 ~ "Magnitude Effect"),
         group = case_when(
           effect != "Magnitude Effect" ~ d_a,
           effect == "Magnitude Effect" ~ m_a
         ),
         model = factor(model,levels=names(models),labels=c('Hyperbolic','Exponential','Hyperboloid',
                                                            'Proportional Difference','Tradeoff','ITCH',
                                                            'Constant Sensitivity','Generalized Hyperbolic',
                                                            'Generalized Hyperbola','Double Exponential',
                                                            'Additive Utility')),
         model = factor(model,levels=c("Exponential","Hyperbolic","Hyperboloid",
                                       "Generalized Hyperbolic","Generalized Hyperbola",
                                       "Proportional Difference","Constant Sensitivity",
                                       "Additive Utility","Double Exponential","ITCH",
                                       "Tradeoff")))

#delay duration - shorter delay - 2-40
#magnitude - small amount 2 - 40

data = pp_data_e2 %>%
  filter(model == "Exponential")  %>%
  filter(sample == min(sample)) %>%
  group_by(effect,group) %>%
  summarise(prob_m = mean(choose_a),
            prob_u = NA,
            prob_l = NA)

sim = pp_data_e2 %>%
  group_by(effect,group,model,sample) %>%
  summarise(choose_a = mean(choose_a_pp)) %>%
  group_by(effect,group,model) %>%
  summarise(prob_m = mean(choose_a),
            prob_u = quantile(choose_a,0.975), 
            prob_l = quantile(choose_a,0.025))

pp_fig_dd = ggplot(filter(sim,effect=="Delay Duration Effect"),aes(x=group)) +
  #geom_ribbon(aes(ymin=prob_l,ymax=prob_u),colour="skyblue") +
  geom_line(aes(y=prob_m,colour=model,linetype=model)) +
  geom_point(data=filter(data,effect=="Delay Duration Effect"),aes(y=prob_m),colour="black",alpha=0.5) +
  # geom_errorbar(data=data,aes(ymin=prob_l,ymax=prob_u),colour="red") +
  labs(x="Larger-Later Delay (Days)",y="Proportion Choosing Larger-Later Option",colour="Model") +
  scale_colour_manual(name="Model",values = m_colours) +
  scale_linetype_manual(name="Model",values = m_linetypes) +
  coord_cartesian(ylim=c(0,1)) +
  #scale_x_continuous(breaks=seq(2,40,2)) +
  ggtitle("Delay Duration Effect") +
  theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5))

pp_fig_cd = ggplot(filter(sim,effect=="Common Difference Effect"),aes(x=group)) +
  #geom_ribbon(aes(ymin=prob_l,ymax=prob_u),colour="skyblue") +
  geom_line(aes(y=prob_m,colour=model,linetype=model)) +
  geom_point(data=filter(data,effect=="Common Difference Effect"),aes(y=prob_m),colour="black",alpha=0.5) +
  # geom_errorbar(data=data,aes(ymin=prob_l,ymax=prob_u),colour="red") +
  labs(x="Larger-Later Delay (Days)",y="Proportion Choosing Larger-Later Option",colour="Model") +
  theme(legend.position = "bottom") +
  scale_colour_manual(name="Model",values = m_colours) +
  scale_linetype_manual(name="Model",values = m_linetypes) +
  coord_cartesian(ylim=c(0,1)) +
  #scale_x_continuous(breaks=seq(2,40,2)) +
  ggtitle("Common Difference Effect") +
  theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(colour="white"))

pp_fig_m = ggplot(filter(sim,effect=="Magnitude Effect"),aes(x=group)) +
  #geom_ribbon(aes(ymin=prob_l,ymax=prob_u),colour="skyblue") +
  geom_line(aes(y=prob_m,colour=model,linetype=model)) +
  geom_point(data=filter(data,effect=="Magnitude Effect"),aes(y=prob_m),colour="black",alpha=0.5) +
  # geom_errorbar(data=data,aes(ymin=prob_l,ymax=prob_u),colour="red") +
  labs(x="Larger-Later Amount ($)",y="Proportion Choosing Larger-Later Option",colour="Model") +
  theme(legend.position = "bottom") +
  scale_colour_manual(name="Model",values = m_colours) +
  scale_linetype_manual(name="Model",values = m_linetypes) +
  coord_cartesian(ylim=c(0,1)) +
  #scale_x_continuous(breaks=seq(2,40,2)) +
  ggtitle("Magnitude Effect") +
  theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(colour="white"))


library(gridExtra)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

legend = g_legend(pp_fig_dd)

pp_fig = grid.arrange(
    arrangeGrob(  pp_fig_dd + theme(legend.position="none"),
                  pp_fig_cd + theme(legend.position="none"),
                  pp_fig_m + theme(legend.position="none"),
                  nrow = 1),
    legend,
    nrow=2,
    heights=c(5,1)
)

ggsave(file="figures/predictives_e2.pdf",plot=pp_fig,width=9,height=5)

