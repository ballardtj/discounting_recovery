#TODO: For some reason, iteration 8 for the tradeoff model did not finish. So need to rerun that one.

rm(list=ls())

library(tidyverse)
library(rstan)
library(grid)
library(gridExtra)

#load model details
source("models/model_details.R")

#mcclure only returned samples for some models
mcclure = c(1:11,13,17,21,30:31,33:34,38,47:48,52:54,80:100)
models[[11]]$parms = c("beta","lambda","sigma") 

#samples used for correlation
corr_samples = sample(1:4000,100)

#loop through models
grob_list=list()
for(i in 1:11){
  print(i)
  
  #construct matrix for storing parameters
  #theta1.true theta2.true theta1.lower theta1.upper 
  n_parms = length(models[[i]]$parms)
  results = matrix(NA,100,n_parms*4 + 1)
  results[,1] = 1:100
  colnames(results) = c('number',paste0(models[[i]]$parms,".true"),
                        paste0(models[[i]]$parms,".lower"),
                        paste0(models[[i]]$parms,".median"),
                        paste0(models[[i]]$parms,".upper"))
  
  corr_list = list()
  ctr=0
  for(j in 1:100){
  
    print(j)
    #load fit object
    if(file.exists(paste0("data/derived/recovery/recovery_fit_",names(models)[i],"_",j,".RData"))) {
      
      if( i!=10 | j %in% mcclure ){
        
        load(paste0("data/derived/recovery/recovery_fit_",names(models)[i],"_",j,".RData"))
        samples = rstan::extract(fit)
        
        #loop through parameter values
        load(paste0("data/derived/recovery/recovery_generating_parms_",names(models)[i],"_",j,".RData"))
        results[j,paste0(models[[i]]$parms,".true")] = parms[j,]
        
        #loop through parameters
        for(k in models[[i]]$parms){
          ctr=ctr+1
          results[j,c(paste0(k,".lower"),paste0(k,".median"),paste0(k,".upper"))] = apply ( apply( samples[[k]] , 2, quantile , probs = c(0.025,0.5,0.975) ) , 1 , median) 
          
          #mean_each_iter = apply( samples[[k]] , 1, mean)
          #results[j,c(paste0(k,".lower"),paste0(k,".median"),paste0(k,".upper"))] = quantile(mean_each_iter,probs=c(0.025,0.5,0.975))  
          
          corr_list[[ctr]] = as_tibble(samples[[k]]) %>%
            dplyr::mutate(parm = k,
                   iter=1:n(),
                   true=unname(results[j,paste0(k,".true")] )) %>%
            filter(iter %in% corr_samples) %>%
            gather(key,value,V1:V100)
          
        }
      }
    }
  }
  
  plot_data = as.tibble(results) %>%
    #mutate(number = 1:n()) %>%
    gather(key="key",value="value",2:(1+n_parms*4)) %>%
    separate(col=key,into=c('parm','key')) %>%
    spread(key=key,value=value) #%>%
  # 
  # plot_data %>%
  #   group_by(parm) %>%
  #   
 
  corrs=bind_rows(corr_list) %>%
    group_by(iter,parm) %>%
    summarise(corr = cor(true,value)) %>%
    group_by(parm) %>%
    summarise(corr_lower = quantile(corr,0.025),
              corr_upper = quantile(corr,0.975))
            
  
  
  #loop through parameters
  
  
  
  
  # plot_list = list()
  #for(k in 1:length(models[[i]]$parms)){
  
  #tmp = filter(plot_data,parm==models[[i]]$parms[k])
  
  #replace "eps" with "epsilon" where relevant so greek character can be parsed
  plot_data$parm[plot_data$parm=="eps"] <- "epsilon"
  corrs$parm[corrs$parm=="eps"] <- "epsilon"
  
  
  plot_tmp = ggplot(plot_data,aes(x=true,y=median)) +
    geom_point() +
    geom_errorbar(aes(ymin=lower,ymax=upper),size=0.25) +
    geom_abline(linetype="dashed",size=1,col="blue") +
    #scale_x_continuous(name="Generating Value",limits=range(c(plot_data$lower,tmp$upper),na.rm=T),breaks= pretty(c(tmp$lower,tmp$upper),n=5)) +
    # scale_y_continuous(name="Recovered Value",limits=range(c(tmp$lower,tmp$upper),na.rm=T),breaks= pretty(c(tmp$lower,tmp$upper),n=5)) +
    # coord_fixed(ratio = 1) +
    facet_wrap(.~parm,scales = "free",labeller= label_parsed)+
    labs(x="Generating Value",y="Recovered\nValue",title=models[[i]]$label) +
    #annotate(geom = 'text', label = "italic(R) == 0.75",
    #         parse = TRUE, x = -Inf, y = Inf, hjust = 0, vjust = 1) +
    geom_text(data=corrs,mapping = aes(x = -Inf, y = Inf, label=paste0("R = [",round(corr_lower,2)," ",round(corr_upper,2),"]")),
              hjust   = -0.05,vjust   = 1.5,parse=FALSE,size=3) +
    theme(plot.margin = unit(c(2,2,2,2),"points"),
          plot.title = element_text(hjust = 0.5,face="bold"),
          axis.text = element_text(size=6))
  
  # if( tmp$parm[1] %in% c("sigma","delta","lambda","beta","gamma","theta","tau","omega","alpha","kappa","epsilon") ){
  #   
  #   plot_list[[k]] =plot_tmp + facet_wrap(~parm,labeller= label_parsed) 
  #   
  # } else {
  #   
  #   plot_list[[k]] = plot_tmp  + facet_wrap(~parm) + theme(strip.text=element_text(face="italic")) 
  #   
  # } 
  
  # }
  
  # grob_list[[i]] = arrangeGrob(grobs = plot_list,
  #                              ncol = 6,
  #                              nrow = 1,
  #                              respect = T,
  #                              left=textGrob("Recovered Value",gp=gpar(fontsize=8),rot=90),
  #                              bottom=textGrob("Generating Value",gp=gpar(fontsize=8)),
  #                              top = models[[i]]$label,gp=gpar(fontsize=10))
  
  grob_list[[i]] = plot_tmp
  
}

#manually arrange grid of plots

#TODO: Structure of figure

#TODO: Reduce spacing between panels?

#TODO: Names above groups of plots


# layout = rbind(c(1,1,1,1,1,1,NA,2,2,2,2,2,2,NA,3,3,3,3,3,3,3,3,3),
#                c(NA,4,4,4,4,4,4,4,4,4,NA,NA,NA,5,5,5,5,5,5,5,5,5,NA),
#                c(NA,6,6,6,6,6,6,NA,NA,NA,7,7,7,7,7,7,7,7,7,7,7,7,NA),
#                c(NA,8,8,8,8,8,8,8,8,8,NA,NA,NA,9,9,9,9,9,9,9,9,9,NA),
#                c(NA,NA,NA,NA,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,NA,NA,NA,NA),
#                c(NA,NA,NA,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,NA,NA,NA)
# )


layout = rbind(c(1,1,1,1,1,1,1,NA,8,8,8,8,8,8,8),
               c(2,2,2,2,2,2,2,NA,9,9,9,9,9,9,9),
               c(3,3,3,3,3,3,3,NA,9,9,9,9,9,9,9),
               c(4,4,4,4,4,4,4,NA,10,10,10,10,10,10,10),
               c(5,5,5,5,5,5,5,NA,10,10,10,10,10,10,10),
               c(6,6,6,6,6,6,6,NA,NA,NA,NA,NA,NA,NA,NA),
               c(7,7,7,7,7,7,7,NA,NA,NA,NA,NA,NA,NA,NA)
)


# fig = arrangeGrob(
#   grob_list[[2]],
#   grob_list[[1]],
#   grob_list[[3]]
# )

fig = arrangeGrob(
  grob_list[[2]],
  grob_list[[1]],
  grob_list[[3]],
  grob_list[[8]], #Gen Hyperbolic
  grob_list[[9]], #Gen Hyperbola
  grob_list[[4]], #Proportional Difference
  grob_list[[7]], #Constant Sensitivity
  grob_list[[11]], #Additive Utility
  grob_list[[10]], #Double Exponential
  #grob_list[[6]], #ITCH
  grob_list[[5]], #Tradeoff
  ncol=1,
  layout_matrix = layout#,
)

# png(filename = "presentations/jdm-2019/experiment1.png",width=5,height=10,units="in")
#   
# dev.off()

ggsave("experiment1.png",grid.arrange(fig),width=10,height=12)
ggsave("presentations/jdm-2019/experiment1.png",grid.arrange(fig),width=10,height=12)


# 
# 
# 
# 
# 
# 
# 
# layout_1 = t(matrix(c(1,1,1,1,1,1,NA,2,2,2,2,2,2,NA,3,3,3,3,3,3,3,3,3)))
# 
# 
# row_1 = arrangeGrob(
#   grob_list[[2]], #Exponential
#   grob_list[[1]], #Hyperbolic
#   grob_list[[3]], #Hyperboloid
#   layout_matrix = layout_1#,
#   #respect = T
#   #widths=1:5
# )
# 
# layout_2 = t(matrix(c(1,1,1,1,1,1,NA,2,2,2,2,2,2)))
# row_2 = arrangeGrob(
#   grob_list[[8]], #Gen Hyperbolic
#   grob_list[[9]], #Gen Hyperbola
#   layout_matrix = layout_2#,
# )
# 
# layout_3 = t(matrix(c(1,1,1,1,NA,2,2,2,2,2,2,2,2)))
# row_3 = arrangeGrob(
#   grob_list[[4]], #Proportional Difference 
#   grob_list[[10]], #Double Exponential
#   layout_matrix = layout_3
# )
# 
# layout_4 = t(matrix(c(1,1,1,1,1,1,NA,2,2,2,2,2,2)))
# row_4 = arrangeGrob(
#   grob_list[[7]], #Constant Sensitivity
#   grob_list[[11]], #Additive Utility
#   layout_matrix = layout_4
# )
# grid.arrange(row_2,row_3,row_4)
# 
# row_5 = arrangeGrob(grob_list[[6]])
# row_6 = arrangeGrob(grob_list[[5]],nrow=1,ncol=1)
# 
# grid.arrange(row_2,row_3,row_4,row_5,ncol=1)
# 
# x=arrangeGrob(
#   row_1,row_2,row_3,row_4,row_5,row_6,
#   nrow=6,ncol=1
# )
# 
# 
# 
# 
# # grid.arrange(x)
# # get_length = function(i){
# #   return(length(models[[i]]$parms))
# # }
# # 
# # recovery_fig = arrangeGrob(grobs = grob_list, nrow = length(grob_list),ncol=1)
# # 
# # grid.arrange(recovery_fig)
# 
# ggsave("figures/recovery_fig3.eps",plot=x,height=20,width=12)
# 
# 
# # recovery_fig = arrangeGrob(grobs = grob_list[1:2], nrow = length(grob_list[1:2]),ncol=1)
# # 
# # grid.arrange(recovery_fig)
