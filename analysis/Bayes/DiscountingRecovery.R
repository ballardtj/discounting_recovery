rm(list=ls())
library(rstan)
library(dplyr)
library(shinystan)

#Simulate Model

d=expand.grid(R=1,k=seq(0,5,by=.1),D1=0:10,D2=0:10)
d$V1=d$R/(1+d$k*d$D1)
d$V2=d$R/(1+d$k*d$D2)
d$m=d$V1-d$V2
d$m2=d$V1/d$V2

ggplot(data=d,aes(x=D1,y=V1)) + geom_line()

ggplot(data=d,aes(x=D1,y=D2)) + geom_raster(aes(fill=m)) + facet_wrap(~k)


#Recovery at 0 should be pretty poor because data looks the same as if k is very high (e.g., no difference in options)


df=d%>%filter(D2==7,D1==1)
ggplot(data=df,aes(x=k,y=m)) + geom_line() + ylab("Value Preference for Sooner Option")


df=d%>%filter(D2==7,D1==1)
ggplot(data=df,aes(x=k,y=m2)) + geom_line() + ylab("Value Preference for Sooner Option")

#df=d%>%filter(D2==5) 
#ggplot(data=df,aes(x=k,y=m,group=D1,colour=D1)) + geom_line()




#simulate data - Kirby paradigm 
#    - fine when you assume deterministic decisions 
#    - but recovery goes to shit when you introduce noise
#    - worth seeing whether the fit with noisey data can be attenutated with more observations or a different noise model (e.g., a cumulative distribution functions)
sir=c(34,54,78,28,47,80,22,54,67,25,49,69,19,40,55,24,34,54,14,27,41,15,25,33,11,20,31)
ldr=c(35,55,80,30,50,85,25,60,75,30,60,85,25,55,75,35,50,80,25,50,75,35,60,80,30,55,85)
delay=c(186,117,162,179,160,157,136,111,119,80,89,91,53,62,61,29,30,30,19,21,20,13,14,14,7,7,7)

r=1
ks = c(0, 0.01, 0.05, .1, .5, 1, 5)
alphas=c(.01, .05, .1, .2, .4)
nsamples=
samples.list=list()
ctr=0
for(k in ks){
  for(alpha in alphas){
    ctr=ctr+1
    d=data.frame(r_a=rep(ldr,r),r_b=rep(sir,r),d_a=rep(delay,r),d_b=0,k=k,alpha=alpha)
    d<-d %>% mutate(v_a = r_a/(1+k*d_a),
                v_b = r_b/(1+k*d_b),
                v_diff =v_a-v_b,
                theta = (v_diff>0)*(1-alpha) + (v_diff<0)*alpha + (v_diff==0)*0.5,
                y=(runif((dim(d)[1]))<theta)*1 )

#prep data
data.list = list(Ntotal = dim(d)[1],
                 y=d$y,
                 d_a = d$d_a,
                 d_b = d$d_b,
                 r_a = d$r_a,
                 r_b = d$r_b)

fit=stan(file="discount.stan",
         data=data.list,
         iter=10000,
         warmup=2000,
         cores=4,
         chains=4)




#launch_shinystan(fit)

samples.list[[ctr]]=data.frame(k=extract(fit,"k")$k,
                   k.true=k,
                   alpha.true=alpha)

  }
}

samples.frame=bind_rows(samples.list)                   
 
#Histograms of samples (kind of hard to interpret)                  
ggplot(data=samples.frame) + geom_histogram(aes(x=k),bins=100) + #geom_vline(xintercept = k.true[1],colour='red') +
  facet_grid(k.true~alpha.true)

samples.ci=samples.frame %>%
  group_by(k.true,alpha.true) %>%
  summarise(k.lower = quantile(k,probs=.025),
            k.mid = quantile(k,probs=.5),
            k.upper = quantile(k,probs=.975),
            k.error=mean(k.true-k.mid))

#Histograms of samples (kind of hard to interpret)                  
ggplot(data=samples.ci) + 
  geom_errorbar(aes(x=as.factor(k.true),ymin=k.lower,ymax=k.upper)) +
  geom_errorbar(aes(x=as.factor(k.true),ymin=k.true,ymax=k.true),color="red") +
  geom_point(aes(x=as.factor(k.true),y=k.mid)) + 
  facet_wrap(~alpha.true)
  
  
#summary(fit, pars = c("k", "alpha"))$summary
#traceplot(fit,pars = c("k","alpha"))

#hist(extract(fit,"k")$k)

#plot(fit, show_density = TRUE, ci_level = 0.95, fill_color = "lightblue",pars="k")
#plot(fit, plotfun = "hist", pars = "k", include = T,bins=100)



#TODO - #Plot posteriors as a function of k and noise for different models.
#k = 0, 0.01, 0.05, .1, .5, 1, 5,
#alpha = .01, .05, .1, .2, .4

#try different paradigm (e.g., Junyi, 2017) where delays are compared to one another.


sr = c(10,10,10,15,15,20,40,40,40,45,45,50,10,10,25,25)
sd = c(1,1,1,2,2,3,7,7,7,8,8,9,1,1,4,4)
lr =  c(15,20,25,20,25,25,45,50,55,50,55,55,40,55,40,55)
ld = c(2,3,4,3,4,4,8,9,10,9,10,10,7,10,7,10)


r=1
ks = c(0, 0.01, 0.05, .1, .5, 1, 5)
alphas=c(.01, .05, .1, .2, .4)
nsamples=
  samples.list=list()
ctr=0
for(k in ks){
  for(alpha in alphas){
    ctr=ctr+1
    d=data.frame(r_a=rep(lr,r),r_b=rep(sr,r),d_a=rep(ld,r),d_b=rep(sd,r),k=k,alpha=alpha)
    d<-d %>% mutate(v_a = r_a/(1+k*d_a),
                    v_b = r_b/(1+k*d_b),
                    v_diff =v_a-v_b,
                    theta = (v_diff>0)*(1-alpha) + (v_diff<0)*alpha + (v_diff==0)*0.5,
                    y=(runif((dim(d)[1]))<theta)*1 )
    
    #prep data
    data.list = list(Ntotal = dim(d)[1],
                     y=d$y,
                     d_a = d$d_a,
                     d_b = d$d_b,
                     r_a = d$r_a,
                     r_b = d$r_b)
    
    fit=stan(file="discount.stan",
             data=data.list,
             iter=10000,
             warmup=2000,
             cores=4,
             chains=4)
    
    
    
    
    #launch_shinystan(fit)
    
    samples.list[[ctr]]=data.frame(k=extract(fit,"k")$k,
                                   k.true=k,
                                   alpha.true=alpha)
    
  }
}

samples.frame=bind_rows(samples.list)                   

#Histograms of samples (kind of hard to interpret)                  
ggplot(data=samples.frame) + geom_histogram(aes(x=k),bins=1000) + #geom_vline(xintercept = k.true[1],colour='red') +
  facet_grid(k.true~alpha.true)

samples.ci=samples.frame %>%
  group_by(k.true,alpha.true) %>%
  summarise(k.lower = quantile(k,probs=.025),
            k.mid = quantile(k,probs=.5),
            k.upper = quantile(k,probs=.975),
            k.error=mean(k.true-k.mid))

#Histograms of samples (kind of hard to interpret)                  
ggplot(data=samples.ci) + 
  geom_errorbar(aes(x=as.factor(k.true),ymin=k.lower,ymax=k.upper)) +
  geom_errorbar(aes(x=as.factor(k.true),ymin=k.true,ymax=k.true),color="red") +
  geom_point(aes(x=as.factor(k.true),y=k.mid)) + 
  facet_wrap(~alpha.true)



