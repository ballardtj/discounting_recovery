
#Leowenstein, G. F. & Prelec, D (1992) Anomolies in intertemporal choice: evidence and interpretation. Quarterly Journal of Economics


  likelihood=function(pars,dat){
    
    m_a = dat$m_a    #magnitude of option a
    d_a = dat$d_a    #delay associated with option a
    m_b = dat$m_b    #magnitude of option b
    d_b = dat$d_b    #delay associated with option b
    alpha = pars[1]   #rate of departure from constant discounting, alpha > 0 
    beta = pars[2]    #rate of discounting, beta >0
    sigma = pars[3]  #scale parameter 
    
  #Note: assumes linear value function
  u_a = m_a * (1+alpha*d_a)^(-beta/alpha) #utility of option a
  u_b = m_b * (1+alpha*d_b)^(-beta/alpha) #utility of option b
  p_a = 1/(1+exp(- (u_a-u_b)*sigma)) #probability of selecting option a
  
  return(p_a)
}