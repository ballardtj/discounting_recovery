#Ebert, J., & Prelec, D. (2007). The Fragility of Time: Time-Insensitivity and Valuation of the near and Far Future. 
#Management Science, 53(9), 1423-1438. Retrieved from http://www.jstor.org/stable/20122300

#Note - they don't explicitly define a value function or utility function.

likelihood=function(pars,dat){
  
  m_a = dat$m_a    #magnitude of option a
  d_a = dat$d_a    #delay associated with option a
  m_b = dat$m_b    #magnitude of option b
  d_b = dat$d_b    #delay associated with option b
  alpha = pars[1]  #rate of departure from constant discounting, alpha > 0 
  beta = pars[2]   #rate of discounting, beta >0
  sigma = pars[3]  #scale parameter
  
  u_a = m_a * exp(-(alpha*d_a)^beta)   #utility of option a
  u_b = m_b * exp(-(alpha*d_b)^beta)   #utility of option b
  p_a = 1/(1+exp(- (u_a-u_b)*sigma)) #probability of selecting option a
  
  return(p_a)
}
