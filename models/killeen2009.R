#Killeen, P. R. (2009). An additive-utility model of delay discounting. 
#Psychological Review, 116(3), 602-619.

likelihood=function(pars,dat){

  m_a = dat$m_a    #magnitude of option a
  d_a = dat$d_a    #delay associated with option a
  m_b = dat$m_b    #magnitude of option b
  d_b = dat$d_b    #delay associated with option b
  alpha = pars[1]  #curvature of value function (> 0)
  beta = pars[2]   #shape of delay function (> 0)
  lambda = pars[3] #weight of delay (> 0)
  sigma = pars[4]  #scale parameter                           
                              
  #Note: assumes value function v(m) = m^alpha is linear (alpha = 1)
  u_a = m_a^alpha - lambda*d_a^beta  #utility of option a
  u_b = m_b^alpha - lambda*d_b^beta  #utility of option b
  p_a = 1/(1+exp(- (u_a-u_b)*sigma)) #probability of selecting option a
  
  return(p_a)
}

