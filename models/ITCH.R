likelihood=function(pars,dat){
  
  m_a = dat$m_a    #magnitude of option a
  d_a = dat$d_a    #delay associated with option a
  m_b = dat$m_b    #magnitude of option b
  d_b = dat$d_b    #delay associated with option b
  B1 = pars[1]   #intercept 
  BxA = pars[2]  #weight given to absolute reward
  BxR = pars[3]  #weight given to relative reward
  BtA = pars[4]  #weight given to absolute delay
  BtR = pars[5]  #weight given to relative reward

  #option _a is later, _b is sooner
  mean_x = max(((m_a+m_b)/2),1e-10)
  mean_t = max(((d_a+d_b)/2),1e-10)
  
  weighted_xA = BxA*(m_a - m_b)
  weighted_xR = BxR*((m_a - m_b)/ mean_x )
  weighted_tA = BtA*(d_a - d_b)
  weighted_tR = BtR*((d_a - d_b)/ mean_t )
  
  logit = B1 + weighted_xA + weighted_xR + weighted_tA + weighted_tR
  p_a = 1/(1+exp(-logit))
  return(p_a)
}
