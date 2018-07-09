ITCH_likelihood=function(m_a,    #magnitude of option a
                               d_a,    #delay associated with option a
                               m_b,    #magnitude of option b
                               d_b,    #delay associated with option b
                               B1,      #intercept term
                               BxA,     #weight given to absolute reward
                               BxR,     #weight given to relative reward
                               BtA,     #weight given to absolute delay
                               BtR      #weight given to relative reward
                              ){
  #option _a is later, _b is sooner
  mean_x = max(((m_a+m_b)/2),1e-10)
  mean_t = max(((d_a+d_b)/2),1e-10)
  
  weighted_xA = Bxa*(m_a - m_b)
  weighted_xR = BxR*((m_a - m_b)/ mean_x )
  weighted_tA = Bta*(d_a - d_b)
  weighted_tR = BtR*((d_a - d_b)/ mean_t )
  
  logit = B1 + weighted_xA + weighted_xR + weighted_tA + weighted_tR
  p_a = 1/(1+exp(-logit))
  return(p_a)
}

ITCH_wrapper=function(pars,dat){
  #extract pars
  B1=pars[1]
  BxA=pars[2]
  BxR=pars[3]
  BtA=pars[4]
  BtR=pars[5]

  p_a = ITCH_likelihood(dat$m_a,dat$d_a,dat$m_b,dat$d_b,
                             B1,BxA,BxR,BtA,BtR)
  p_a = pmin(pmax(p_a,0.001),0.999)  
  neglnLs = -log(dat$choose_a*p_a + (1-dat$choose_a)*(1-p_a)) 
  return(sum(neglnLs))
}
