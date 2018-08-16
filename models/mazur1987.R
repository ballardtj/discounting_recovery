#Mazur 1987
likelihood=function(pars,dat){
    
    m_a = dat$m_a    #magnitude of option a
    d_a = dat$d_a    #delay associated with option a
    m_b = dat$m_b    #magnitude of option b
    d_b = dat$d_b    #delay associated with option b
    k = pars[1]  #discount factor, k > 0
    s = pars[2] #discount exponent, s > 0
    sigma = pars[3]  #scale parameter
    
  u_a = m_a/(1+k*d_a^s) #utility of option a
  u_b = m_b/(1+k*d_b^s) #utility of option b
  p_a = 1/(1+exp(- (u_a-u_b)*sigma)) #probability of selecting option a
  return(p_a)
}