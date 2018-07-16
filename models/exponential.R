exponential_likelihood=function(m_a,    #magnitude of option a
                               d_a,    #delay associated with option a
                               m_b,    #magnitude of option b
                               d_b,    #delay associated with option b
                               k,      #discount factor
                               sigma   #scale parameter
                              ){
  u_a = m_a*exp(-d_a*k) #utility of option a
  u_b = m_b*exp(-d_b*k) #utility of option b
  p_a = pnorm((u_a-u_b)/sigma) #probability of selecting option a
  return(p_a)
}

exponential_wrapper=function(pars,dat){
  #extract pars
  k=pars[1]
  sigma=pars[2]

  p_a = exponential_likelihood(dat$m_a,dat$d_a,dat$m_b,dat$d_b,
                             k,sigma)
  p_a = pmin(pmax(p_a,0.001),0.999)  
  neglnLs = -log(dat$choose_a*p_a + (1-dat$choose_a)*(1-p_a)) 
  return(sum(neglnLs))
}
