hyperbolic_gm_likelihood=function(m_a,    #magnitude of option a
                               d_a,    #delay associated with option a
                               m_b,    #magnitude of option b
                               d_b,    #delay associated with option b
                               k,      #discount factor
                               s,      #discount exponent
                               sigma   #scale parameter
                              ){
  u_a = m_a/(1+k*d_a^s) #utility of option a
  u_b = m_b/(1+k*d_b^s) #utility of option b
  p_a = pnorm((u_a-u_b)/sigma) #probability of selecting option a
  return(p_a)
}

hyperbolic_gm_wrapper=function(pars,dat){
  #extract pars
  k=pars[1]
  s=pars[2]
  sigma=pars[3]

  p_a = hyperbolic_gm_likelihood(dat$m_a,dat$d_a,dat$m_b,dat$d_b,
                             k,s,sigma)
  p_a = pmin(pmax(p_a,0.001),0.999)  
  neglnLs = -log(dat$choose_a*p_a + (1-dat$choose_a)*(1-p_a)) 
  return(sum(neglnLs))
}
