#Killeen, P. R. (2009). An additive-utility model of delay discounting. 
#Psychological Review, 116(3), 602-619.

killeen_likelihood=function(m_a,    #magnitude of option a
                            d_a,    #delay associated with option a
                            m_b,    #magnitude of option b
                            d_b,    #delay associated with option b
                            beta,   #shape of delay function (> 0)
                            lambda, #weight of delay (> 0)
                            sigma   #scale parameter
                            ){
  #Note: assumes value function v(m) = m^alpha is linear (alpha = 1)
  u_a = m_a - lambda*d_a^beta  #utility of option a
  u_b = m_b - lambda*d_b^beta  #utility of option b
  p_a = pnorm((u_a-u_b)/sigma) #probability of selecting option a
  return(p_a)
}

killeen_wrapper=function(pars,dat){
  #extract pars
  beta=pars[1]
  lambda=pars[2]
  sigma=pars[3]

  p_a = killeen_likelihood(dat$m_a,dat$d_a,dat$m_b,dat$d_b,
                             beta,lambda,sigma)
  p_a = pmin(pmax(p_a,0.001),0.999)  
  neglnLs = -log(dat$choose_a*p_a + (1-dat$choose_a)*(1-p_a)) 
  return(sum(neglnLs))
}
