leow_likelihood=function(m_a,    #magnitude of option a
                               d_a,    #delay associated with option a
                               m_b,    #magnitude of option b
                               d_b,    #delay associated with option b
                               alpha,  #rate of departure from constant discounting, alpha > 0 
                               beta,   #rate of discounting, beta >0
                               sigma   #scale parameter
                              ){
  #Note: assumes linear value function
  u_a = m_a * (1+alpha*d_a)^(-beta/alpha) #utility of option a
  u_b = m_b * (1+alpha*d_b)^(-beta/alpha) #utility of option b
  p_a = pnorm((u_a-u_b)/sigma) #probability of selecting option a
  return(p_a)
}

leow_wrapper=function(pars,dat){
  #extract pars
  alpha=pars[1]
  beta=pars[2]
  sigma=pars[3]

  p_a = leow_likelihood(dat$m_a,dat$d_a,dat$m_b,dat$d_b,
                             alpha,beta,sigma)
  p_a = pmin(pmax(p_a,0.001),0.999)  
  neglnLs = -log(dat$choose_a*p_a + (1-dat$choose_a)*(1-p_a)) 
  return(sum(neglnLs))
}
