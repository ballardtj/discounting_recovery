#Ebert, J., & Prelec, D. (2007). The Fragility of Time: Time-Insensitivity and Valuation of the near and Far Future. 
#Management Science, 53(9), 1423-1438. Retrieved from http://www.jstor.org/stable/20122300

#Note - they don't explicitly define a value function or utility function.

const_sens_likelihood=function(m_a,    #magnitude of option a
                               d_a,    #delay associated with option a
                               m_b,    #magnitude of option b
                               d_b,    #delay associated with option b
                               alpha,  #rate of departure from constant discounting, alpha > 0 
                               beta,   #rate of discounting, beta >0
                               sigma   #scale parameter
                              ){
  u_a = m_a * exp(-(alpha*d_a)^beta)   #utility of option a
  u_b = m_b * exp(-(alpha*d_b)^beta)   #utility of option b
  p_a = pnorm((u_a-u_b)/sigma) #probability of selecting option a
  return(p_a)
}

const_sens_wrapper=function(pars,dat){
  #extract pars
  alpha=pars[1]
  beta=pars[2]
  sigma=pars[3]

  p_a = const_sens_likelihood(dat$m_a,dat$d_a,dat$m_b,dat$d_b,
                        alpha,beta,sigma)
  p_a = pmin(pmax(p_a,0.001),0.999)  
  neglnLs = -log(dat$choose_a*p_a + (1-dat$choose_a)*(1-p_a)) 
  return(sum(neglnLs))
}
