#Time Discounting for Primary Rewards
#Samuel M. McClure, Keith M. Ericson, David I. Laibson, George Loewenstein, Jonathan D. Cohen
#Journal of Neuroscience 23 May 2007, 27 (21) 5796-5804; DOI: 10.1523/JNEUROSCI.4246-06.2007

double_exp_likelihood=function(m_a,    #magnitude of option a
                               d_a,    #delay associated with option a
                               m_b,    #magnitude of option b
                               d_b,    #delay associated with option b
                               omega,  #relative influence of beta vs delta system, omega is between 0 and 1 
                               beta,   #discount factor associated with beta system, beta >0
                               delta,  #discount factor associated with delta system, delta > 0
                               sigma   #scale parameter
                              ){
  u_a = m_a * (omega*exp(-beta*d_a) + (1-omega)*exp(-delta*d_a))    #utility of option a
  u_b = m_b * (omega*exp(-beta*d_b) + (1-omega)*exp(-delta*d_b))    #utility of option b
  p_a = pnorm((u_a-u_b)/sigma) #probability of selecting option a
  return(p_a)
}

double_exp_wrapper=function(pars,dat){
  #extract pars
  omega=pars[1]
  beta=pars[2]
  delta=pars[3]
  sigma=pars[4]

  p_a = double_exp_likelihood(dat$m_a,dat$d_a,dat$m_b,dat$d_b,
                        omega,beta,delta,sigma)
  p_a = pmin(pmax(p_a,0.001),0.999)  
  neglnLs = -log(dat$choose_a*p_a + (1-dat$choose_a)*(1-p_a)) 
  return(sum(neglnLs))
}
