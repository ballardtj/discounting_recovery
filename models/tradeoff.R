tradeoff_likelihood=function(m_a,   #magnitude of option a
                             d_a,    #delay associated with option a
                             m_b,    #magnitude of option b
                             d_b,    #delay associated with option b
                             gamma,  #diminishing absolute sensitivity to outcomes (linear at 0)
                             tau,    #diminishing absolute sensitivity to time
                             theta,  #superadditivity (1 = no supperadditivity)
                             kappa,  #time sensitivity (greater analogous to more discounting)
                             alpha,
                             eps     #noise ( -> 0 = no noise)
                              ){
  #subjective magnitudes
  sm_a = (1/gamma)*log1p(gamma*m_a)
  sm_b = (1/gamma)*log1p(gamma*m_b)
  Qm = sm_a - sm_b
  
  #subjective delays
  sd_a = (1/tau)*log1p(tau*d_a)
  sd_b = (1/tau)*log1p(tau*d_b)
  Qd = (kappa/alpha * log1p(alpha*((sd_a-sd_b)/theta)^theta))
  
  #choice probability
  p_a = Qm^(1/eps) / (Qm^(1/eps) + Qd^(1/eps))
  return(p_a)
}

tradeoff_wrapper=function(pars,dat){
  #extract pars
  gamma=pars[1]
  tau=pars[2]
  theta=pars[3]
  kappa=pars[4]
  alpha=pars[5]
  eps=pars[6]
  
  p_a = tradeoff_likelihood(dat$m_a,dat$d_a,dat$m_b,dat$d_b,
                             gamma,tau,theta,kappa,alpha,eps)
  p_a = pmin(pmax(p_a,0.001),0.999)  
  neglnLs = -log(dat$chooseLL*p_a + (1-dat$choose_a)*(1-p_a)) 
  return(sum(neglnLs))
}