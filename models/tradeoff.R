likelihood=function(pars,dat){
    
  m_a = dat$m_a    #magnitude of option a
  d_a = dat$d_a    #delay associated with option a
  m_b = dat$m_b    #magnitude of option b
  d_b = dat$d_b    #delay associated with option b
  gamma = pars[1]  #diminishing absolute sensitivity to outcomes (linear at 0)
  tau = pars[2]    #diminishing absolute sensitivity to time
  theta = pars[3]  #superadditivity (1 = no supperadditivity)
  kappa = pars[4]  #time sensitivity (greater analogous to more discounting)
  alpha = pars[5]
  eps = pars[6]   #noise ( -> 0 = no noise)
  
  #subjective magnitudes
  sm_a = (1/gamma)*log1p(gamma*m_a)
  sm_b = (1/gamma)*log1p(gamma*m_b)
  Qm = sm_a - sm_b
  
  #subjective delays
  sd_a = (1/tau)*log1p(tau*d_a)
  sd_b = (1/tau)*log1p(tau*d_b)
  Qd = kappa/alpha * log1p( alpha*((sd_a-sd_b)/theta)^theta ) #log(1+x)
  
  #choice probability
  p_a = Qm^(1/eps) / (Qm^(1/eps) + Qd^(1/eps))
  return(p_a)
}