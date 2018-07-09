tradeoff_likelihood=function(m_a,   #magnitude of option a
                             d_a,    #delay associated with option a
                             m_b,    #magnitude of option b
                             d_b,    #delay associated with option b
                             delta,  #decision threshold
                             sigma   #scale parameter
                              ){
  #subjective values
  sm_a = (1/gamma)*log(1+gamma*m_a)
  sm_b = (1/gamma)*log(1+gamma*m_b)
  sm_diff = sm_a - sm_b
  
  sd_a = (1/tau)*log(1+tau+d_a)
  sd_b = (1/tau)*log(1+tau+d_b)
  sd_diff = sd_a - sd_b
  
  (kappa/alpha * log(1+alpha*((sd_a-sd_b)/theta)^theta))
  
  alpha_m = pmax(abs(m_a),abs(m_b),1e-10)
  alpha_d = pmax(abs(d_a),abs(d_b),1e-10)
  
  pi_m = (pmax(abs(m_a),abs(m_b)) - pmin(abs(m_a),abs(m_b))) / alpha_m) #max minus min of each pair divided by max of all options
  pi_d = (pmax(abs(d_a),abs(d_b)) - pmin(abs(d_b),abs(d_b))) / alpha_d)
  
  d = pi_m-pi_d
  pLL = pnorm((d-delta) / sigma)
  return(pLL)
}

tradeoff_wrapper=function(pars,dat){
  #extract pars
  gamma=pars[1]
  eps=pars[2]
  tau=pars[3]
  theta=pars[4]
  kappa=pars[5]
  alpha=pars[6]
  
  pLLs = tradeoff_likelihood(dat$m_a,dat$d_a,dat$m_b,dat$d_b,
                             gamma,eps,tau,theta,kappa,alpha)
  pLLs = pmin(pmax(pLLs,0.001),0.999)  
  neglnLs = -log(dat$chooseLL*pLLs + (1-dat$chooseLL)*(1-pLLs)) 
  return(sum(neglnLs))
}


tradeoff.likelihood=function(ll_m,ll_d,ss_m,ss_d,
                             gamma,eps,tau,theta,kappa,alpha){
  #Equation 5
  term1 = ( (1/gamma)*log(1+gamma*ll_m) - (1/gamma)*log(1+gamma*ss_m) )^(1/eps) 
  term2 = ( ( (1/tau)*log(1+tau+ll_d) - (1/tau)*log(1+tau*ss_d) )/theta )^theta
  term3 = (kappa/alpha)*log(1+alpha*term2)^(1/eps)
  pLL = term1/(term1+term3)
  
  #This happens at some parameter values and causes NaN
  pLL[term1==0 & term3==0] <- 0.5
  
  return(pLL)
}

