likelihood=function(pars,dat){
    
  m_a = dat$m_a    #magnitude of option a
  d_a = dat$d_a    #delay associated with option a
  m_b = dat$m_b    #magnitude of option b
  d_b = dat$d_b    #delay associated with option b
  delta = pars[1]  #decision threshold
  sigma = pars[2]  #scale parameter
    
  #Equation 1
  alpha_m = pmax(abs(m_a),abs(m_b),1e-10)
  alpha_d = pmax(abs(d_a),abs(d_b),1e-10)
  
  pi_m = (pmax(abs(m_a),abs(m_b)) - pmin(abs(m_a),abs(m_b))) / alpha_m #max minus min of each pair divided by max of all options
  pi_d = (pmax(abs(d_a),abs(d_b)) - pmin(abs(d_a),abs(d_b))) / alpha_d
  
  d = pi_m-pi_d
  
  p_a = pnorm((d-delta)*sigma)    

  return(p_a)
}