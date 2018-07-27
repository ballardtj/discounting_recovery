prop_diff_likelihood=function(m_a,   #magnitude of option a
                             d_a,    #delay associated with option a
                             m_b,    #magnitude of option b
                             d_b,    #delay associated with option b
                             delta,  #decision threshold
                             sigma   #scale parameter
                              ){
  #Equation 1
  alpha_m = pmax(abs(m_a),abs(m_b),1e-10)
  alpha_d = pmax(abs(d_a),abs(d_b),1e-10)
  
  pi_m = (pmax(abs(m_a),abs(m_b)) - pmin(abs(m_a),abs(m_b))) / alpha_m #max minus min of each pair divided by max of all options
  pi_d = (pmax(abs(d_a),abs(d_b)) - pmin(abs(d_a),abs(d_b))) / alpha_d
  
  d = pi_m-pi_d
  p_a = pnorm((d-delta) / sigma)
  return(p_a)
}

prop_diff_wrapper=function(pars,dat){
  #extract pars
  delta=pars[1]
  sigma=pars[2]
  
  p_a = prop_diff_likelihood(dat$m_a,dat$d_a,dat$m_b,dat$d_b,
                               delta,sigma)
  p_a = pmin(pmax(p_a,0.001),0.999)  
  neglnLs = -log(dat$chooseLL*p_a + (1-dat$chooseLL)*(1-p_a)) 
  return(sum(neglnLs))
}