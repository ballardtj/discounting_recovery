data {
  int<lower=0> Ntotal;
  int<lower=0> Nsubj;
  int<lower=0> subj[Ntotal];
  int y[Ntotal];
  vector[Ntotal] m_a;
  vector[Ntotal] m_b;
  vector[Ntotal] d_a;
  vector[Ntotal] d_b;
}  

transformed data{
  real alpha_m; 
  real alpha_d;
  real pi_m;
  real pi_d;
  vector[Ntotal] d;
  
  for(i in 1:Ntotal){
    alpha_m = fmax( fmax( fabs(m_a[i]),fabs(m_b[i])), 1e-10);
    alpha_d = fmax( fmax( fabs(d_a[i]),fabs(d_b[i])), 1e-10);
    
    pi_m = (fmax( fabs(m_a[i]),fabs(m_b[i])) - fmin( fabs(m_a[i]),fabs(m_b[i])) ) / alpha_m;
    pi_d = (fmax( fabs(d_a[i]),fabs(d_b[i])) - fmin( fabs(d_a[i]),fabs(d_b[i])) ) / alpha_d;
    
    d[i] = pi_m - pi_d;
  }
}
 
parameters {
  
  real delta_mean;
  real<lower=0> delta_sd;
  vector[Nsubj] delta_raw;

  real<lower=0> sigma_mean;
  real<lower=0> sigma_sd;
  vector<lower=0>[Nsubj] sigma_raw;

}

model {
  //model
   vector[Ntotal] delta = delta_raw[subj]*delta_sd + delta_mean;
   vector[Ntotal] sigma = sigma_raw[subj]*sigma_sd + sigma_mean;

   vector[Ntotal] p_a;
   
   for(i in 1:Ntotal){
     p_a[i] = Phi_approx((d[i]-delta[subj[i]]) * sigma[subj[i]]);
    }
   
  //priors
  delta_mean ~ normal(0,1);
  delta_sd ~ normal(0,1);
  
  sigma_mean ~ normal(0,1);
  sigma_sd ~ normal(0,1);
  
  sigma_raw ~ normal(0,1);
  delta_raw ~ normal(0,1);

  //likelihood
  y~bernoulli(p_a);
}
