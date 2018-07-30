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

parameters {
  
  real<lower=0> omega_a;
  real<lower=0> omega_b;
   vector<lower=0,upper=1>[Nsubj] omega;

  real<lower=0> beta_mean;
  real<lower=0> beta_sd;
   vector<lower=0>[Nsubj] beta_raw;
  
  real<lower=0> delta_mean;
  real<lower=0> delta_sd;
   vector<lower=0>[Nsubj] delta_raw;
  
  real<lower=0> sigma_mean;
  real<lower=0> sigma_sd;
   vector<lower=0>[Nsubj] sigma_raw;
  
}

model {
  //model
  vector[Ntotal] p_a_logit;
  vector[Ntotal] beta = beta_raw[subj]*beta_sd + beta_mean;
  vector[Ntotal] delta = delta_raw[subj]*delta_sd + delta_mean;
  vector[Ntotal] sigma = sigma_raw[subj]*sigma_sd + sigma_mean;
  real u_a;
  real u_b;

  for(i in 1:Ntotal){
    u_a = m_a[i] * (omega[subj[i]] * exp(-beta[i]*d_a[i]) + (1-omega[subj[i]])*exp(-delta[i]*d_a[i]));
    u_b = m_b[i] * (omega[subj[i]] * exp(-beta[i]*d_b[i]) + (1-omega[subj[i]])*exp(-delta[i]*d_b[i]));
  
    p_a_logit[i] = (u_a - u_b) / sigma[i];
  }

  //priors
  omega_a ~ normal(0,1);
  omega_b ~ normal(0,1);
  
  beta_mean ~ normal(0,1);
  beta_sd ~ normal(0,1);
  
  delta_mean ~ normal(0,1);
  delta_sd ~ normal(0,1);
  
  sigma_mean ~ normal(0,1);
  sigma_sd ~ normal(0,1);
  
  omega ~ beta(omega_a,omega_b);
  beta_raw ~ normal(0,1);
  delta_raw ~ normal(0,1);
  sigma_raw ~ normal(0,1);
  
  //likelihood
  y~bernoulli_logit(p_a_logit);
}
