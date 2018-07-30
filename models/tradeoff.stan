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
  
  real<lower=0> gamma_mean;
  real<lower=0> gamma_sd;
   vector<lower=0>[Nsubj] gamma_raw;

  real<lower=0> tau_mean;
  real<lower=0> tau_sd;
   vector<lower=0>[Nsubj] tau_raw;
  
  real<lower=0> theta_mean;
  real<lower=0> theta_sd;
   vector<lower=0>[Nsubj] theta_raw;
  
  real<lower=0> kappa_mean;
  real<lower=0> kappa_sd;
   vector<lower=0>[Nsubj] kappa_raw;
  
  real<lower=0> alpha_mean;
  real<lower=0> alpha_sd;
   vector<lower=0>[Nsubj] alpha_raw;
  
  real<lower=0> eps_mean;
  real<lower=0> eps_sd;
  vector<lower=0>[Nsubj] eps_raw;

}

model {
  //model
  vector[Ntotal] p_a;
  real sm_a;
  real sm_b;
  real Qm;
  real sd_a;
  real sd_b;
  real Qd;
  
  vector[Ntotal] gamma = gamma_raw[subj]*gamma_sd + gamma_mean;
  vector[Ntotal] tau = tau_raw[subj]*tau_sd + tau_mean;
  vector[Ntotal] theta = theta_raw[subj]*theta_sd + theta_mean;
  vector[Ntotal] kappa = kappa_raw[subj]*kappa_sd + kappa_mean;
  vector[Ntotal] alpha = alpha_raw[subj]*alpha_sd + alpha_mean;
  vector[Ntotal] eps = eps_raw[subj]*eps_sd + eps_mean;
  
  for(i in 1:Ntotal){
    sm_a = (1/gamma[i]) * log1p(gamma[i]*m_a[i]);
    sm_b = (1/gamma[i]) * log1p(gamma[i]*m_b[i]);
    Qm = sm_a - sm_b;
    
    sd_a = (1/tau[i]) * log1p(tau[i]*d_a[i]);
    sd_b = (1/tau[i]) * log1p(tau[i]*d_b[i]);
    Qd = (kappa[i]/alpha[i]) * log1p(alpha[i]*((sd_a-sd_b)/theta[i])^theta[i]);
    
    p_a[i] = pow(Qm,1/eps[i]) / (pow(Qm,1/eps[i]) + pow(Qd,1/eps[i]));
  }

  
  //priors
  gamma_mean ~ normal(0,1);
  gamma_sd ~ normal(0,1);
  
  tau_mean ~ normal(0,1);
  tau_sd ~ normal(0,1);
  
  theta_mean ~ normal(0,1);
  theta_sd ~ normal(0,1);
  
  kappa_mean ~ normal(0,1);
  kappa_sd ~ normal(0,1);
  
  alpha_mean ~ normal(0,1);
  alpha_sd ~ normal(0,1);
  
  eps_mean ~ normal(0,1);
  eps_sd ~ normal(0,1);
  
  gamma_raw ~ normal(0,1);
  tau_raw ~ normal(0,1);
  theta_raw ~ normal(0,1);
  kappa_raw ~ normal(0,1);
  alpha_raw ~ normal(0,1);
  eps_raw ~ normal(0,1);

  //likelihood
  y~bernoulli(p_a);
}
