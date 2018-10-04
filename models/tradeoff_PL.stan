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
  
   vector<lower=0>[Nsubj] gamma;
   vector<lower=0>[Nsubj] tau;
   vector<lower=1>[Nsubj] theta;
   vector<lower=0>[Nsubj] kappa;
   vector<lower=0>[Nsubj] alpha;
   vector<lower=0>[Nsubj] eps;

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
  gamma_raw ~ normal(0,5);
  tau_raw ~ normal(0,5);
  theta_raw ~ normal(1,5);
  kappa_raw ~ normal(0,5);
  alpha_raw ~ normal(0,5);
  eps_raw ~ normal(0,5);

  //likelihood
  y~bernoulli(p_a);
}
