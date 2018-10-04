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
    sm_a = (1/gamma[subj[i]]) * log1p(gamma[subj[i]]*m_a[i]);
    sm_b = (1/gamma[subj[i]]) * log1p(gamma[subj[i]]*m_b[i]);
    Qm = sm_a - sm_b;
    
    sd_a = (1/tau[subj[i]]) * log1p(tau[subj[i]]*d_a[i]);
    sd_b = (1/tau[subj[i]]) * log1p(tau[subj[i]]*d_b[i]);
    Qd = (kappa[subj[i]]/alpha[subj[i]]) * log1p(alpha[subj[i]]* pow((sd_a-sd_b)/theta[subj[i]],theta[subj[i]]));
    
    p_a[i] = pow(Qm,1/eps[subj[i]]) / (pow(Qm,1/eps[subj[i]]) + pow(Qd,1/eps[subj[i]]) + 1e-10);
  }

  
  //priors
  gamma ~ normal(0,5);
  tau ~ normal(0,5);
  theta ~ normal(1,5);
  kappa ~ normal(0,5);
  alpha ~ normal(0,5);
  eps ~ normal(0,5);

  //likelihood
  y~bernoulli(p_a);
}
