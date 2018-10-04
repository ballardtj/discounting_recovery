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
  
   vector<lower=0,upper=1>[Nsubj] omega;
   vector<lower=0>[Nsubj] beta;
   vector<lower=0>[Nsubj] delta;
   vector<lower=0>[Nsubj] sigma;
  
}

model {
  //model
  vector[Ntotal] p_a_logit;
  real u_a;
  real u_b;

  for(i in 1:Ntotal){
    u_a = m_a[i] * (omega[subj[i]] * exp(-beta[subj[i]]*d_a[i]) + (1-omega[subj[i]])*exp(-delta[subj[i]]*d_a[i]));
    u_b = m_b[i] * (omega[subj[i]] * exp(-beta[subj[i]]*d_b[i]) + (1-omega[subj[i]])*exp(-delta[subj[i]]*d_b[i]));
  
    p_a_logit[i] = (u_a - u_b) / sigma[subj[i]];
  }

  //priors
  //omega ~ uniform(0,1);
  beta ~ normal(0,5);
  delta ~ normal(0,5);
  sigma ~ normal(0,5);
  
  //likelihood
  y~bernoulli_logit(p_a_logit);
}
