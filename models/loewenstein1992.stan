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
  
  real<lower=0> alpha_mean;
  real<lower=0> alpha_sd;
  vector<lower=0>[Nsubj] alpha_raw;

  real<lower=0> beta_mean;
  real<lower=0> beta_sd;
  vector<lower=0>[Nsubj] beta_raw;

  real<lower=0> sigma_mean;
  real<lower=0> sigma_sd;
  vector<lower=0>[Nsubj] sigma_raw;
  
}


model {
  //model
   vector[Ntotal] p_a_logit;
   vector[Ntotal] alpha = alpha_raw[subj]*alpha_sd + alpha_mean;
   vector[Ntotal] beta = beta_raw[subj]*beta_sd + beta_mean;
   vector[Ntotal] sigma = sigma_raw[subj]*sigma_sd + sigma_mean;
   real u_a;
   real u_b;

   for(i in 1:Ntotal){
     u_a = m_a[i] * pow(1+alpha[i] * d_a[i], -beta[i] / alpha[i] ) ; //utility of option a
     u_b = m_b[i] * pow(1+alpha[i] * d_b[i], -beta[i] / alpha[i] ) ; //utility of option b
     p_a_logit[i] = (u_a-u_b) * sigma[i]; //probability of selecting option a
   }

  //priors
  alpha_mean ~ normal(0,1);
  alpha_sd ~ normal(0,1);
  
  beta_mean ~ normal(0,1);
  beta_sd ~ normal(0,1);
  
  sigma_mean ~ normal(0,1);
  sigma_sd ~ normal(0,1);

  alpha_raw ~ normal(0,1);
  beta_raw ~ normal(0,1);
  sigma_raw ~ normal(0,1);
  
  //likelihood
  y~bernoulli_logit(p_a_logit);
}
