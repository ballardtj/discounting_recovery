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
  
  real<lower=0> alpha_a;
  real<lower=0> alpha_b;
  vector<lower=0,upper=1>[Nsubj] alpha;
  
  real<lower=0> beta_a;
  real<lower=0> beta_b;
  vector<lower=0,upper=1>[Nsubj] beta;
  
  real<lower=0> lambda_mean;
  real<lower=0> lambda_sd;
  vector<lower=0>[Nsubj] lambda_raw;

  real<lower=0> sigma_mean;
  real<lower=0> sigma_sd;
  vector<lower=0>[Nsubj] sigma_raw;
  
}


model {
  //model
   vector[Ntotal] p_a_logit;
   vector[Ntotal] lambda = lambda_raw[subj]*lambda_sd + lambda_mean;
   vector[Ntotal] sigma = sigma_raw[subj]*sigma_sd + sigma_mean;
   real u_a;
   real u_b;

   for(i in 1:Ntotal){
     u_a = pow(m_a[i],alpha[subj[i]]) - lambda[i]*pow(d_a[i],beta[subj[i]]); //utility of option a
     u_b = pow(m_b[i],alpha[subj[i]]) - lambda[i]*pow(d_b[i],beta[subj[i]]); //utility of option b
     p_a_logit[i] = (u_a-u_b) * sigma[i]; //probability of selecting option a
   }

  //priors
  alpha_a ~ normal(0,1);
  alpha_b ~ normal(0,1);
  
  beta_a ~ normal(0,1);
  beta_b ~ normal(0,1);
  
  lambda_mean ~ normal(0,1);
  lambda_sd ~ normal(0,1);
  
  sigma_mean ~ normal(0,1);
  sigma_sd ~ normal(0,1);

  alpha ~ beta(alpha_a,alpha_b);
  beta ~ beta(beta_a,beta_b);
  lambda_raw ~ normal(0,1);
  sigma_raw ~ normal(0,1);
  
  //likelihood
  y~bernoulli_logit(p_a_logit);
}
