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
  
  real<lower=0> k_mean;
  real<lower=0> k_sd;
  vector<lower=0>[Nsubj] k_raw;

  real<lower=0> s_mean;
  real<lower=0> s_sd;
  vector<lower=0>[Nsubj] s_raw;

  real<lower=0> sigma_mean;
  real<lower=0> sigma_sd;
  vector<lower=0>[Nsubj] sigma_raw;
  
 

}

model {
  //model
   vector[Ntotal] p_a_logit;
   vector[Ntotal] k = k_raw[subj]*k_sd + k_mean;
   vector[Ntotal] s = s_raw[subj]*s_sd + s_mean;
   vector[Ntotal] sigma = sigma_raw[subj]*sigma_sd + sigma_mean;
   real u_a;
   real u_b;

   for(i in 1:Ntotal){
     u_a = m_a[i] / pow(1+ k[i] * d_a[i],s[i] ); //utility of option a
     u_b = m_b[i] / pow(1+ k[i] * d_b[i],s[i] ); //utility of option b
     p_a_logit[i] = (u_a-u_b) .* sigma[i]; //probability of selecting option a
   }

  //priors
  k_mean ~ normal(0,1);
  k_sd ~ normal(0,1);
  
  s_mean ~ normal(0,1);
  s_sd ~ normal(0,1);
  
  sigma_mean ~ normal(0,1);
  sigma_sd ~ normal(0,1);

  k_raw ~ normal(0,1);
  s_raw ~ normal(0,1);
  sigma_raw ~ normal(0,1);
  
  //likelihood
  y~bernoulli_logit(p_a_logit);
}
