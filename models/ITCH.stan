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
  vector[Ntotal] mean_x; 
  vector[Ntotal] mean_t;
  vector[Ntotal] xA; 
  vector[Ntotal] xR;
  vector[Ntotal] tA; 
  vector[Ntotal] tR;
  
  for(i in 1:Ntotal){
    mean_x[i] = fmax( (m_a[i] + m_b[i]) / 2 , 1e-10 );
    mean_t[i] = fmax( (d_a[i] + d_b[i]) / 2 , 1e-10 );
  }
  
  xA = m_a - m_b;
  xR = (m_a - m_b) ./ mean_x;
  tA = d_a - d_b;
  tR = (d_a - d_b) ./ mean_t;
  
}
 
parameters {
  
  real B1_mean;
  real<lower=0> B1_sd;
  vector[Nsubj] B1_raw;

  real<lower=0> BxA_mean;
  real<lower=0> BxA_sd;
  vector<lower=0>[Nsubj] BxA_raw;

  real<lower=0> BxR_mean;
  real<lower=0> BxR_sd;
  vector<lower=0>[Nsubj] BxR_raw;
  
  real<upper=0> BtA_mean;
  real<lower=0> BtA_sd;
  vector<upper=0>[Nsubj] BtA_raw;

  real<upper=0> BtR_mean;
  real<lower=0> BtR_sd;
  vector<upper=0>[Nsubj] BtR_raw;
  
}

model {
  //model
   vector[Ntotal] B1 = B1_raw[subj]*B1_sd + B1_mean;
   vector[Ntotal] BxA = BxA_raw[subj]*BxA_sd + BxA_mean;
   vector[Ntotal] BxR = BxR_raw[subj]*BxR_sd + BxR_mean;
   vector[Ntotal] BtA = BtA_raw[subj]*BtA_sd + BtA_mean;
   vector[Ntotal] BtR = BtR_raw[subj]*BtR_sd + BtR_mean;
   vector[Ntotal] p_a_logit = B1 + BxA .*xA + BxR .*xR + BtA .*tA + BxR .*xR;
   
  //priors
  B1_mean ~ normal(0,1);
  B1_sd ~ normal(0,1);
  
  BxA_mean ~ normal(0,1);
  BxA_sd ~ normal(0,1);
  
  BxR_mean ~ normal(0,1);
  BxR_sd ~ normal(0,1);
  
  BtA_mean ~ normal(0,1);
  BtA_sd ~ normal(0,1);
  
  BtR_mean ~ normal(0,1);
  BtR_sd ~ normal(0,1);

  B1_raw ~ normal(0,1);
  BxA_raw ~ normal(0,1);
  BxR_raw ~ normal(0,1);
  BtA_raw ~ normal(0,1);
  BtR_raw ~ normal(0,1);
  
  //likelihood
  y~bernoulli_logit(p_a_logit);
}
