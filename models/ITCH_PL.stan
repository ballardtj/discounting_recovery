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
  
  vector[Nsubj] B1;
  vector<lower=0>[Nsubj] BxA;
  vector<lower=0>[Nsubj] BxR;
  vector<upper=0>[Nsubj] BtA;
  vector<upper=0>[Nsubj] BtR;
  
}

model {
  //model
   vector[Ntotal] p_a_logit = B1 + BxA[subj].*xA + BxR[subj].*xR + BtA[subj].*tA + BxR[subj].*xR;
   
  //priors
  B1 ~ normal(0,5);
  BxA ~ normal(0,5);
  BxR ~ normal(0,5);
  BtA ~ normal(0,5);
  BtR ~ normal(0,5);
  
  //likelihood
  y~bernoulli_logit(p_a_logit);
}
