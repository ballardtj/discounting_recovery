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
  
  vector<lower=0>[Nsubj] beta;
  vector<lower=0>[Nsubj] lambda;
  vector<lower=0>[Nsubj] sigma;
  
}


model {
  //model
   vector[Ntotal] p_a_logit;
   real u_a;
   real u_b;

   for(i in 1:Ntotal){
     u_a = m_a[i] - lambda[i]*pow(d_a[i],beta[i]); //utility of option a
     u_b = m_b[i] - lambda[i]*pow(d_b[i],beta[i]); //utility of option b
     p_a_logit[i] = (u_a-u_b) * sigma[i]; //probability of selecting option a
   }

  //priors
  beta ~ normal(0,5);
  lambda ~ normal(0,5);
  sigma ~ normal(0,5);
  
  //likelihood
  y~bernoulli_logit(p_a_logit);
}
