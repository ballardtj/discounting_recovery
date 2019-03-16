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
  
  vector<lower=0,upper=1>[Nsubj] alpha;
  vector<lower=0,upper=1>[Nsubj] beta;
  vector<lower=0>[Nsubj] lambda;
  vector<lower=0>[Nsubj] sigma;
  
}


model {
  //model
   vector[Ntotal] p_a_logit;
   real u_a;
   real u_b;

   for(i in 1:Ntotal){
     u_a = pow(m_a[i],alpha[subj[i]]) - lambda[subj[i]]*pow(d_a[i],beta[subj[i]]); //utility of option a
     u_b = pow(m_b[i],alpha[subj[i]]) - lambda[subj[i]]*pow(d_b[i],beta[subj[i]]); //utility of option b
     p_a_logit[i] = (u_a-u_b) * sigma[subj[i]]; //probability of selecting option a
   }

  //priors
  lambda ~ normal(0,5);
  sigma ~ normal(0,5);
  
  //likelihood
  y~bernoulli_logit(p_a_logit);
}
