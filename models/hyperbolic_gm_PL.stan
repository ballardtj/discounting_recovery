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
  
  vector<lower=0>[Nsubj] k;
  vector<lower=0>[Nsubj] s;
  vector<lower=0>[Nsubj] sigma;
  
}


model {
  //model
   vector[Ntotal] p_a_logit;
   real u_a;
   real u_b;

   for(i in 1:Ntotal){
     u_a = m_a[i] / pow(1+ k[subj[i]] * d_a[i],s[subj[i]] ); //utility of option a
     u_b = m_b[i] / pow(1+ k[subj[i]] * d_b[i],s[subj[i]] ); //utility of option b
     p_a_logit[i] = (u_a-u_b) * sigma[subj[i]]; //probability of selecting option a
   }

  //priors
  k ~ gamma(0,0.5);
  s ~ gamma(0,0.5);
  sigma ~ gamma(0,0.5);
  
  //likelihood
  y~bernoulli_logit(p_a_logit);
}
