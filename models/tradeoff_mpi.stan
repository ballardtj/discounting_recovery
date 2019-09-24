functions {
  
    vector likelihood( vector phi,               //the sequence of parameters shared across shards,
                       vector theta2,             //the sequence of parameters specific to this shard,
                       real[] real_data,         //sequence of real-valued data
                       int[] int_data            //sequence of integer data
                      ){

          //get subject propoerties and choices
          int Nvalid = int_data[1];
          int Nplaces = int_data[2]; //total number of elements in array (including 
         
          //uncenter parameters (only those that are normally distributed)
          real gamma = phi[1] + phi[2]*theta2[1];
          real tau = phi[3] + phi[4]*theta2[2];
          real theta = phi[5] + phi[6]*theta2[3];
          real kappa = phi[7] + phi[8]*theta2[4];
          real alpha = phi[9] + phi[10]*theta2[5];
          real eps = phi[11] + phi[12]*theta2[6];
          
          //unpack data
          
          int y[Nvalid];
          vector[Nvalid] p_a;
          real m_a;
          real m_b;
          real d_a;
          real d_b;
          real sm_a;
          real sm_b;
          real Qm;
          real Qd;
          real sd_a;
          real sd_b;
          real lp;

          for(i in 1:Nvalid){
            m_a = real_data[i];
            m_b = real_data[Nplaces+i];
            d_a = real_data[2*Nplaces+i];
            d_b = real_data[3*Nplaces+i];
            
            sm_a = (1/gamma) * log1p(gamma*m_a);
            sm_b = (1/gamma) * log1p(gamma*m_b);
            Qm = sm_a - sm_b;
    
            sd_a = (1/tau) * log1p(tau*d_a);
            sd_b = (1/tau) * log1p(tau*d_b);
            Qd = (kappa/alpha) * log1p(alpha* pow((sd_a-sd_b)/theta,theta));
    
            p_a[i] = pow(Qm,1/eps) / (pow(Qm,1/eps) + pow(Qd,1/eps) + 1e-10);
            
            y[i] = int_data[2+i];
          }
          
          lp = bernoulli_lpmf(y | p_a);

	        return [lp]';
  }
}

data {
 int Nsubj;
 int Max_obs_per_subj;
 real real_data[Nsubj,Max_obs_per_subj*4];
 int int_data[Nsubj,Max_obs_per_subj+2]; 
}  

parameters {
  
  real<lower=0> gamma_mean;
  real<lower=0> gamma_sd;
   vector<lower=0>[Nsubj] gamma_raw;

  real<lower=0> tau_mean;
  real<lower=0> tau_sd;
   vector<lower=0>[Nsubj] tau_raw;
  
  real<lower=1> theta_mean;
  real<lower=0> theta_sd;
   vector<lower=0>[Nsubj] theta_raw;
  
  real<lower=0> kappa_mean;
  real<lower=0> kappa_sd;
   vector<lower=0>[Nsubj] kappa_raw;
  
  real<lower=0> alpha_mean;
  real<lower=0> alpha_sd;
   vector<lower=0>[Nsubj] alpha_raw;
  
  real<lower=0> eps_mean;
  real<lower=0> eps_sd;
  vector<lower=0>[Nsubj] eps_raw;

}




transformed parameters {
    
    vector[12] phi;
    vector[6] theta2[Nsubj];

    //insert hyperpriors into phi vector
    phi[1] = gamma_mean;
    phi[2] = gamma_sd;
    phi[3] = tau_mean;
    phi[4] = tau_sd;
    phi[5] = theta_mean;
    phi[6] = theta_sd;
    phi[7] = kappa_mean;
    phi[8] = kappa_sd;
    phi[9] = alpha_mean;
    phi[10] = alpha_sd;
    phi[11] = eps_mean;
    phi[12] = eps_sd;
  
    //insert priors into theta array of vectors
    for(subj in 1:Nsubj){
      theta2[subj,1] = gamma_raw[subj];
      theta2[subj,2] = tau_raw[subj];
      theta2[subj,3] = theta_raw[subj];
      theta2[subj,4] = kappa_raw[subj];
      theta2[subj,5] = alpha_raw[subj];
      theta2[subj,6] = eps_raw[subj];
    }
}


model {
 
  //priors
  gamma_mean ~ normal(0,1);
  gamma_sd ~ normal(0,1);
  
  tau_mean ~ normal(0,1);
  tau_sd ~ normal(0,1);
  
  theta_mean ~ normal(1,1);
  theta_sd ~ normal(0,1);
  
  kappa_mean ~ normal(0,1);
  kappa_sd ~ normal(0,1);
  
  alpha_mean ~ normal(0,1);
  alpha_sd ~ normal(0,1);
  
  eps_mean ~ normal(0,1);
  eps_sd ~ normal(0,1);
  
  gamma_raw ~ normal(0,1);
  tau_raw ~ normal(0,1);
  theta_raw ~ normal(0,1);
  kappa_raw ~ normal(0,1);
  alpha_raw ~ normal(0,1);
  eps_raw ~ normal(0,1);

  
  
  //likelihood
  target += sum(map_rect(likelihood,phi,theta2,real_data,int_data));
}
