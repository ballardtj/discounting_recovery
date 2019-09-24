functions {
  
    vector likelihood( vector phi,               //the sequence of parameters shared across shards,
                       vector theta,             //the sequence of parameters specific to this shard,
                       real[] real_data,         //sequence of real-valued data
                       int[] int_data            //sequence of integer data
                      ){

          //get subject propoerties and choices
          int Nvalid = int_data[1];
          int Nplaces = int_data[2]; //total number of elements in array (including 
         
          //uncenter parameters (only those that are normally distributed)
          real delta = phi[1] + phi[2]*theta[1];
          real sigma = phi[3] + phi[4]*theta[2];
      
          //unpack data
          int y[Nvalid];
          vector[Nvalid] p_a;
          real m_a;
          real m_b;
          real d_a;
          real d_b;
          real alpha_m;
          real alpha_d;
          real pi_m;
          real pi_d;
          real d;
          real lp;

          for(i in 1:Nvalid){
            m_a = real_data[i];
            m_b = real_data[Nplaces+i];
            d_a = real_data[2*Nplaces+i];
            d_b = real_data[3*Nplaces+i];
            
            alpha_m = fmax( fmax( fabs(m_a),fabs(m_b)), 1e-10);
            alpha_d = fmax( fmax( fabs(d_a),fabs(d_b)), 1e-10);
            
            pi_m = (fmax( fabs(m_a),fabs(m_b)) - fmin( fabs(m_a),fabs(m_b)) ) / alpha_m;
            pi_d = (fmax( fabs(d_a),fabs(d_b)) - fmin( fabs(d_a),fabs(d_b)) ) / alpha_d;
    
            d = pi_m - pi_d;
            
            
            
            p_a[i] = Phi_approx((d-delta) * sigma);
            
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
  
  real delta_mean;
  real<lower=0> delta_sd;
  vector[Nsubj] delta_raw;

  real<lower=0> sigma_mean;
  real<lower=0> sigma_sd;
  vector<lower=0>[Nsubj] sigma_raw;

}




transformed parameters {
    
    vector[4] phi;
    vector[2] theta[Nsubj];

    //insert hyperpriors into phi vector
    phi[1] = delta_mean;
    phi[2] = delta_sd;
    phi[3] = sigma_mean;
    phi[4] = sigma_sd;
  
    //insert priors into theta array of vectors
    for(subj in 1:Nsubj){
      theta[subj,1] = delta_raw[subj];
      theta[subj,2] = sigma_raw[subj];
    }
}


model {
 
  //priors
  delta_mean ~ normal(0,1);
  delta_sd ~ normal(0,1);
  
  sigma_mean ~ normal(0,1);
  sigma_sd ~ normal(0,1);
  
  sigma_raw ~ normal(0,1);
  delta_raw ~ normal(0,1);
  
  
  //likelihood
  target += sum(map_rect(likelihood,phi,theta,real_data,int_data));
}
