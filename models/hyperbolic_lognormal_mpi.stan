functions {
  
    vector likelihood( vector phi,               //the sequence of parameters shared across shards,
                       vector theta,             //the sequence of parameters specific to this shard,
                       real[] real_data,         //sequence of real-valued data
                       int[] int_data            //sequence of integer data
                      ){

          //get subject propoerties and choices
          int Nvalid = int_data[1];
          int Nplaces = int_data[2]; //total number of elements in array (including 
         
          //uncenter parameters
          real k = theta[1];
          real sigma = phi[3] + phi[4]*theta[2];
          
          //unpack data
          
          int y[Nvalid];
          vector[Nvalid] p_a_logit;
          real m_a;
          real m_b;
          real d_a;
          real d_b;
          real u_a;
          real u_b;
          real lp;

          for(i in 1:Nvalid){
            m_a = real_data[i];
            m_b = real_data[Nplaces+i];
            d_a = real_data[2*Nplaces+i];
            d_b = real_data[3*Nplaces+i];
            
            u_a = m_a / (1+ k * d_a); //utility of option a
            u_b = m_b / (1+ k * d_b); //utility of option b
            p_a_logit[i] = (u_a-u_b) * sigma; //probability of selecting option a
            
            y[i] = int_data[2+i];
          }
          
          lp = bernoulli_logit_lpmf(y | p_a_logit);

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
  
  real k_shape;
  real<lower=0> k_scale;
  vector<lower=0>[Nsubj] k;

  real<lower=0> sigma_mean;
  real<lower=0> sigma_sd;
  vector<lower=0>[Nsubj] sigma_raw;

}

transformed parameters {
    
    vector[4] phi;
    vector[2] theta[Nsubj];

    //insert hyperpriors into phi vector
    phi[1] = k_shape;
    phi[2] = k_scale;
    phi[3] = sigma_mean;
    phi[4] = sigma_sd;
    
    //insert priors into theta array of vectors
    for(subj in 1:Nsubj){
      theta[subj,1] = k[subj];
      theta[subj,2] = sigma_raw[subj];
    }
}


model {
 
  //priors
  k_shape ~ normal(0,1);
  k_scale ~ normal(0,1);
  
  sigma_mean ~ normal(0,1);
  sigma_sd ~ normal(0,1);

  k ~ lognormal(k_shape,k_scale);
  sigma_raw ~ normal(0,1);
  
  //likelihood
  target += sum(map_rect(likelihood,phi,theta,real_data,int_data));
}
