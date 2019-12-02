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
          real B1 = phi[1] + phi[2]*theta[1];
          real BxA = phi[3] + phi[4]*theta[2];
          real BxR = phi[5] + phi[6]*theta[3];
          real BtA = phi[7] + phi[8]*theta[4];
          real BtR = phi[9] + phi[10]*theta[5];
          
          //unpack data
          
          int y[Nvalid];
          vector[Nvalid] p_a_logit;
          real m_a;
          real m_b;
          real d_a;
          real d_b;
          real mean_x;
          real mean_t;
          real xA;
          real xR;
          real tA;
          real tR;
          real lp;

          for(i in 1:Nvalid){
            m_a = real_data[i];
            m_b = real_data[Nplaces+i];
            d_a = real_data[2*Nplaces+i];
            d_b = real_data[3*Nplaces+i];
            
            mean_x = fmax( (m_a + m_b) / 2 , 1e-10 );
            mean_t = fmax( (d_a + d_b) / 2 , 1e-10 );
            
            xA = m_a - m_b;
            xR = (m_a - m_b) / mean_x;
            tA = d_a - d_b;
            tR = (d_a - d_b) / mean_t;
            
            p_a_logit[i] = B1 + BxA*xA + BxR*xR + BtA*tA + BtR*tR;
            
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



transformed parameters {
    
    vector[10] phi;
    vector[5] theta[Nsubj];

    //insert hyperpriors into phi vector
    phi[1] = B1_mean;
    phi[2] = B1_sd;
    phi[3] = BxA_mean;
    phi[4] = BxA_sd;
    phi[5] = BxR_mean;
    phi[6] = BxR_sd;
    phi[7] = BtA_mean;
    phi[8] = BtA_sd;
    phi[9] = BtR_mean;
    phi[10] = BtR_sd;
    
    //insert priors into theta array of vectors
    for(subj in 1:Nsubj){
      theta[subj,1] = B1_raw[subj];
      theta[subj,2] = BxA_raw[subj];
      theta[subj,3] = BxR_raw[subj];
      theta[subj,4] = BtA_raw[subj];
      theta[subj,5] = BtR_raw[subj];
    }
}


model {
 
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
  target += sum(map_rect(likelihood,phi,theta,real_data,int_data));
}
