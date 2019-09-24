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
          real delta = theta[1];
          real sigma = theta[2];
      
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
  
  vector[Nsubj] delta;
  vector<lower=0>[Nsubj] sigma;

}



transformed parameters {
    
    vector[2] theta[Nsubj];
    vector[1] phi;
    phi[1] = 1;

    //insert priors into theta array of vectors
    for(subj in 1:Nsubj){
      theta[subj,1] = delta[subj];
      theta[subj,2] = sigma[subj];
    }
}


model {
 
  //priors
  sigma ~ normal(0,5);
  delta ~ normal(0,5);
  
  
  //likelihood
  target += sum(map_rect(likelihood,phi,theta,real_data,int_data));
}
