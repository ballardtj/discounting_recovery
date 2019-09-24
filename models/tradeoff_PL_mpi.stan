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
          real gamma = theta2[1];
          real tau = theta2[2];
          real theta = theta2[3];
          real kappa = theta2[4];
          real alpha = theta2[5];
          real eps = theta2[6];
          
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
  
   vector<lower=0>[Nsubj] gamma;
   vector<lower=0>[Nsubj] tau;
   vector<lower=1>[Nsubj] theta;
   vector<lower=0>[Nsubj] kappa;
   vector<lower=0>[Nsubj] alpha;
   vector<lower=0>[Nsubj] eps;

}

transformed parameters {
    
    vector[6] theta2[Nsubj];
    vector[1] phi;
    phi[1] = 1;

    //insert priors into theta array of vectors
    for(subj in 1:Nsubj){
      theta2[subj,1] = gamma[subj];
      theta2[subj,2] = tau[subj];
      theta2[subj,3] = theta[subj];
      theta2[subj,4] = kappa[subj];
      theta2[subj,5] = alpha[subj];
      theta2[subj,6] = eps[subj];
    }
}


model {
 
 //priors
  gamma ~ normal(0,5);
  tau ~ normal(0,5);
  theta ~ normal(1,5);
  kappa ~ normal(0,5);
  alpha ~ normal(0,5);
  eps ~ normal(0,5);

  //likelihood
  target += sum(map_rect(likelihood,phi,theta2,real_data,int_data));
}
