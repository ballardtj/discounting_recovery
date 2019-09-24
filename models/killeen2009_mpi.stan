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
          real alpha = theta[1];
          real beta = theta[2];
          real lambda = phi[1] + phi[2]*theta[3];
          real sigma = phi[3] + phi[4]*theta[4];
          
          
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
            
            u_a = pow(m_a,alpha) - lambda*pow(d_a,beta); //utility of option a
            u_b = pow(m_b,alpha) - lambda*pow(d_b,beta); //utility of option b
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
  
  real<lower=0> alpha_a;
  real<lower=0> alpha_b;
  vector<lower=0,upper=1>[Nsubj] alpha;
  
  real<lower=0> beta_a;
  real<lower=0> beta_b;
  vector<lower=0,upper=1>[Nsubj] beta;
  
  real<lower=0> lambda_mean;
  real<lower=0> lambda_sd;
  vector<lower=0>[Nsubj] lambda_raw;

  real<lower=0> sigma_mean;
  real<lower=0> sigma_sd;
  vector<lower=0>[Nsubj] sigma_raw;
  
}




transformed parameters {
    
    vector[4] phi;
    vector[4] theta[Nsubj];

    //insert hyperpriors into phi vector
    phi[1] = lambda_mean;
    phi[2] = lambda_sd;
    phi[3] = sigma_mean;
    phi[4] = sigma_sd;
  
    //insert priors into theta array of vectors
    for(subj in 1:Nsubj){
      theta[subj,1] = alpha[subj];
      theta[subj,2] = beta[subj];
      theta[subj,3] = lambda_raw[subj];
      theta[subj,4] = sigma_raw[subj];
    }
}


model {
 
  //priors
  alpha_a ~ normal(0,1);
  alpha_b ~ normal(0,1);
  
  beta_a ~ normal(0,1);
  beta_b ~ normal(0,1);
  
  lambda_mean ~ normal(0,1);
  lambda_sd ~ normal(0,1);
  
  sigma_mean ~ normal(0,1);
  sigma_sd ~ normal(0,1);

  alpha ~ beta(alpha_a,alpha_b);
  beta ~ beta(beta_a,beta_b);
  lambda_raw ~ normal(0,1);
  sigma_raw ~ normal(0,1);
  
  //likelihood
  target += sum(map_rect(likelihood,phi,theta,real_data,int_data));
}
