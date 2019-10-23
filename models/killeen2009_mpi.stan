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
          real lambda = theta[3];
          real sigma = theta[4];
          
          
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
  
  real<lower=0,upper=1> alpha_mean;
  real<lower=0> alpha_prec;
  vector<lower=0,upper=1>[Nsubj] alpha;
  
  real<lower=0,upper=1> beta_mean;
  real<lower=0> beta_prec;
  vector<lower=0,upper=1>[Nsubj] beta;
  
  real<lower=0> lambda_shape;
  real<lower=0> lambda_scale;
  vector<lower=0>[Nsubj] lambda;

  real<lower=0> sigma_shape;
  real<lower=0> sigma_scale;
  vector<lower=0>[Nsubj] sigma;
  
}




transformed parameters {
    
    vector[4] phi;
    vector[4] theta[Nsubj];

    //insert hyperpriors into phi vector
    phi[1] = 1;
    
    //insert priors into theta array of vectors
    for(subj in 1:Nsubj){
      theta[subj,1] = alpha[subj];
      theta[subj,2] = beta[subj];
      theta[subj,3] = lambda[subj];
      theta[subj,4] = sigma[subj];
    }
}


model {
 
  //priors
  //alpha_mean ~ uniform(0,1);
  alpha_prec ~ gamma(1,20);
  
  //beta_u ~ uniform(0,1);
  beta_prec ~ gamma(1,20);
  
  lambda_shape ~ normal(0,1);
  lambda_scale ~ normal(0,1);
  
  sigma_shape ~ normal(0,1);
  sigma_scale ~ normal(0,1);

  alpha ~ beta(alpha_mean*alpha_prec , (1-alpha_mean)*alpha_prec);
  beta ~  beta(beta_mean*beta_prec , (1-beta_mean)*beta_prec);
  lambda ~ gamma(lambda_shape,inv(lambda_scale));
  sigma ~ gamma(sigma_shape,inv(sigma_scale));
  
  //likelihood
  target += sum(map_rect(likelihood,phi,theta,real_data,int_data));
}
