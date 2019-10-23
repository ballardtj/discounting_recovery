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
          real beta_on_alpha = theta[2];
          real sigma = theta[3];
          
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
            
            u_a = m_a * pow(1+alpha * d_a, -beta_on_alpha) ; //utility of option a
            u_b = m_b * pow(1+alpha * d_b, -beta_on_alpha) ; //utility of option b
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
  
  vector<lower=0>[Nsubj] alpha;
  vector<lower=0>[Nsubj] beta_on_alpha;
  vector<lower=0>[Nsubj] sigma;
  
}




transformed parameters {
    
    vector[3] theta[Nsubj];
    vector[1] phi;
    phi[1] = 1;

    //insert priors into theta array of vectors
    for(subj in 1:Nsubj){
      theta[subj,1] = alpha[subj];
      theta[subj,2] = beta_on_alpha[subj];
      theta[subj,3] = sigma[subj];
    }
}


model {
 
//priors
alpha ~ gamma(1,0.5);
beta_on_alpha ~ gamma(1,0.5);
sigma ~ gamma(1,0.5);
  
  
  
  //likelihood
  target += sum(map_rect(likelihood,phi,theta,real_data,int_data));
}
