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
          real s = theta[2];
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
            
            u_a = m_a / (1+ k * pow(d_a,s ) ); //utility of option a
            u_b = m_b / (1+ k * pow(d_b,s ) ); //utility of option b
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
  
  vector<lower=0>[Nsubj] k;
  vector<lower=0>[Nsubj] s;
  vector<lower=0>[Nsubj] sigma;
  
}

transformed parameters {
    
    vector[3] theta[Nsubj];
    vector[1] phi;
    phi[1] = 1;

    //insert priors into theta array of vectors
    for(subj in 1:Nsubj){
      theta[subj,1] = k[subj];
      theta[subj,2] = s[subj];
      theta[subj,3] = sigma[subj];
    }
}


model {
 
   //priors
  k ~ gamma(1,0.5);
  s ~ normal(1,0.5);
  sigma ~ normal(1,0.5);
  
  //likelihood
  target += sum(map_rect(likelihood,phi,theta,real_data,int_data));
}
