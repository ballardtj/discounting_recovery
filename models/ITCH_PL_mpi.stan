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
          real B1 = theta[1];
          real BxA = theta[2];
          real BxR = theta[3];
          real BtA = theta[4];
          real BtR = theta[5];
          
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
            
            p_a_logit[i] = B1 + BxA*xA + BxR*xR + BtA*tA + BxR*xR;
            
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
  
  vector[Nsubj] B1;
  vector<lower=0>[Nsubj] BxA;
  vector<lower=0>[Nsubj] BxR;
  vector<upper=0>[Nsubj] BtA;
  vector<upper=0>[Nsubj] BtR;
  
}



transformed parameters {
    
    vector[5] theta[Nsubj];
    vector[1] phi;
    phi[1] = 1;

    //insert priors into theta array of vectors
    for(subj in 1:Nsubj){
      theta[subj,1] = B1[subj];
      theta[subj,2] = BxA[subj];
      theta[subj,3] = BxR[subj];
      theta[subj,4] = BtA[subj];
      theta[subj,5] = BtR[subj];
    }
}


model {
 
  //priors
  B1 ~ normal(0,5);
  BxA ~ normal(0,5);
  BxR ~ normal(0,5);
  BtA ~ normal(0,5);
  BtR ~ normal(0,5);
  
  //likelihood
  target += sum(map_rect(likelihood,phi,theta,real_data,int_data));
}
