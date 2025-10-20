// This is a Stan implementation of the first-difference model
// that does not share any information across strata - intercepts and annual differences are independent among strata

// Consider moving annual index calculations outside of Stan to 
// facilitate the ragged array issues



data {
  int<lower=1> nsites;
  int<lower=1> nstrata;
  int<lower=1> ncounts;
  int<lower=1> nyears;
  int<lower=1> fixed_year; //middle year of the time-series scaled to ~(n_years/2)
  int<lower=1> neffort_preds; //middle year of the time-series scaled to ~(n_years/2)

  array[ncounts] int<lower=0> count;              // count observations
  array[ncounts] int<lower=1> strat;               // strata indicators
  array[ncounts] int<lower=1> year; // year index
  array[ncounts] int<lower=1> site; // site index
  vector[neffort_preds] effort_preds;
  
 array[nstrata] int<lower=1> nsites_strata; // number of sites in each stratum
 int<lower=0> maxnsites_strata; //largest value of nsites_strata

  array[nstrata,maxnsites_strata] int<lower=0> ste_mat; //matrix identifying which sites are in each stratum
  // above is actually a ragged array, but filled with 0 values so that it works
  // but throws an error if an incorrect strata-site combination is called
 
 
   // extra data to support the first-difference time-series implementation, which is centered at the mid-year of the available time
  // data to center the abundance estimate
  int nIy1; //indexing vector dimension - number of years before fixed_year
  int nIy2; //indexing vector dimension - number of years after fixed_year
  array[nIy1] int Iy1;//indexing vector - indices of years before fixed_year
  array[nIy2] int Iy2;//indexing vector - indices of years after fixed_year
  // a vector of zeros to fill fixed beta values for fixed_year
  vector[nstrata] zero_betas;


  // CBC effort values - party_hours scaled to the mean across all surveys (party_hours/mean(party_hours))
  vector[ncounts] hours;


  // circle inclusion scaling value
  array[nstrata] real nonzeroweight;

}


transformed data{
       int<lower=1> n_years_m1 = nyears-1; 
       //BETA_raw(differences) are 1-fewer than n_years
}

parameters {
  vector[nstrata] strata_raw;

  vector[nsites] ste_raw;   // 
  real<lower=0> sdnoise;    // sd of over-dispersion
  real<lower=0> sdste;    // sd of site effects
  vector<lower=0>[nstrata] sdbeta;   // sd of annual changes within strata
  real<lower=0> sdstrata;    // sd of intercepts

  vector[nstrata] b_raw; //stratum effort slopes
  vector[nstrata] p_raw; //stratum effort exponents

  matrix[nstrata,n_years_m1] beta_raw;         // strata level parameters

}

transformed parameters { 
  vector[ncounts] E;           // log_scale additive likelihood
  matrix[nstrata,nyears] beta;         // strata-level mean differences (0-centered deviation from continental mean BETA)
  matrix[nstrata,nyears] yeareffect;  // matrix of estimated annual values of trajectory
  real<lower=0> phi; //transformed sdnoise if use_pois == 0 (and therefore Negative Binomial)

  phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations


 beta[,fixed_year] = zero_betas; //fixed at zero
  yeareffect[,fixed_year] = zero_betas; //fixed at zero

for(s in 1:nstrata){
// first half of time-series - runs backwards from fixed_year
  for(t in Iy1){
    beta[s,t] = (sdbeta[s] * beta_raw[s,t]);// + BETA[t];
    yeareffect[s,t] = yeareffect[s,t+1] - beta[s,t];
  }
// second half of time-series - runs forwards from fixed_year
   for(t in Iy2){
    beta[s,t] = (sdbeta[s] * beta_raw[s,t-1]);// + BETA[t-1];//t-1 indicators to match dimensionality
    yeareffect[s,t] = yeareffect[s,t-1] + beta[s,t];
  }
}

// intercepts and slopes

  


  for(i in 1:ncounts){
    real strata = (sdstrata*strata_raw[strat[i]]);
    real ste = sdste*ste_raw[site[i]]; // site intercepts
    real b = b_raw[strat[i]];  // effort slope for stratum i
    real p = p_raw[strat[i]];  //effort coefficient for stratum i
    real effort_effect = (b*((hours[i]^p)-1))/p; //effort correction for count i
    
    E[i] =  strata + yeareffect[strat[i],year[i]] + effort_effect + ste;
  }
  
  }
  
  
  
model {
  sdnoise ~ student_t(3,0,1); //prior on scale of extra Poisson log-normal variance
  sdste ~ student_t(3,0,1); //prior on sd of site effects
  sdbeta ~ student_t(3,0,0.2); // prior on sd of differences among strata
  sdstrata ~ student_t(3,0,1); //prior on sd of intercept variation

  
  b_raw ~ std_normal();//effort slopes by stratum
 
  p_raw ~ std_normal();//effort exponents by stratum


  ste_raw ~ std_normal();//site effects
  sum(ste_raw) ~ normal(0,0.001*nsites);
 

//spatial iCAR intercepts and differences by strata
for(s in 1:(nstrata)){
    beta_raw[s,] ~ std_normal();
}

   strata_raw ~ std_normal();

  count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation


}

 generated quantities {

   array[nstrata,nyears] real<lower=0> n; //full annual indices
   array[nstrata,neffort_preds] real effort_vis;
  // vector[ncounts] log_lik; // alternative value to track the observervation level log-likelihood
  // potentially useful for estimating loo-diagnostics, such as looic
  
  // for(i in 1:ncounts){
  // log_lik[i] = neg_binomial_2_log(count[i] | E[i],phi);
  // }
  
for(s in 1:nstrata){
    real b = b_raw[s];  // effort slope for stratum i
    real p = p_raw[s];  //effort coefficient for stratum i
 
  for(i in 1:neffort_preds){

    effort_vis[s,i] = (b*((effort_preds[i]^p)-1))/p;
  
      }
  }
for(y in 1:nyears){

      for(s in 1:nstrata){

  array[nsites_strata[s]] real n_t;
  real strata = (sdstrata*strata_raw[s]);

  

        for(t in 1:nsites_strata[s]){

  real ste = sdste*ste_raw[ste_mat[s,t]]; // site intercepts


      n_t[t] = exp(strata + yeareffect[s,y] + ste);
        }
        n[s,y] = mean(n_t) * nonzeroweight[s]; //mean of exponentiated predictions across sites in a stratum
        //using the mean of hte exponentiated values, instead of including the log-normal
        // retransformation factor (0.5*sdste^2), because this retransformation makes 2 questionable assumptions:
          // 1 - assumes that sites are exchangeable among strata - e.g., that sdste is equal among all strata
          // 2 - assumes that the distribution of site-effects is normal



    }
  }



 }

