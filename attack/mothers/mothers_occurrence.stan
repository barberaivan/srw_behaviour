
data {
  int N;
  int N_fol;
  int N_years;
  int K;
  
  int am_marg[N];
  vector[N] am_marg_prev;
  matrix[N, K] X_behavm;
  int year_id[N];
  int follow_id[N];

  // priors (_p for probability, l for lambda)
  real prior_intercepts_sd_p;
  real prior_years_sd_p;
  real prior_fol_sd_p;

}

parameters {
  // coefficients
  vector[K] alpha_raw; // for p | am_marg[t-1] = 0
  vector[K] beta_raw;  // for p | am_marg[t-1] = 1

  // random effects
  vector[N_fol] e_fol_p_raw;
  vector[N_years] e_year_p_raw;

  // sds random effects
  real<lower=0> sd_fol_p_raw;
  real<lower=0> sd_year_p_raw;

}

transformed parameters {
  //// sds random effects
  real<lower=0> sd_fol_p = sd_fol_p_raw * prior_fol_sd_p;
  real<lower=0> sd_year_p = sd_year_p_raw * prior_years_sd_p;

  //// coefficients
  vector[K] alpha = alpha_raw * prior_intercepts_sd_p; // for p | am_marg[t-1] = 0
  vector[K] beta = beta_raw * prior_intercepts_sd_p;  // for p | am_marg[t-1] = 1

  //// random effects
  // follow 
  vector[N_fol] e_fol_p = e_fol_p_raw * sd_fol_p; 
  // year
  vector[N_years] e_year_p = e_year_p_raw * sd_year_p;


}

model {
  
  //// Priors
  
  // coefficients
  alpha_raw ~ std_normal();
  beta_raw ~ std_normal();

  // random effects
  e_fol_p_raw ~ std_normal();
  e_year_p_raw ~ std_normal(); 
  
  // sds random effects
  sd_fol_p_raw ~ std_normal();
  sd_year_p_raw ~ std_normal();
  
  //// Likelihood
  
  {
    real logit_p;
  
    for(i in 1:N) {
      logit_p = X_behavm[i, ] * alpha * (1 - am_marg_prev[i]) + // there was no attack before
                X_behavm[i, ] * beta * am_marg_prev[i] +        // there was attack before
                e_fol_p[follow_id[i]] + 
                e_year_p[year_id[i]];
      am_marg[i] ~ bernoulli_logit(logit_p);
    }
 
  }
}

generated quantities {
  real marginal_sd_p = sqrt(sd_fol_p ^ 2 + sd_year_p ^ 2);
}

