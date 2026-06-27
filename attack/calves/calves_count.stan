
data {
  int N;
  int N_fol;
  int N_years;
  int K;
  //int B;
  
  int nac[N];
  matrix[N, K] X_behavc;
  int year_id[N];
  int follow_id[N];
  //matrix[N_years, B] spline_mfit;
  
  // priors (_p for probability, l for lambda)
  real prior_intercepts_sd_l;
  //real prior_lin_sd_l;
  //real prior_spline_sd_l;
  real prior_years_sd_l;
  real prior_fol_sd_l;
}

/*
transformed data {
  matrix[B, N_years] spline_mfit_t = spline_mfit';
}
*/

parameters {
  // spline coefficients
  // matrix[K, B] gamma_coef_raw; 
  // behaviour effects
  vector[K] gamma_raw;
  
  // random effects
  vector[N_fol] e_fol_l_raw;
  vector[N_years] e_year_l_raw;
  
  // sds spline
  //real<lower=0> sd_gamma_raw;
  
  // sds random effects
  real<lower=0> sd_fol_l_raw;
  real<lower=0> sd_year_l_raw;

}

transformed parameters {
  //// sds splines
  //real<lower=0> sd_gamma = sd_gamma_raw * prior_spline_sd_l;

  //// sds random effects
  real<lower=0> sd_fol_l = sd_fol_l_raw * prior_fol_sd_l;
  real<lower=0> sd_year_l = sd_year_l_raw * prior_years_sd_l;

  //// spline coefficients
  //matrix[K, B] gamma;
  vector[K] gamma;
  
  //// spline evaluated at observed years
  //matrix[K, N_years] gamma; 


  //// random effects
  vector[N_fol] e_fol_l = e_fol_l_raw * sd_fol_l; 
  vector[N_years] e_year_l = e_year_l_raw * sd_year_l;
 
  //// Compute spline params
  // spline coefficients for non linear terms
  //gamma_coef[, 2:(B-1)] = gamma_coef_raw[, 2:(B-1)] * sd_gamma;
  // spline intercepts
  //gamma_coef[, 1] = gamma_coef_raw[, 1] * prior_intercepts_sd_l;
  // spline slopes (linear terms)
  //gamma_coef[, B] = gamma_coef_raw[, B] * prior_lin_sd_l;

  //// Evaluate splines in observed years
  gamma = gamma_raw * prior_intercepts_sd_l;

}

model {
  
  //// Priors
  
  // spline coefficients
  //to_vector(gamma_coef_raw) ~ std_normal();
  gamma_raw ~ std_normal(); 
  
  // random effects
  e_fol_l_raw ~ std_normal();
  e_year_l_raw ~ std_normal(); 
  
  // sds splines
  //sd_gamma_raw ~ std_normal();
  
  // sds random effects
  sd_fol_l_raw ~ std_normal();
  sd_year_l_raw ~ std_normal();
  
  //// Likelihood
  
  {
    real lambda;

    // Attack number
    for(i in 1:N) {
      lambda = exp(X_behavc[i, ] * gamma + 
                   e_fol_l[follow_id[i]] + 
                   e_year_l[year_id[i]]);
      nac[i] ~ poisson(lambda) T[1, ]; // includes 1
    }   
    
  }
}

generated quantities {
  vector[K] lambda_not_trunc; // lambda from the untruncated distribution
  vector[K] means;            // mean for the truncated distribution
  real marginal_sd_l = sqrt(sd_fol_l ^ 2 + sd_year_l ^ 2);
  
  lambda_not_trunc = exp(gamma + 0.5 * marginal_sd_l ^ 2);
  for(i in 1:K) 
    means[i] = lambda_not_trunc[i] / (1 - exp(-lambda_not_trunc[i]));
    
}

