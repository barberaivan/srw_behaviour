
data {
  int N;
  int N_fol;
  int N_years;
  int K;
  int B;
  
  int am[N];
  vector[N] am_prev;
  matrix[N, K] X_behavm;
  int year_id[N];
  int follow_id[N];

  matrix[N_years, B] spline_mfit;
  
  // priors (_p for probability, l for lambda)
  real prior_intercepts_sd_p;
  real prior_lin_sd_p;
  real prior_spline_sd_p;
  real prior_years_sd_p;
  real prior_fol_sd_p;

}

transformed data {
  matrix[B, N_years] spline_mfit_t = spline_mfit';
}

parameters {
  // spline coefficients
  matrix[K, B] alpha_coef_raw; // for p | am[t-1] = 0
  matrix[K, B] beta_coef_raw;  // for p | am[t-1] = 1

  // random effects
  vector[N_fol] e_fol_p_raw;
  vector[N_years] e_year_p_raw;
  
  // sds splines
  real<lower=0> sd_alpha_raw;
  real<lower=0> sd_beta_raw;

  // sds random effects
  real<lower=0> sd_fol_p_raw;
  real<lower=0> sd_year_p_raw;

}

transformed parameters {
  //// sds splines
  real<lower=0> sd_alpha = sd_alpha_raw * prior_spline_sd_p;
  real<lower=0> sd_beta = sd_beta_raw * prior_spline_sd_p;

  //// sds random effects
  real<lower=0> sd_fol_p = sd_fol_p_raw * prior_fol_sd_p;
  real<lower=0> sd_year_p = sd_year_p_raw * prior_years_sd_p;

  //// spline coefficients
  matrix[K, B] alpha_coef; // for p | am[t-1] = 0
  matrix[K, B] beta_coef;  // for p | am[t-1] = 1

  //// spline evaluated at observed years
  matrix[K, N_years] alpha; // for p | am[t-1] = 0
  matrix[K, N_years] beta;  // for p | am[t-1] = 1

  //// random effects
  // follow independent for p
  vector[N_fol] e_fol_p = e_fol_p_raw * sd_fol_p; 
  // year
  vector[N_years] e_year_p = e_year_p_raw * sd_year_p;
 
  //// Compute spline params
  // spline coefficients for non linear terms
  alpha_coef[, 2:(B-1)] = alpha_coef_raw[, 2:(B-1)] * sd_alpha;
  beta_coef[, 2:(B-1)] = beta_coef_raw[, 2:(B-1)] * sd_beta;
  // spline intercepts
  alpha_coef[, 1] = alpha_coef_raw[, 1] * prior_intercepts_sd_p;
  beta_coef[, 1] = beta_coef_raw[, 1] * prior_intercepts_sd_p;
  // spline slopes (linear terms)
  alpha_coef[, B] = alpha_coef_raw[, B] * prior_lin_sd_p;
  beta_coef[, B] = beta_coef_raw[, B] * prior_lin_sd_p;

  //// Evaluate splines in observed years
  alpha = alpha_coef * spline_mfit_t;
  beta = beta_coef * spline_mfit_t;

}

model {
  
  //// Priors
  
  // spline coefficients
  to_vector(alpha_coef_raw) ~ std_normal();
  to_vector(beta_coef_raw) ~ std_normal();

  // random effects
  e_fol_p_raw ~ std_normal();
  e_year_p_raw ~ std_normal(); 
  
  // sds splines
  sd_alpha_raw ~ std_normal();
  sd_beta_raw ~ std_normal();
  
  // sds random effects
  sd_fol_p_raw ~ std_normal();
  sd_year_p_raw ~ std_normal();
  
  //// Likelihood
  
  {
    real logit_p;
  
    for(i in 1:N) {
      logit_p = X_behavm[i, ] * alpha[, year_id[i]] * (1 - am_prev[i]) + // there was no attack before
                X_behavm[i, ] * beta[, year_id[i]] * am_prev[i] +        // there was attack before
                e_fol_p[follow_id[i]] + 
                e_year_p[year_id[i]];
      am[i] ~ bernoulli_logit(logit_p);
    }
 
  }
}

generated quantities {
  real marginal_sd_p = sqrt(sd_fol_p ^ 2 + sd_year_p ^ 2);
}

