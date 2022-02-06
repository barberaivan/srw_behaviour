/* 
  In this model there is no year : behav_t-1 interaction, so I add a year
  effect to the complete column in the transition matrix. 
*/

data{
  // Ns and indexes
  int<lower = 1> N;         // full number of obervation intervals
  int<lower = 1> N_cons;    // number of observations with observed current and previous behavior (consecutive)
  int<lower = 2> K;         // number of behaviors
  int<lower = 2> b;         // Basis dimension for splines, WITHOUT INTERCEPT (lowercase means that)
  
  int<lower = 1> N_fol;                // number of follows
  int<lower = 1> int_start[N_fol];     // follows begginings
  int<lower = 1> int_end[N_fol];       // follows ends
  
  int<lower = 1> N_years;              // number of years
  
  // index matrix to parameterize transition matrix with zeroes in the diagonal
  int ind_matrix[K-1, K];    
 
  // Years-related data
  matrix[N_years, b] spline_mfit;      // spline design matrix without intercept
  int<lower = 1> year_start[N_years];  
  int<lower = 1> year_end[N_years];
  // year_start and end identify indexes where each year begin or ends 
  // in the subsetted design matrix, which contains only consecutive observed
  // behaviors
  int<lower = 1> year_length[N_years];
  // number of intervals with consecutive observed behaviors by year
  int<lower = 1> period_id[N]; // period identifier to compute z (1995 or the remaining years)

  // Response variable 
  int<lower = 1, upper = K+1> y_vec[N];   // behavior
  matrix[N, K] y_mat;                     // behavior in matrix form (binary elements)
  int cons_id[N_cons];       // index identifying observations with observed behavior at current and previous interval
  int cons_id_prev[N_cons];  // cons_id - 1 (previous behavior index)  
  
  // attack variables                                  
  vector[N] a;        // attack presence to mother or calf
  vector[N] am_only;  // attack presence only to mother
  vector[N] am_marg;  // attack presence to mother independent of what happens to calf
  vector[N] ac_only;  // attack presence only to mother
  vector[N] ac_marg;  // attack presence to mother independent of what happens to calf
  vector[N] nam;      // number of attacks to mother
  vector[N] nac;      // number of attacks to calf
  

  //// Priors
  
  // Transition matrix 
  real<lower = 0> prior_alpha_intercept_sd;
  real<lower = 0> prior_beta_intercept_sd;
  
  // Standard deviation of transition matrix parameters (spline)
  real<lower = 0> prior_sd_spline_tmat_sd;
  // Standard deviation for linear spline term
  real<lower = 0> prior_sd_spline_lin_tmat_sd;
  
  // Standard deviation for year random effect
  real<lower = 0> prior_sd_years_sd;

  // z parameters prior sd
  real<lower = 0> prior_zpar_sd;
}

transformed data{
  int y_vec_cons[N_cons] = y_vec[cons_id];
  matrix[N_cons, K] y_mat_cons = y_mat[cons_id_prev, ];
  // Behavior at t-1 for consecutive observed behaviors (it's cons_id - 1)
}

parameters{
  //// Transition matrices (off-diagonal elements; diagonal is fixed at 0)
  
  // intercepts
  matrix[K-1, K] alpha_intercept_raw;
  matrix[K-1, K] beta_intercept_raw;
  
  // column-wise random effects for alpha
  matrix[N_years, K] alpha_raneff_raw;
  // alpha raneff (year) sd 
  real<lower = 0> sd_years_raw;
  
  // spline coefficients (without intercept, b = 4)
  matrix[b, K] alpha_spline_coef_raw; 
  matrix[b, K] beta_spline_coef_raw;
  // Standard deviation for tmat spline 
  real<lower = 0> sd_spline_tmat_raw;
  //  matrix[N_years, b] spline_mfit;      // spline design matrix  
  
  //// Latent variable initial value
  vector<lower = 0, upper = 1>[N_fol] z_0;
  
  //// z dynamics with number of attacks
  vector[2] decay_raw;
  
  vector[2] kappa_raw;                    // intercept when both are attacked
  vector<lower = 0>[2] delta_m_raw;       // difference from kappa when only mother is attacked
  real<lower = 0> delta_c_raw;            // difference from kappa when only calf is attacked
  //vector<lower = 0>[2] lambda_m;        // slope for number of attacks to mother
  real<lower = 0> lambda_c_raw;           // slope for number of attacks to calf  
}

transformed parameters{
  
  // Variables definition --------------------------------------------------
  
  //// Transition matrices (off-diagonal elements; diagonal is fixed at 0)
  
  // N_years transition matrices without reference row
  matrix[K-1, K] alpha_sub [N_years];
  matrix[K-1, K] beta_sub [N_years];
  // N_years complete transition matrices
  matrix[K, K] alpha [N_years];
  matrix[K, K] beta [N_years];
  
  // intercepts
  matrix[K-1, K] alpha_intercept = alpha_intercept_raw * prior_alpha_intercept_sd;
  matrix[K-1, K] beta_intercept = beta_intercept_raw * prior_beta_intercept_sd;
  
  // alpha raneff (year) sd 
  real<lower = 0> sd_years = sd_years_raw * prior_sd_years_sd;
  // column-wise random effects for alpha
  matrix[N_years, K] alpha_raneff = alpha_raneff_raw * sd_years;
  
  // Standard deviation for tmat spline 
  real<lower = 0> sd_spline_tmat = sd_spline_tmat_raw * prior_sd_spline_tmat_sd;
  // spline coefficients (without intercept, b = 4)
  matrix[b, K] alpha_spline_coef; 
  matrix[b, K] beta_spline_coef;
  // splines evaluated at the target years
  matrix[N_years, K] alpha_splines;// = spline_mfit * alpha_spline_coef; 
  matrix[N_years, K] beta_splines;// = spline_mfit * beta_spline_coef; 

  //// z dynamics parameters 
  vector<lower = 0, upper = 1>[N] z;
  vector[2] decay_logit = decay_raw * prior_zpar_sd;
  vector[2] decay = inv_logit(decay_raw * prior_zpar_sd);
  
  vector[2] kappa = kappa_raw * prior_zpar_sd;
  vector[2] delta_m = delta_m_raw * prior_zpar_sd;
  real delta_c = delta_c_raw * prior_zpar_sd;
  real lambda_c = lambda_c_raw * prior_zpar_sd;

  vector[2] increase_b = inv_logit(kappa);
  vector[2] increase_m = inv_logit(kappa - delta_m);
  vector[2] increase_c = inv_logit(kappa - delta_c);
  vector[2] increase_b3 = inv_logit(kappa + 2 * lambda_c); // 3 attacks to calf
  vector[2] increase_c3 = inv_logit(kappa - delta_c + lambda_c * 2);  
  
  //// Design matrix elements
  matrix[N_cons, K] y_mat_cons_z; // half desing matrix
  matrix[N_cons, K*2] X; // full design matrix: cbind(y_mat_cons, y_mat_cons_z)
  
  
  // Transition probabilities parameters -------------------------------------
  
  // Compute actual tmat spline coefficients
  // non linear terms (random effects)
  alpha_spline_coef[1:(b-1), ] = alpha_spline_coef_raw[1:(b-1), ] * sd_spline_tmat; 
  beta_spline_coef[1:(b-1), ] = beta_spline_coef_raw[1:(b-1), ] * sd_spline_tmat; 
  // linear terms (fixed effects)
  alpha_spline_coef[b, ] = alpha_spline_coef_raw[b, ] * prior_sd_spline_lin_tmat_sd; 
  beta_spline_coef[b, ] = beta_spline_coef_raw[b, ] * prior_sd_spline_lin_tmat_sd; 
  
  // Compute splines (N_years x K matrices)
  alpha_splines = spline_mfit * alpha_spline_coef; 
  beta_splines = spline_mfit * beta_spline_coef;

  // Compute alpha_sub and beta_sub
  for(y in 1:N_years) {
    for(k in 1:K) {
      alpha_sub[y, , k] = alpha_intercept[, k] + 
                          alpha_raneff[y, k] +        //matrix[N_years, K]
                          alpha_splines[y, k];        //matrix[N_years, K]
      
      beta_sub[y, , k] = beta_intercept[, k] + 
                         beta_splines[y, k];          //matrix[N_years, K]
    }
  }
  
  
  // Complete matrix of actual parameters by year (zero-diagonal matrices)
  for(y in 1:N_years) {
    for(k in 1:K) {
      // Fill off-diagonal elements
      alpha[y, ind_matrix[, k], k] = alpha_sub[y, , k]; 
      beta[y, ind_matrix[, k], k] = beta_sub[y, , k];       
    
      // Fill the diagonal with zeroes
      alpha[y, k, k] = 0;
      beta[y, k, k] = 0;
    }
  }
  
  
  //// z dynamics ------------------------------------------------------------
  
  /*
     Initialize z in every follow
     the estimated z_0 is the value before the follow, so attack at t = 1
     helps estimating the z at t = 1, the first z value used to compute the
     the likelihood of an observed behavior.
  
     z[t] = f(z[t-1], attacks[t]), so
     transition probability from time [t-1] to [t] = f(z[t-1]).
     This way, 
     transition probability from time [t-1] to [t] = f(attack[t-1]).
     This is because the behaviour is recorded in the beginning of each interval,
     while attacks are recorded during the whole 5 min interval.
     It's like saying that z[t] is the disturbance state at the end of the [t]
     interval.
     
     The design matrix for behaviour[t] is made with 
     behaviour[t-1] and z[t-1].
  */
  
  
  for(f in 1:N_fol) {
    z[int_start[f]] = z_0[f] 
      // increase
      + a[int_start[f]] * (1 - z_0[f]) * inv_logit(
        kappa[period_id[int_start[f]]]                                // intercept 1 attack both
        - delta_m[period_id[int_start[f]]] * am_only[int_start[f]]    // intercept 1 attack mother
        - delta_c * ac_only[int_start[f]]                             // intercept 1 attack calf
        + lambda_c * (nac[int_start[f]] - 1) * ac_marg[int_start[f]]  // effect of number of attacks to calf 
      ) 
      // decay
      - (1 - a[int_start[f]]) * decay[period_id[int_start[f]]] * z_0[f];
  }         
  /* 
    a = attack presence to mother or calf
    am_only = attack presence only to mother
    am_marg = attack presence to mother independent of what happens to calf
    ac_only = attack presence only to calf
    ac_marg = attack presence to calf independent of what happens to mother
    nac = number of attacks to calf
    
    To compute effects of number of attacks, nam is substracted 1, 
    so the intercepts means logit(increase) when there is 1 attack, not 0, 
    because if nac is 0 and am_only = 0, the whole increase term shuts down.
  */  
    
  // The remaining z's
  for(f in 1:N_fol) {
    for(t in (int_start[f] + 1):int_end[f]) {
      z[t] = z[t-1] 
        // increase
        + a[t] * (1 - z[t-1]) * inv_logit(
          kappa[period_id[t]]                       // intercept 1 attack both
          - delta_m[period_id[t]] * am_only[t]      // intercept 1 attack mother
          - delta_c * ac_only[t]                    // intercept 1 attack calf
          + lambda_c * (nac[t] - 1) * ac_marg[t]    // effect of number of attacks to calf 
        ) 
        // decay
        - (1 - a[t]) * decay[period_id[t]] * z[t-1];
    }
  }     
  
  // Design matrix --------------------------------------------------------
  
  for(k in 1:K)
    y_mat_cons_z[, k] = y_mat_cons[, k] .* z[cons_id_prev];
  // (z has to be subsetted as y_mat was because z[t] affects transition
  // from y[t] to y[t+1])
  
  X = append_col(y_mat_cons, y_mat_cons_z);
  
}

model{
  ////  Priors  ////

  //// z params
  // z_0 are given implicit flat priors in [0, 1]
  kappa_raw ~ std_normal();
  delta_m_raw ~ std_normal();
  delta_c_raw ~ std_normal();
  lambda_c_raw ~ std_normal();
  decay_raw ~ std_normal();
  
  //// Transition matrix parameters
  
  to_vector(alpha_intercept_raw) ~ std_normal();
  to_vector(beta_intercept_raw) ~ std_normal();
  
  to_vector(alpha_spline_coef_raw) ~ std_normal();
  to_vector(beta_spline_coef_raw) ~ std_normal();

  to_vector(alpha_raneff_raw) ~ std_normal();
  
  sd_spline_tmat_raw ~ std_normal();
  sd_years_raw ~ std_normal();
 
  ////  Likelihood  //// 
  
  { // Local scope to define variables
    matrix[K*2, K] coef; // stacked alpha and beta matrices for year y
    
    // Loop over years to avoid large matrix multiplications
    for(y in 1:N_years) {
      matrix[year_length[y], K] eta;
      int y_vec_local[year_length[y]];
      
      coef = append_row(alpha[y], beta[y]); 
      //print(coef);
      eta = X[year_start[y] : year_end[y], ] * coef;
      y_vec_local = y_vec_cons[year_start[y] : year_end[y]];
      
      for(n in 1:year_length[y]) 
        y_vec_local[n] ~ categorical_logit(eta[n, ]');
      
    }
  }
  
}
// No generated quantities to avoid C stack size issue

