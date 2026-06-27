/* 
  In this model there is no year : behav_t-1 interaction, so I add a year
  effect to the complete column in the transition matrix. 
*/

data {
  // Ns and indexes
  int<lower = 1> N;         // full number of obervation intervals
  int<lower = 1> N_cons;    // number of observations with observed current and previous behavior (consecutive)
  int<lower = 2> K;         // number of behaviors
  int<lower = 2> b;         // Basis dimension for tmat splines, WITHOUT INTERCEPT (lowercase means that)
  int<lower = 2> Bz;        // Basis dimension for zpar splines, WITH INTERCEPT (uppercase means that)
  
  int<lower = 1> N_fol;                // number of follows
  int<lower = 1> int_start[N_fol];     // follows begginings
  int<lower = 1> int_end[N_fol];       // follows ends
  
  int<lower = 1> N_years;              // number of years
  
  // index matrix to parameterize transition matrix with zeroes in the diagonal
  int ind_matrix[K-1, K];    
 
  // Years-related data
  matrix[N_years, b] spline_mfit;      // tmat spline design matrix without intercept
  matrix[N_years, Bz] spline_zfit;     // zpar spline design matrix with intercept
  int<lower = 1> year_start[N_years];  
  int<lower = 1> year_end[N_years];
  // year_start and end identify indexes where each year begin or ends 
  // in the subsetted design matrix, which contains only consecutive observed
  // behaviors
  int<lower = 1> year_length[N_years];
  // number of intervals with consecutive observed behaviors by year
  
  int<lower = 1, upper = N_years> year[N]; 
  // year identifier in [1:N_years] to compute z

  // Response variable 
  int<lower = 1, upper = K+1> y_vec[N];   // behavior
  matrix[N, K] y_mat;                     // behavior in matrix form (binary elements)
  int cons_id[N_cons];       // index identifying observations with observed behavior at current and previous interval
  int cons_id_prev[N_cons];  // cons_id - 1 (previous behavior index)  
  
  // attack variables                                  
  vector[N] a;        // attack presence to mother or calf
  vector[N] ac_only;  // attack presence only to calf
  vector[N] am_only;  // attack presence only to mother
  vector[N] ac_marg;  // attack presence to calf independent of what happens to mother
  vector[N] am_marg;  // attack presence to mother independent of what happens to calf
  vector[N] ab;       // attack presence to mother AND calf
  vector[N] nac;      // number of attacks to calf
  
  // Variables to compute likelihood of jumps and UW behaviours
  // Jumps data
  int<lower = 1> N_jumps;              // number of jumps
  int<lower = 1> jump_start[N_jumps];  // first NA observation in jump
  int<lower = 1> jump_end[N_jumps];    // first observed behaviour after jump
  int<lower = 1> jump_year[N_jumps];   // year of jump
  
  // Endings data
  int<lower = 1> N_uw_ends;            // number of UW endigns 
  int<lower = 1> end_start[N_uw_ends]; // their starting positions
  int<lower = 1> end_end[N_uw_ends];   
  int<lower = 1> end_year[N_uw_ends];  
  
  // UW data
  int<lower = 0> uw[N];                // missing behaviour under water
  int half_K;
  int<lower = 1> s_cols[half_K];       // 1, 3, 5: columns corresponding to surface behaviours
  
  
  
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

  // z parameters priors sd
  vector<lower = 0>[Bz] prior_dec_sds;     // decay
  vector<lower = 0>[Bz] prior_inc_sds;     // calf, both, and mother
  real<lower = 0> prior_lambdac_sd;
  
  // Lower limits for decay and increase
  real<lower = 0, upper = 1> L_inc;
  real<lower = 0, upper = 1> L_dec;
}

transformed data {
  int y_vec_cons[N_cons] = y_vec[cons_id];
  matrix[N_cons, K] y_mat_cons = y_mat[cons_id_prev, ];
  // Behavior at t-1 for consecutive observed behaviors (it's cons_id - 1)
}

parameters {
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
  
  //// z dynamics with number of attacks and varying by year with splines
  /* 
     spline coefficients. 4 columns because there are 4 z parameters varying
     by year, in this order:
     decay, increase_m, increase_b, increase_c
  */
  matrix[Bz, 4] zpar_spline_coef_raw;
  real<lower = 0> lambda_c_raw;           // slope for number of attacks to calf  
  // lambda_c does not vary among years
}

transformed parameters {
  
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
  real lambda_c = lambda_c_raw * prior_lambdac_sd;
  matrix[Bz, 4] zpar_spline_coef;
  matrix[N_years, 4] zpar_spline_raw;
  
  // z increase parameters at logit scale
  vector[N_years] kappa_c; // attack to calf
  vector[N_years] kappa_b; // attack to both
  vector[N_years] kappa_m; // attack to mother
  
  // decay
  vector[N_years] decay;       // in [L_dec, 1] 
  
  // increase (p scale)
  vector[N_years] increase_c;  // in [L_inc, 1]
  vector[N_years] increase_b;
  vector[N_years] increase_m;
  
  // derived quantities to access easily:
  vector[N_years] increase_b3; // some attack to mother, 3 to calf
  vector[N_years] increase_c3; // 3 attacks to calf
  
  //// Design matrix elements
  matrix[N_cons, K] y_mat_cons_z; // half desing matrix
  matrix[N_cons, K*2] X; // full design matrix: cbind(y_mat_cons, y_mat_cons_z)
  
  
  //// Transition probabilities parameters ----------------------------------
  
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
  
  // Compute z parameters
  
  /* 
     Spline coefficients for z params have varying priors, because they mean 
     different things. For decay and increase_m, raw splines are logits that will
     be restricted in [lower, 1]. As deltas live in the logit scale and must be 
     >0, we define their splines at log scale:
     deltas = exp(spline_raw_delta).
     
     Columns order [1:4]: decay, increase_c, increase_b, increase_m
     rows order [1:Bz]: quadratic, linear, intercept.
     (priors vary between rows and columns)
  */
  zpar_spline_coef[, 1] = zpar_spline_coef_raw[, 1] .* prior_dec_sds;
  for (j in 2:4) {
    zpar_spline_coef[, j] = zpar_spline_coef_raw[, j] .* prior_inc_sds;
  }
  
  // compute raw splines
  zpar_spline_raw = spline_zfit * zpar_spline_coef;
  
  // compute z parameters with proper restrictions:
  decay = inv_logit(zpar_spline_raw[, 1]) * (1 - L_dec) + L_dec;
  increase_c = inv_logit(zpar_spline_raw[, 2]) * (1 - L_inc) + L_inc;
  increase_b = inv_logit(zpar_spline_raw[, 3]) * (1 - L_inc) + L_inc;
  increase_m = inv_logit(zpar_spline_raw[, 4]) * (1 - L_inc) + L_inc;
  
  // define the kappa at logit scale but metting the constraint above.
  kappa_c = logit(increase_c); 
  kappa_b = logit(increase_b); 
  kappa_m = logit(increase_m); 
  
  // derived quantities to access easily:
  increase_c3 = inv_logit(kappa_c + lambda_c * 2); // 3 attacks to calf
  increase_b3 = inv_logit(kappa_b + lambda_c * 2); // 3 attacks to calf, some to mother
  
  /*
     Initialize z in every follow
     the estimated z_0 is the value before the follow, so attack at t = 1
     helps estimating z at t = 1, the first z value used to compute the
     the likelihood of an observed behavior.
  
     z[t] = f(z[t-1], attacks[t]), and
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
        
        // attack calf (only)
        kappa_c[year[int_start[f]]] * ac_only[int_start[f]] +                          
        
        // attack both
        kappa_b[year[int_start[f]]] * ab[int_start[f]] +
        
        // attack mother (only)
        kappa_m[year[int_start[f]]] * am_only[int_start[f]] +
        
        // if there is more than 1 attack to calf, the inc proportion rises
        lambda_c * (nac[int_start[f]] - 1) * ac_marg[int_start[f]]  
      ) 
      // decay
      - (1 - a[int_start[f]]) * decay[year[int_start[f]]] * z_0[f];
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
        
        // attack calf (only)
        kappa_c[year[t]] * ac_only[t] +                          
        
        // attack both
        kappa_b[year[t]] * ab[t] +
        
        // attack mother (only)
        kappa_m[year[t]] * am_only[t] +
        
        // if there is more than 1 attack to calf, the inc proportion rises
        lambda_c * (nac[t] - 1) * ac_marg[t]
      ) 
      // decay
      - (1 - a[t]) * decay[year[t]] * z[t-1];
    }
  }
  
   
  // Design matrix --------------------------------------------------------
  
  for(k in 1:K)
    y_mat_cons_z[, k] = y_mat_cons[, k] .* z[cons_id_prev];
  // (z has to be subsetted as y_mat was because z[t] affects transition
  // from y[t] to y[t+1])
  
  X = append_col(y_mat_cons, y_mat_cons_z);
  
}

model {
  ////  Priors  ////

  //// z params
  // z_0 are given implicit flat priors in [0, 1]
  to_vector(zpar_spline_coef_raw) ~ std_normal();
  lambda_c_raw ~ std_normal();

  
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
  
    //// Likelihood for non consecutive observed behaviors ("jumps") 
  
  /*
    Jumps are sequences of unobserved behaviours (including UW) bounded by observed
    behaviours within a follow. 
    The jump start is the first NA observation, and
    the jump end is the first observed data after the NA sequence.
    We need to compute the probabilities vector from the first NA (start)
    up to the jump end, which is observed (for t in start:end).
    The last probability vector is the likelihood for the jump-end observation.
    We take into account that UW missing observations can't be surface behaviours 
    (R_S, ST_S or HE_S).
  */
  
  {  // Local scope to define some variables that don't need to be saved
    int start;
    int end;
    int year_local;
    
    row_vector[K] probs_prev;
    row_vector[K] probs_focal;
    row_vector[K] probs_norm; // normalized under water probabilities
    
    matrix[K, K] tmat_logit;
    matrix[K, K] tmat;
    
    // Loop over jumps
    for(j in 1:N_jumps) {
      start = jump_start[j];
      end = jump_end[j];
      year_local = jump_year[j];
      
      // Initialize probabilities row-vector. As the first is observed, 
      // it's a binary vector.
      probs_prev = y_mat[start-1, ]; // start is NA, the observed is start - 1
      
      // Advance probability vector until start + length
      for(t in start:end) {
        
        // Compute transition matrix (it's a function of z)
        tmat_logit = alpha[year_local] + beta[year_local] * z[t - 1]; // z at t-1 affects transition probabilities from t-1 to t.
        // Apply softmax over rows (r)
        for(r in 1:K) {
          tmat[r, ] = softmax(tmat_logit[r, ]')';
        }
    
        // Compute focal probability vector
        probs_focal = probs_prev * tmat;
        
        if (uw[t]) { // uw == 1 when there is an UW undetermined behaviour
          // set surface behaviours probabilities to zero
          probs_focal[s_cols] = rep_row_vector(0, half_K);      // s_cols = c(1, 3, 5)
          // normalize under water behaviours probabilities
          probs_focal = probs_focal / sum(probs_focal);
          
          // add uw behaviours log_likelihood to the log_posterior
          target += log(sum(probs_focal));
        }
        
        // Turn forcal probabilities into prev to advance loop
        probs_prev = probs_focal;
      }
      
      // When we arrive to the jump end, probs_focal carries the likelihood
      // for the ending observed behaviour
      y_vec[end] ~ categorical(probs_focal');
    } // Close loop over jumps
  } // Close local scope
  
  
  //// Likelihood for unobserved endings containing UW
  
  
  /* 
    These are follows that ended in unobserved behaviours containing at least
    one UW in the ending sequence.
    Examples:
      R_UW, R_UW, NA, UW, NA
      R_UW, R_UW, UW, UW, UW
      R_UW, R_UW, NA, NA, UW
    In these cases, if t starts at 1 and it's the first follow in the data 
    end_start = 3 and end_end = 5. Subsequent end_start and end_end match their
    corresponding row in the dataset.
    Only the log_likelihood of UW will be added to the log_posterior.
  */
  
  
  { // Local scope to define some variables that don't need to be saved
    int start;
    int end;
    int year_local;
    
    row_vector[K] probs_prev;
    row_vector[K] probs;
    row_vector[K] probs_norm; // normalized under water probabilities
    
    matrix[K, K] tmat_logit;
    matrix[K, K] tmat;
    
    // Loop over endings
    for(j in 1:N_uw_ends) {
      start = end_start[j];
      end = end_end[j];
      year_local = end_year[j];
      
      // Initialize probabilities row-vector. As the first is observed, 
      // it's a binary vector.
      probs_prev = y_mat[start - 1, ];
      
      for(t in start:end) {
        
        // Compute transition matrix (it's a function of z)
        // (transition from t-1 to t)
        tmat_logit = alpha[year_local] + beta[year_local] * z[t - 1]; // get previous z
        // Apply softmax over rows (r)
        for(r in 1:K) {
          tmat[r, ] = softmax(tmat_logit[r, ]')';
        }
    
        // Compute probability vector for behaviour at time t
        probs = probs_prev * tmat;
        
        if (uw[t]) { // uw == 1 when there is an UW undetermined behaviour
          
          // set surface behaviours probabilities to zero
          probs[s_cols] = rep_row_vector(0, half_K);      // s_cols = c(1, 3, 5)
          // normalize under water behaviours probabilities
          probs_norm = probs / sum(probs);
          // add uw behaviours log_likelihood to the log_posterior
          target += log(sum(probs_norm));
          // move on
          probs = probs_norm;
        }
        
        probs_prev = probs;
      } 
    } // Close loop over ends
  } // Close local scope
  
    
}
// No generated quantities to avoid C stack size issue

