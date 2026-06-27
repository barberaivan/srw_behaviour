
data {
  int N;
  int K;
  int nam[N];
  matrix[N, K] X_behavm;
  real prior_intercepts_sd;
  real prior_sigmas_sd;
}

parameters {
  vector[K] gamma_raw;
  vector<lower = 0>[K] sigma_raw;
  vector[N] error_raw;
}

transformed parameters {
  vector[K] gamma = gamma_raw * prior_intercepts_sd;
  vector<lower = 0>[K] sigma = sigma_raw * prior_sigmas_sd;
}

model {
  gamma_raw ~ std_normal(); 
  sigma_raw ~ std_normal();
  error_raw ~ std_normal(); 
  
  //// Likelihood
  {
    real lambda;

    // Attack number
    for(i in 1:N) {
      lambda = exp(X_behavm[i, ] * gamma + 
                   error_raw[i] * X_behavm[i, ] * sigma);
                   // second line adds a obs-level random effect with 
                   // sd varying by behaviour.
      nam[i] ~ poisson(lambda) T[1, ]; // includes 1
    }   
    
  }
}

generated quantities {
  vector[K] lambda_not_trunc; // lambda from the untruncated distribution
  vector[K] means;            // mean for the truncated distribution

  lambda_not_trunc = exp(gamma + 0.5 * sigma .* sigma);
  for(i in 1:K) 
    means[i] = lambda_not_trunc[i] / (1 - exp(-lambda_not_trunc[i]));
    
}

