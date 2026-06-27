data {
  int N;
  int Nyear;
  int Nfol;
  int follow[N];
  int year[N];
  int N_with[N];
  int N_int[N];
}

parameters {
  real alpha_raw;
  vector[Nfol] error_fol;
  vector[Nyear] error_year;
  real<lower = 0> sd_fol_raw;
  real<lower = 0> sd_year_raw;
  real<lower = 0> nu_fol_raw;
  real<lower = 0> nu_year_raw;
}

transformed parameters {
  vector[N] p_obs;
  
  real alpha = alpha_raw * 3;
  
  real<lower = 0> sd_fol = sd_fol_raw * 3;
  real<lower = 0> sd_year = sd_year_raw * 3;
  real<lower = 0> nu_fol = nu_fol_raw * 3;
  real<lower = 0> nu_year = nu_year_raw * 3;
  
  real marginal_sd = sqrt(sd_fol ^ 2 + sd_year ^ 2);
  for(n in 1:N) 
    p_obs[n] = inv_logit(alpha + error_fol[follow[n]] + error_year[year[n]]);
  
}

model {
  alpha_raw ~ std_normal();
  error_year ~ student_t(nu_year, 0, sd_year);
  error_fol ~ student_t(nu_fol, 0, sd_fol);
  sd_year_raw ~ std_normal();
  sd_fol_raw ~ std_normal();
  nu_year_raw ~ std_normal();
  nu_fol_raw ~ std_normal();
  
  for(n in 1:N)
    N_with[n] ~ binomial(N_int[n], p_obs[n]);
    
}

