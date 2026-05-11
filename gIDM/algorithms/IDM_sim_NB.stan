data {
  
  // ---- Dimensions ----
  int<lower=1> n;                     // number of sites
  int<lower=1> p;                     // number of covariates 
  int<lower=1> n_z;                   // number of presence/absence sites
  int<lower=1> n_N;                   // number of abundance sites

  // ---- Matrices ----
  matrix[n, n] D;                     // distance matrix 
  matrix[n, p] X;                     // covariate matrix 

  // ---- Data ----
  array[n_z]int<lower=0, upper=1> z;  // presence/absence observations 
  array[n_N]int<lower=0> N;           // abundance observations 
  
  // ---- Fixed Parameters ----
  real<lower=0> phi;                  // spatial range (treating as known to facilitate comparison)

}

transformed data {
  
  // geostatistical covariance 
  matrix[n, n] R = exp(-D / phi);
  matrix[n, n] L_R = cholesky_decompose(R);

}


parameters {
  
  // regression coefficients 
  vector[p] beta;
  // spatial parameter
  real<lower=0> sigma2;
  // dispersion parameter 
  real<lower=0> nu_inv;
  // linear predictor 
  vector[n] log_lambda;
  
}

transformed parameters {
  
  // concentration parameter (nu --> infinity, N ~ Poisson)
  real<lower=0> nu = 1.0 / nu_inv;
  
  // intensity 
  vector[n_N+n_z] lambda;
  
  // probability of occupancy 
  vector[n_N+n_z] pi;
  
  for (i in 1:(n_N+n_z)) {
    
    // abundance 
    lambda[i] = exp(log_lambda[i]);
    
    // presence/absence
    // pi[i] = 1-pow(nu/(nu+lambda[i]), nu);
    real log_p0 = nu * log(nu) - nu * log(nu + lambda[i]);
    pi[i] = 1 - exp(log_p0); // more stable than the above
    
  }

}

model {
  
  //
  // Priors 
  //
  
  // variance parameters 
  sigma2 ~ inv_gamma(2, 1.5^2);

  // regression coefficients 
  beta ~ normal(0, 1.5);

  // negative binomial (1 / dispersion) parameter 
  nu_inv ~ cauchy(0, 1);

  //
  // Process
  //
  
  log_lambda ~ multi_normal_cholesky(X * beta, sqrt(sigma2) * L_R);
  
  //
  // Likelihood
  // 
  
  // abundance likelihood (we treat first n_z observations as PA data)
  for (i in 1:n_N) {
    target += neg_binomial_2_lpmf(N[i] | lambda[i+n_z], nu);
  }
  
  // presence/absence likelihood
  for (i in 1:n_z) {
    target += bernoulli_lpmf(z[i] | pi[i]);
  }
  

}

generated quantities {
  
  // spatial random effect (or residual)
  vector[n] eta;
  for(i in 1:n){
    eta[i] = log_lambda[i] - X[i,] * beta;
  }
  
}
