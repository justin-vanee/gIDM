functions {
  real weighted_nb_loglik(row_vector pi, real mu, real nu) {
    int M = num_elements(pi) - 1;
    real log_lik = negative_infinity();

    for (k in 0:M) {
      if (pi[k + 1] > 0) {
        real log_nb = neg_binomial_2_lpmf(k | mu, nu);
        real log_weighted = log(pi[k + 1]) + log_nb;
        log_lik = log_sum_exp(log_lik, log_weighted);
      }
    }
    return log_lik;
  }
}

data {
  int<lower=1> n;              
  int<lower=1> n_z;            
  int<lower=1> n_N;            
  int<lower=1> T;              
  int<lower=1> p_X;            
  int<lower=1> M;              

  array[n_z] int z;            
  matrix[n_N, M + 1] pi;       
  matrix[n, p_X] X;            
  matrix[n, T] year_mat;       
  vector[n] S;                 
}

parameters {
  vector[p_X] beta;
  vector[T] xi_raw;
  real<lower=0> tau;
  // dispersion parameter 
  real<lower=0> nu_inv;
}

transformed parameters {
  // concentration parameter (nu --> infinity, N ~ Poisson)
  real<lower=0> nu = 1.0 / nu_inv;
  vector[n] log_lambda;
  vector[n] lambda;
  vector[n] prob;
  vector[T] xi = xi_raw - mean(xi_raw);  

  for (i in 1:n) {
    // mean intensity 
    log_lambda[i] = X[i] * beta + year_mat[i] * xi_raw;
    lambda[i] = S[i] * exp(log_lambda[i]);
    // presence/absence
    // prob[i] = 1-pow(nu/(nu+lambda[i]), nu);
    real log_p0 = nu * log(nu) - nu * log(nu + lambda[i]);
    prob[i] = 1 - exp(log_p0); // more stable than the above
  }
  

}

model {
  // Priors
  tau ~ inv_gamma(2, square(1.5));
  beta ~ normal(0, 1.5);
  xi_raw ~ normal(0, sqrt(tau));

  // negative binomial (1 / dispersion) parameter 
  nu_inv ~ cauchy(0, 1);

  // presence/absence likelihood
  for (i in 1:n_z) {
    target += bernoulli_lpmf(z[i] | prob[i]);
  }

  // Weighted NB likelihood for count data
  for (i in 1:n_N) {
    target += weighted_nb_loglik(pi[i], lambda[i + n_z], nu);
  }
}
