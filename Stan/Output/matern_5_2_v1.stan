functions {
matrix cov_matern_5_2(vector x1, vector x2, real alpha, real rho,
                      real a1, real b1, real a2, real b2,
                      real psi1, real psi2, real sigma) {
    int N = num_elements(x1);
    matrix[N, N] K;
    vector[N] wx1;
    vector[N] wx2;
    for(i in 1:N){
    wx1[i] = beta_cdf(x1[i], a1, b1);
    }
    for(i in 1:N){
    wx2[i] = beta_cdf(x2[i], a2, b2);
    }
    for (i in 1:(N - 1)) {
      K[i,i] = alpha + sigma;
      for (j in (i + 1):N) {
        real dsq = psi1*(wx1[i]-wx1[j])^2 + psi2*(wx2[i]-wx2[j])^2;
        K[i, j] = alpha*(1+sqrt(5*dsq)/rho+5*dsq/(3*rho^2))*exp(-sqrt(5*dsq)/rho);
        K[j, i] = K[i, j];
      }
    }
    K[N,N] = alpha + sigma;
    K = K + diag_matrix(rep_vector(1e-9, N));
    return cholesky_decompose(K);
  }

}

data {
  int<lower=1> N; // sample size
  int<lower=1> N1; // training sample size
  int<lower=1> N2; // prediction sample size

  vector[N] x1; // development lag
  vector[N] x2; // accident year

  vector[N1] y_input; // target (standardized losses)
}

transformed data {
  vector[N] mu; // mean vector for GP prior
  mu = rep_vector(0, N);
}

parameters{
  real<lower=0> alpha; // covariance paramters
  real<lower=0> rho; // covariance paramters
  real<lower=0> a1; // input warping parameters
  real<lower=0> b1; // input warping parameters
  real<lower=0> a2; // input warping parameters
  real<lower=0> b2; // input warping parameters
  real<lower=0> psi1; // bandwith parameters
  real<lower=0> psi2; // bandwith parameters
  real<lower=0> sigma; // noise parameters
  
  vector[N2] y_missing; // prediction vector
}

transformed parameters {
  vector[N] y;

  for (n in 1:N1) y[n] = y_input[n];
  for (n in 1:N2) y[N1 + n] = y_missing[n];

}

model{
  matrix[N,N] K_exp = cov_matern_5_2(x1, x2, alpha, rho, a1, b1, a2, b2, psi1, psi2, sigma) + diag_matrix(rep_vector(1e-9, N));
  
  alpha ~ student_t(nu, mu, sigma);
  sigma ~ student_t(nu, mu, sigma);
  rho ~ student_t(nu, mu, sigma);

  a1 ~ lognormal(mu, sigma);
  b1 ~ lognormal(mu, sigma);
  a2 ~ lognormal(mu, sigma);
  b2 ~ lognormal(mu, sigma);

  psi1 ~ gamma(alpha, beta);
  psi2 ~ gamma(alpha, beta);
  
  // GP
  y ~ multi_normal_cholesky(mu, K_exp);
}
