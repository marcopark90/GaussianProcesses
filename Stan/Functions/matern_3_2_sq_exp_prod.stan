functions {
  matrix matern_3_2_sq_exp_prod(vector x1, vector x2, real alpha, real rho,
                          real a1, real b1, real a2, real b2,
                          real psi1, real psi2, real sigma) {
    int N = num_elements(x1);
    matrix[N, N] K;
    vector[N] wx1;
    vector[N] wx2;
    for(i in 1:N){
    wx1[i] = beta_cdf(x1[i] | a1, b1);
    }
    for(i in 1:N){
    wx2[i] = beta_cdf(x2[i] | a2, b2);
    }
    for (i in 1:(N - 1)) {
      K[i,i] = alpha + sigma;
      for (j in (i + 1):N) {
        real dsq = psi1*(wx1[i]-wx1[j])^2 + psi2*(wx2[i]-wx2[j])^2;
        K[i, j] = alpha*(1+sqrt(3*dsq)/rho)*exp(-sqrt(3*dsq)/rho)*exp(-dsq/(2*rho^2));
        K[j, i] = K[i, j];
      }
    }
    K[N,N] = alpha + sigma;
    K = K + diag_matrix(rep_vector(1e-9, N));
    return K;
  }

}

data {
    int<lower=1> N;
    int<lower=1> N_upper;
    int<lower=1> N_lower;

    vector[N] x1;
    vector[N] x2;
    vector[N] y_input;
    vector[N_upper] y_input_upper;
  }

  parameters {

    real<lower=0> alpha;
    real<lower=0> rho;
    real<lower=0> a1;
    real<lower=0> b1;
    real<lower=0> a2;
    real<lower=0> b2;
    real<lower=0> psi1;
    real<lower=0> psi2;
    real<lower=0> sigma;
    
    cov_matrix[N] KW;

    vector[N_lower] y_lower;
  }

  transformed parameters {
    vector[N] y = append_row(y_input_upper, y_lower);
    matrix[N, N] K = cov_matern_3_2_sq_exp_prod(x1, x2, alpha, rho, a1, b1, a2, b2, psi1, psi2, sigma);
  }

  model {

  // Parameters

  // GP
   KW ~ inv_wishart(N + 2, K);
   y ~ multi_normal(rep_vector(0, N), KW);
  }

  generated quantities {
    vector[N] y_new = multi_normal_rng(y_input, KW);
    vector[N] y_post_pred = multi_normal_rng(rep_vector(0, N), KW);
  }
