data {
  int<lower=0> T; // No. of time periods observed
  int<lower=0> M; // No. of variables in y_t
  int<lower=1> P; // Lag order
  vector[M] Y[T]; // Obs of state variables
  int nU; // No. of covariates
  vector[nU] U[T]; // Covariates
}
transformed data {
  vector[M] Y_obs[T-P];
  for (t in 1:(T-P)) 
    Y_obs[t] <- Y[t + P];
  
}
parameters {
  matrix[M, M] B[P];
  cholesky_factor_corr[M] L_corr_noise;
  vector<lower=0>[M] sd_noise;
  vector[M] A;
  matrix[M, nU] C;
}
transformed parameters {
  matrix[M,M] L_sigma;
  L_sigma <- diag_pre_multiply(sd_noise, L_corr_noise);
}
model {
  vector[M] mus[T-P]; 
  for (t in 1:(T-P)) {
    mus[t] <- A + C * U[t+P];
    for (p in 1:P) 
      mus[t] <- mus[t] + B[p] * Y[t+p-1];
  }
  L_corr_noise ~ lkj_corr_cholesky(2.0);
  sd_noise ~ normal(0,1);
  A ~ normal(0,1);

  for (p in 1:P)
    to_vector(B[p]) ~ normal(0, 1);
  Y_obs ~ multi_normal_cholesky(mus,L_sigma);
}
generated quantities {
  matrix[M,M] Sigma;
  Sigma <- L_sigma * L_sigma';
}