data {
	int<lower=0> T;         // No. of time periods observed
	int<lower=0> M;         // No. of variables in y_t
	int<lower=1> P;         // Lag order
	vector[M] Y[T];    
}
transformed data {
	vector[M] Y_obs[T-P];
	for (t in 1:(T-P)) 
		Y_obs[t] <- Y[t + P];
	
}
parameters {
	matrix[M, M] A[P];
	cholesky_factor_corr[M] L_corr_noise;
	vector<lower=0>[M] sd_noise;
	vector[M] intercept;
}
transformed parameters {
	matrix[M,M] L_sigma;
  L_sigma <- diag_pre_multiply(sd_noise, L_corr_noise);
}
model {
	vector[M] mus[T-P]; 
	for (t in 1:(T-P)) {
		mus[t] <- intercept;
		for (p in 1:P) 
			mus[t] <- mus[t] + A[p] * Y[t+p-1];
	}
	L_corr_noise ~ lkj_corr_cholesky(2.0);
	sd_noise ~ normal(0,1);
	intercept ~ normal(0,1);

	for (p in 1:P)
		to_vector(A[p]) ~ normal(0, 1);
	Y_obs ~ multi_normal_cholesky(mus,L_sigma);
}
generated quantities {
	matrix[M,M] Sigma;
	Sigma <- L_sigma * L_sigma';
}