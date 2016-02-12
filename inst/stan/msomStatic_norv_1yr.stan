
// =============
// = Functions =
// =============
functions {

	// these functions were part of my failed attempt to vectorize

  // ---- log probability functions ----

  // lp for species that are observed
  vector lp_obs(int[] x, int K,  vector lil_lp, row_vector logit_theta) {
    return lil_lp + binomial_logit_log(x, K, logit_theta);
  }

  // lp for species that aren't observed, but known to exist
  real lp_unobs(int K, vector lil_lp, row_vector logit_theta, vector l1mil_lp) {
    int D;
    D <- cols(logit_theta);
    return log_sum_exp(sum(lp_obs(rep_array(0, D), K, lil_lp, logit_theta)), sum(l1mil_lp));
  }

  // lp that works as either lp_obs or lp_unobs, depending on isUnobs values
  real lp_exists(int[] x, int K, vector lil_lp, row_vector logit_theta, vector isUnobs, vector l1mil_lp) {
    return log_sum_exp(sum(lp_obs(x, K, lil_lp, logit_theta)), sum(l1mil_lp .* isUnobs));
  }

  // lp for species that were never observed
  real lp_never_obs(int K, vector lil_lp, row_vector logit_theta, real Omega, vector l1mil_lp) {
    real lp_unavailable;
    real lp_available;
    int D;
    D <- cols(logit_theta);
    lp_unavailable <- bernoulli_log(0, Omega)*D;
    lp_available <- bernoulli_log(1, Omega)*D + lp_unobs(K, lil_lp, logit_theta, l1mil_lp);
    return log_sum_exp(lp_unavailable, lp_available);
  }


}


// ========
// = Data =
// ========
data {
  int<lower=1> Kmax; // max samples in a site
  int<lower=1> Jmax; // number of sites
  int<lower=0, upper=Kmax> nK[Jmax]; // number of samples (replicates/ hauls)
  int nU; // number of Covariates for psi (presence)
  int nV; // number of Covariates for theta (detection)
  
  int N; // total number of observed species (anywhere, ever)
  int<lower=1> nS; // size of super population, includes unknown species
  vector[nS] isUnobs[Jmax]; // was a species unobserved (across all K) each site?
	
	
  matrix[Jmax, nU] U; // presence covariates
  matrix[Jmax, nV] V; // detection covariates
	
  
  int X[Jmax,nS]; // species abundances
}


// ====================
// = Transformed Data =
// ====================
transformed data {
	int<lower=0> n0;
	
	n0 <- nS - N;
}


// ==============
// = Parameters =
// ==============
parameters { 
  real<lower=0, upper=1> Omega;
  
  vector[nU] alpha_mu; // hyperparameter mean
  vector<lower=0>[nU] alpha_sd; // hyperparameter sd
  matrix[nU,nS] alpha_raw; // non-centered presence coefficient
  
  vector[nV] beta_mu; // hyperparameter mean
  vector<lower=0>[nV] beta_sd; // hyperparameter sd
  matrix[nV,nS] beta_raw; // detection coefficient
  
}


// ==========================
// = Transformed Parameters =
// ==========================
transformed parameters {

  // ---- declare ----
  matrix[nU,nS] alpha; // presence coefficient
  matrix[nV,nS] beta; // detection coefficient
	
	matrix[Jmax, nS] logit_psi; // logit presence probability
	matrix[Jmax, nS] logit_theta; // logit detection probability  
  
  // ---- define ---- 
  // coefficients
  for (u in 1:nU) {
    alpha[u] <- alpha_mu[u] + alpha_sd[u]*alpha_raw[u];
  }
  for (v in 1:nV) {
    beta[v] <- beta_mu[v] + beta_sd[v]*beta_raw[v];
  }
	
	// psi and theta
	logit_psi <- U*alpha;
	logit_theta <- V*beta;
	 
}


// =========
// = Model =
// =========
model {
  
  // ---- Priors for Hyperparameters ----
	alpha_mu ~ cauchy(0, 1);
	alpha_sd ~ cauchy(0, 2);
  beta_mu ~ cauchy(0, 1);
  beta_sd ~ cauchy(0, 2);
  
  // ---- Priors for Parameters ----
  Omega ~ beta(2,2);
  
  for (u in 1:nU) {
    alpha_raw[u] ~ normal(0, 1); // implies alpha ~ normal(alpha_mu, alpha_sd)
  }
  for (v in 1:nV) {
    beta_raw[v] ~ normal(0, 1); // implies beta ~ normal(beta_mu, beta_sd)
  }
  
  
  // ---- Increment Omega for Known Species (Vectorized) ----
	increment_log_prob(bernoulli_log(1, Omega) * N);
	
	// ---- Species that are Known to Exist: Iterative Form ----
	for (n in 1:N) {
		for (j in 1:Jmax) {
			if (X[j,n] > 0) {
				increment_log_prob(log_inv_logit(logit_psi[j,n]) + binomial_logit_log(X[j,n], nK[j], logit_theta[j,n]));
			} else {
				increment_log_prob(log_sum_exp(bernoulli_logit_log(1, logit_psi[j,n]) + binomial_logit_log(0, nK[j], logit_theta[j,n]), bernoulli_logit_log(0, logit_psi[j,n])));
			}
		}
	}

	// ---- Never-Observed Species: Iterative Form ----
	for (s in (N+1):nS) {
		real lp_unavailable; // unavail part of never obs prob
		real lp_available; // available part of never obs prob
		vector[Jmax] lp_available_pt1; // stores 'present but undetected' part of 'never obs but avail' term

		for (j in 1:Jmax) {
			lp_available_pt1[j] <- log_sum_exp(log_inv_logit(logit_psi[j,s]) + binomial_logit_log(0, nK[j], logit_theta[j,s]), log1m_inv_logit(logit_psi[j,s]));
		}
		lp_unavailable <- bernoulli_log(0, Omega);
		lp_available <- bernoulli_log(1, Omega) + sum(lp_available_pt1);
		increment_log_prob(log_sum_exp(lp_unavailable, lp_available));
	}
	
}
	
// ========================
// = Generated Quantities =
// ========================
// generated quantities {
// }
