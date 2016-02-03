
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
  int<lower=1> nT; // number of years
  int<lower=1> Kmax; // max samples in a site-year
  int<lower=1> Jmax; // number of sites
  int<lower=0, upper=Jmax> nJ[nT]; // number of sites in each year
  int<lower=0, upper=Kmax> nK[nT,Jmax]; // number of samples (replicates/ hauls)
  int nU; // number of Covariates for psi (presence)
  int nV; // number of Covariates for theta (detection)
  
  int N; // total number of observed species (anywhere, ever)
  int<lower=1> nS; // size of super population, includes unknown species
  vector[nS] isUnobs[nT,Jmax]; // was a species unobserved (across all K) each site-year?
	
	
  matrix[Jmax, nU] U[nT]; // presence covariates
  matrix[Jmax, nV] V[nT]; // detection covariates
	
  
  int X[nT,Jmax,nS]; // species abundances
}


// ====================
// = Transformed Data =
// ====================
transformed data {
	int<lower=0, upper=(Jmax*nT)> nJ_sum;
	int<lower=0> n0;
	
	n0 <- nS - N;
	nJ_sum <- sum(nJ);
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
	
	matrix[Jmax, nS] logit_psi[nT]; // logit presence probability
	matrix[Jmax, nS] logit_theta[nT]; // logit detection probability  
  
  // ---- define ---- 
  // coefficients
  for (u in 1:nU) {
    alpha[u] <- alpha_mu[u] + alpha_sd[u]*alpha_raw[u];
  }
  for (v in 1:nV) {
    beta[v] <- beta_mu[v] + beta_sd[v]*beta_raw[v];
  }
	
	// psi and theta
	for (t in 1:nT) {
		logit_psi[t] <- U[t]*alpha;
		logit_theta[t] <- V[t]*beta;
	}
	 
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
		for (t in 1:nT) {
			for (j in 1:Jmax) {
				if (nK[t,j] > 0) {
					if (X[t,j,n] > 0) {
						increment_log_prob(log_inv_logit(logit_psi[t][j,n]) + binomial_logit_log(X[t,j,n], nK[t,j], logit_theta[t][j,n]));
					} else {
						increment_log_prob(log_sum_exp(bernoulli_logit_log(1, logit_psi[t][j,n]) + binomial_logit_log(0, nK[t,j], logit_theta[t][j,n]), bernoulli_logit_log(0, logit_psi[t][j,n])));
					}
				}
			}
		}
	}

	// ---- Never-Observed Species: Iterative Form ----
	for (s in (N+1):nS) {
		real lp_unavailable; // unavail part of never obs prob
		real lp_available; // available part of never obs prob
		vector[nJ_sum] lp_available_pt1; // stores 'present but undetected' part of 'never obs but avail' term
		int pos;

		pos <- 1;

		for (t in 1:nT) {
			for (j in 1:Jmax) {
				if (nK[t,j] > 0) {
					lp_available_pt1[pos] <- log_sum_exp(log_inv_logit(logit_psi[t][j,s]) + binomial_logit_log(0, nK[t,j], logit_theta[t][j,s]), log1m_inv_logit(logit_psi[t][j,s]));
					pos <- pos + 1;
				}
			}
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
