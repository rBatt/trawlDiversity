// ===============
// = Conventions =
// ===============
// Data Structures:
	// Arrays: last index is fastest (most specific, smallest)
	// Matrices: first index (row) is fastest
// Notation:
	// Scalar capital letters are the maximum value for their lowercase counterpart (t = 1, 2, ..., T)
	// Greek letters are parameters
	// State variables and covariates are capital Roman letters (X, Z, Y, U)

// =============
// = Functions =
// =============
functions {
	
	// log probability functions
	
	// lp for species that are observed
	real lp_obs(int x, real logit_psi, real logit_theta) {
		// could probably just drop x; 
		// taken from function that used binomial
		
		// log(p(present) * p(detected))
		return log_inv_logit(logit_psi) + bernoulli_logit_log(x, logit_theta);
	}
	
	// lp for species that aren't observed, but known to exist
	real lp_unobs(real logit_psi, real logit_theta) {
		// log(p(observed)) = log(p(present) * p(detected))
		// log(p(unobs)) = log(p(observed) * 1-p(present))
		
		return log_sum_exp(lp_obs(0, logit_psi, logit_theta), log1m_inv_logit(logit_psi));
	}
	
	// real lp_exists(int x[N], int N, real logit_psi, real logit_theta) {
	// 	// could probably just drop x;
	// 	// taken from function that used binomial
	//
	// 	int is0[N];
	//
	// 	for(i in 1:N){
	// 		is0[i] <- x[i] == 0;
	// 	}
	//
	// 	// log(p(present) * p(detected))
	// 	return log_sum_exp(log_inv_logit(logit_psi) + bernoulli_logit_log(x, logit_theta), log1m_inv_logit(logit_psi)*is0);
	// }
	
	real lp_never_obs(real logit_psi, real logit_theta, real Omega) {
		real lp_unavailable;
		real lp_available;
		lp_unavailable <- bernoulli_log(0, Omega);
		lp_available <- bernoulli_log(1, Omega) + lp_unobs(logit_psi, logit_theta);
		return log_sum_exp(lp_unavailable, lp_available);
	}
	
}


// ========
// = Data =
// ========
data {
    int<lower=1> nT; // number of years
    int<lower=1> Kmax;
    int<lower=1> Jmax; // number of sites
    int<lower=0, upper=Kmax> nK[nT,Jmax]; // number of 'replicates' (hauls)
    int nU; // number of Covariates for psi (presence)
    int nV; // number of Covariates for theta (detection)
    
    matrix[Kmax,nU] U[nT,Jmax]; // covariates for psi (presence)
    matrix[Kmax,nV] V[nT,Jmax]; // covariates for theta (detectability)
    
    int N; // total number of observed species (anywhere, ever)
		int<lower=1> nS; // size of super population, includes unknown species
    int X[nT,Jmax,Kmax,nS]; // species observations
}


// ====================
// = Transformed Data =
// ====================
// transformed data {
// 	int n0;
//
// 	n0 <- nS - N
// }


// ==============
// = Parameters =
// ==============
parameters {	
	real<lower=0, upper=1> Omega[nT]; // average availability
	
	vector[nU] alpha_mu; // hyperparameter mean
	vector<lower=0>[nU] alpha_sd; // hyperparameter sd
	matrix[nU,nS] alpha; // presence coefficient
	
	vector[nV] beta_mu; // hyperparameter mean
	vector<lower=0>[nV] beta_sd; // hyperparameter sd
	matrix[nV,nS] beta; // detection coefficient
	
}


// ==========================
// = Transformed Parameters =
// ==========================
// transformed parameters {
// 	row_vector[nS] logit_psi[nT,Jmax,Kmax]; // presence probability
// 	row_vector[nS] logit_theta[nT,Jmax,Kmax]; // detection probability
//
//
// 	for (t in 1:nT) {
// 		for (j in 1:Jmax) {
// 			// if the K matrix gets too large,
// 			// it is in here that I should add
// 			// a loop, wherein I subset to U[t,j,k], and
// 			// vectorize over spp only
// 			for (k in 1:nK[t,j]) {
// 				logit_psi[t,j,k] <- U[t,j,k]*alpha;
// 				logit_theta[t,j,k] <- V[t,j,k]*beta;
// 			}
// 		}
// 	}
//
// }


// =========
// = Model =
// =========
model {
	// Priors for hyperparameters
	alpha_mu ~ cauchy(0, 1);
	alpha_sd ~ cauchy(0, 2);
	beta_mu ~ cauchy(0, 1);
	beta_sd ~ cauchy(0, 2);
	
	Omega ~ beta(2,2);
	
	// Priors for parameters
	// for(s in 1:nS){
		for (u in 1:nU) {
			alpha[u] ~ normal(alpha_mu[u], alpha_sd[u]);
		}
		for (v in 1:nV) {
			beta[v] ~ normal(beta_mu[v], beta_sd[v]);
		}
	// }
	
	
	// ---- Begin Looping down to Point Observations ----
	for (t in 1:nT) { // loop through years
		increment_log_prob(bernoulli_log(1, Omega[t]) * N * sum(nK[t])); // observed, so available
		for (j in 1:Jmax) { // sites
			for (k in 1:nK[t,j]) {// sampling events
				
				row_vector[nS] logit_psi; // presence probability
				row_vector[nS] logit_theta; // detection probability
				
				logit_psi <- U[t,j,k]*alpha;
				logit_theta <- V[t,j,k]*beta;
				
				// 1)
				for (n in 1:N) {// observed species
					
					if (X[t,j,k,n] > 0) { // if present
						// increment_log_prob(lp_obs(X[t,j,k,n], logit_psi[t,j,k,n], logit_theta[t,j,k,n]));
						increment_log_prob(lp_obs(X[t,j,k,n], logit_psi[n], logit_theta[n]));
					} else { // if absent
						// increment_log_prob(lp_unobs(logit_psi[t,j,k,n], logit_theta[t,j,k,n]));
						increment_log_prob(lp_unobs(logit_psi[n], logit_theta[n]));
					}
				} // end first species loop
				
				// 2)
				for (s in (N+1):nS) { // padded species
					// increment_log_prob(lp_never_obs(logit_psi[t,j,k,s], logit_theta[t,j,k,s], Omega[t]));
					increment_log_prob(lp_never_obs(logit_psi[s], logit_theta[s], Omega[t]));
				} // end second species loop
			} // end sample loop
		} // end site loop
	} // end year loop
	
	
    
} // end model


// ========================
// = Generated Quantities =
// ========================
// generated quantities {
//
// }

