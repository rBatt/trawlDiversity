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
    int<lower=1, upper=Kmax> nK[Jmax,nT]; // number of 'replicates' (hauls)
    int nU; // number of Covariates for psi (presence)
    int nV; // number of Covariates for theta (detection)
    
    matrix[Kmax, nU] U[nT,Jmax]; // covariates for psi (presence)
    matrix[Kmax, nV] V[nT,Jmax]; // covariates for theta (detectability)
    
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
	vector[nS] Z[nT,Jmax,Kmax]; // true state of presence; will have to fix some at 0 (ragged J K)
	
	real<lower=0, upper=1> Omega[nT]; // average availability
	// real<lower=0, upper=1> w[nT,nS]; // species availability
	
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
transformed parameters {
	matrix[Kmax,nS] logit_psi[nT,Jmax]; // presence probability
	matrix[Kmax,nS] logit_theta[nT,Jmax]; // detection probability
	
	
	for(t in 1:nT){
		for(j in 1:Jmax){
			int K;
			K <- nK[t,j];
			if(K){ // if K is not 0
				// if the K matrix gets too large, 
				// it is in here that I should add
				// a loop, wherein I subset to U[t,j,k], and
				// vectorize over spp only
				logit_psi[t,j] <- U[t,j]*alpha;
				logit_theta[t,j] <- V[t,j]*beta;
			}
		}
	}
	
}


// =========
// = Model =
// =========
model {
	// Priors for hyperparameters
	alpha_mu ~ cauchy(0, 2.5);
	alpha_sd ~ cauchy(0, 2.5);
	beta_mu ~ cauchy(0, 2.5);
	beta_sd ~ cauchy(0, 2.5);
	
	Omega ~ beta(2,2);
	
	
	// Priors for parameters
	for(s in 1:nS){
		for(u in 1:nU){
			alpha[u,s] ~ normal(alpha_mu[u], alpha_sd[u]);
		}
		for(v in 1:nV){
			beta[v,s] ~ normal(beta_mu[v], beta_sd[v]);
		}
	}
	
	
	// ============================================
	// = Begin Looping down to Point Observations =
	// ============================================
	for(t in 1:nT){ // loop through years
		for(j in 1:Jmax){ // loop through sites
			
			// Define variable for number of sampling events
			// (at this site, during this year)
			int tK;
			tK <- nK[t,j];
			
			if (tK){ // if >0 sampling events, proceed
				for(k in 1:tK){// loop through sampling events
					
					
					// =====================================
					// = Two Species Loops for Likelihoods =
					// =====================================
					
					// 1)
					for(n in 1:N){// loop through observed species
						1 ~ bernoulli(Omega[t]); // observed, so available
						if(X[t,j,k,n] > 0){ // if present
							// increment likelihood
							increment_log_prob(
								lp_obs(X[t,j,k,n], logit_psi[t,j,k,n], logit_theta[t,j,k,n])
							);
						}else{ // if absent
							// increment likelihood
							increment_log_prob(
								lp_unobs(logit_psi[t,j,k,n], logit_theta[t,j,k,n])
							);
						}
					} // end first species loop
					
					// 2)
					for(s in (N+1):nS){ // loop through padded species
						// increment likelihood
						increment_log_prob(
							lp_never_obs(logit_psi[t,j,k,s], logit_theta[t,j,k,s], Omega[t])
						);
					} // end second species loop
					
				} // end sample loop
			} // if 0 samples, do nothing
		} // end site loop
	} // end year loop
	
	
    
} // end model


// ========================
// = Generated Quantities =
// ========================
// generated quantities {
//
// }

