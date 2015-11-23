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
// functions {
//
// }


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
    
    // int N; // total number of observed species (anywhere, ever)
	int<lower=1> nS; // size of super population, includes unknown species
    vector[nS] X[nT,Jmax,Kmax]; // species observations
}


// ====================
// = Transformed Data =
// ====================
// transformed data {
    
// }


// ==============
// = Parameters =
// ==============
parameters {
	vector[nS] Z[nT,Jmax,Kmax]; // true state of presence; will have to fix some at 0 (ragged J K)
	
	real<lower=0, upper=1> Omega[nT]; // average availability
	real<lower=0, upper=1> vector[nS] w[nT]; // species availability
	
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
	
	int K;
	
	for(t in 1:T){
		for(j in 1:J){
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
	
	// for(s in 1:nS){
	// 	w[s] ~ bernoulli(Omega);
	// }
	// for(t in 1:nt){
	// 	w[t] ~ bernoulli(Omega[t]);
	// }
	
	
	// loop through years
	for(t in 1:nT){
		
		// draw parameter for availability (1 or 0)
		w[t] ~ bernoulli(Omega[t]);
		
		// loop through sites
		for(j in 1:Jmax){
			
			// Define variable for number of sampling events
			// (at this site, during this year)
			real tK[1];
			tK <- nK[t,j];
			
			// if >0 sampling events, proceed
			if (tK){
				
				// loop through sampling events
				for(k in 1:tK){
					
					// true presence
					Z[t,j,k] ~ bernoulli_logit(psi_logit[t,j,k] * w[t]);
					
					// observed presence; likelihood
					X[t,j,k] ~ bernoulli_logit_log(theta_logit[t,j,k] * Z[t,j,k]);
					
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

