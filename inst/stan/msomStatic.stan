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
	int<lower=N> nS; // size of super population, includes unknown species
    int X[nT,Jmax,Kmax,nS]; // species observations
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
	X[nT,Jmax,Kmax,nS]; // true state of presence; will have to fix some at 0 (ragged J K)
	
	real<lower=0, upper=1> Omega; // average availability
	real<lower=0, upper=1> w[nS]; // species availability
	
	vector[nU] alpha_mu; // hyperparameter mean
	vector[nU] alpha_sd; // hyperparameter sd
	matrix[nU,nS] alpha; // presence coefficient
	
	vector[nV] beta_mu; // hyperparameter mean
	vector[nV] beta_sd; // hyperparameter sd
	matrix[nV,nS] beta; // detection coefficient
	
}


// ==========================
// = Transformed Parameters =
// ==========================
// transformed parameters {
	matrix[Kmax,nS] logit_psi[nT,Jmax]; // presence probability
	matrix[Kmax,nS] logit_theta[nT,Jmax]; // detection probability
	
	int K;
	
	for(t in 1:T){
		for(j in 1:J){
			K <- nK[t,j];
			if(K){ // if K is not 0
				// if the K matrix gets too large, 
				// it is in here that I should add
				// a loop, whereby I subset to U[t,j,k], and
				// vectorize over spp only
				logit_psi[t,j] <- U[t,j]*alpha;
				logit_theta[t,j] <- V[t,j]*beta;
			}
		}
	}
	
// }


// =========
// = Model =
// =========
model {
    
}


// ========================
// = Generated Quantities =
// ========================
generated quantities {
    
}

