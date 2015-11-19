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
    int<lower=1> T; // number of years
    int<lower=1> Kmax[T];
    int<lower=1> J[T]; // number of sites
    int<lower=1, upper=Kmax> K[J,T]; // number of 'replicates' (hauls)
    int nU_psi; // number of Covariates for psi (presence)
    int nU_theta; // number of Covariates for theta (detection)
    
    vector[nU_psi] [T,S,J]; // covariates for psi (presence)
    vector[nU_theta] U_theta[T,S,J,K]; // covariates for theta (detectability)
    
    int N; // total number of observed species (anywhere, ever)
	int<lower=N> S; // size of super population, included unknown species
    int Y[J,K,S,T]; // species observations
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
	X[J,K,S,T]; // true state of presence
	
	real<lower=0, upper=1> omega; // average availability
	real<lower=0, upper=1> w[S]; // species availability
	
	vector[nU_psi] alpha_mu; // hyperparameter mean
	vector[nU_psi] alpha_sd; // hyperparameter sd
	vector[nU_psi] alpha[S]; // presence coefficient
	
	vector[nU_theta] beta_mu; // hyperparameter mean
	vector[nU_theta] beta_sd; // hyperparameter sd
	vector[nU_theta] beta[T,S]; // detection coefficient
	
}


// ==========================
// = Transformed Parameters =
// ==========================
// transformed parameters {
	array[S,T] logit_psi; // presence probability
	array[S,T] logit_theta; // detection probability
	
	for(t in 1:T){
		for(j in 1:J){
			for(s in 1:S){
				logit_psi[t][s] <- U_psi*alpha[s]
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

