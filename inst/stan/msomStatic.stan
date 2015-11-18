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
    
    array[J,T,nU_psi] U_psi; // covariates for psi (presence)
    array[J,K,T,nU_p] U_theta; // covariates for theta (detectability)
    
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
	X[J,K,S,T];
	
	real<lower=0, upper=1> omega; // average availability
	real<lower=0, upper=1> w[S]; // species availability
	
	
	
    array[nU_psi,S] logit_psi; // presence parameters
	array[nU_p,S,T] logit_theta; // detectability parameters
	
}


// ==========================
// = Transformed Parameters =
// ==========================
// transformed parameters {
//
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

