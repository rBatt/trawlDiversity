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
  
  // ---- log probability functions ----
  
  // lp for species that are observed
  real lp_obs(int x, real logit_psi, real logit_theta) {
    return log_inv_logit(logit_psi) + bernoulli_logit_log(x, logit_theta);
  }
  
  vector lp_obs2(int[] x, vector lil_lp, row_vector logit_theta) {
    return lil_lp + bernoulli_logit_log(x, logit_theta);
  }
  
  // lp for species that aren't observed, but known to exist
  real lp_unobs(real logit_psi, real logit_theta) {
    return log_sum_exp(lp_obs(0, logit_psi, logit_theta), log1m_inv_logit(logit_psi));
  }
  
  // lp that works as either lp_obs or lp_unobs, depending on isUnobs values
  real lp_exists(int[] x, row_vector logit_psi, row_vector logit_theta, vector isUnobs, vector lil_lp) {
    // FYI: both isUnobs and x are from X[t,j], 
    // isUnobs is !max() swept out across X[t,j,1:K,n] for 1:N, x is just X[t,j,k,1:N]
    // Intention is that isUnobs will be 0 for each x that is observed for any k.
    return log_sum_exp(log_sum_exp(lp_obs2(x, lil_lp, logit_theta)), log_sum_exp(lil_lp .* isUnobs));
  }
  
  // lp for species that were never observed
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
    int<lower=1> Kmax; // max samples in a site-year
    int<lower=1> Jmax; // number of sites
    int<lower=0, upper=Jmax> nJ[nT]; // number of sites in each year
    int<lower=0, upper=Kmax> nK[nT,Jmax]; // number of samples (replicates/ hauls)
    int nU; // number of Covariates for psi (presence)
    int nV; // number of Covariates for theta (detection)
    
    int N; // total number of observed species (anywhere, ever)
    int<lower=1> nS; // size of super population, includes unknown species
    vector[nS] isUnobs[nT,Jmax]; // was a species unobserved (across all K) each site-year?
    
    matrix[Jmax,nU] U[nT]; // covariates for psi (presence)
    matrix[Kmax,nV] V[nT,Jmax]; // covariates for theta (detectability)
    
    int X[nT,Jmax,Kmax,nS]; // species observations
}


// ====================
// = Transformed Data =
// ====================
// transformed data {
// 
// }


// ==============
// = Parameters =
// ==============
parameters { 
	// real Omega_mu; // logit availability hyper mean
	// real<lower=0> Omega_sd; // logit availability hyper sd
	// vector[nT] Omega_logit_raw; // logit availability, non-centered
  real<lower=0, upper=1> Omega[nT]; // average availability
  
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

  matrix[nU,nS] alpha; // presence coefficient
  matrix[nV,nS] beta; // detection coefficient
	// vector[nT] Omega; // availability parameter
	
  for (u in 1:nU) {
    alpha[u] <- alpha_mu[u] + alpha_sd[u]*alpha_raw[u];
  }
  for (v in 1:nV) {
    beta[v] <- beta_mu[v] + beta_sd[v]*beta_raw[v];
  }
	
	// for (t in 1:nT)
	// 	// Omega[t] <- inv_logit(Omega_mu + Omega_logit_raw[t]*Omega_sd); // centered (0,1) availability
	// 	Omega[t] <- inv_logit(Omega_logit_raw[t]*2); // centered (0,1) availability
	
	
}


// =========
// = Model =
// =========
model {
  // Priors for hyperparameters
  alpha_mu ~ cauchy(0, 1);
  alpha_sd ~ cauchy(0, 2);
  beta_mu ~ cauchy(0, 1);
  beta_sd ~ cauchy(0, 2);
	
	// omega hyperparameters (for logit scale)
	// in R, try: hist(plogis(rnorm(100, mean=rnorm(100, 0, 1), sd=abs(rcauchy(100, 0, 1)))))
	// Omega_mu ~ normal(0, 1);
	// Omega_sd ~ cauchy(0, 1);
	
	  
  // Priors for parameters
	// Omega_logit_raw ~ normal(0, 1); // implies Omega ~ inv_logit(normal(Omega_mu, Omega_sd))
	Omega ~ beta(1.5,1.5);
	
  for (u in 1:nU) {
    alpha_raw[u] ~ normal(0, 1); // implies alpha ~ normal(alpha_mu, alpha_sd)
  }
  for (v in 1:nV) {
    beta_raw[v] ~ normal(0, 1);// implies beta ~ normal(beta_mu, beta_sd)
  }
  
  
  // ---- Begin Looping down to Point Observations ----
  for (t in 1:nT) { // loop through years
    increment_log_prob(bernoulli_log(1, Omega[t]) * N * sum(nK[t])); // observed, so available
    
    if(nJ[t]){
      matrix[nJ[t],nS] logit_psi; // presence probability
      logit_psi <- block(U[t], 1, 1, nJ[t], nU)*alpha; // block matrix algebra for speed
    
      for (j in 1:nJ[t]) { // sites
      
        if(nK[t,j]){ // if samples in site
          
          row_vector[N] t_logit_psi;
          vector[N] lil_lp;
          matrix[nK[t,j],nS] logit_theta; // detection probability
          
          t_logit_psi <- sub_row(logit_psi, j, 1, N);
          for (nn in 1:N)
            lil_lp[nn] <- log1m_inv_logit(t_logit_psi[nn]);
          logit_theta <- block(V[t,j], 1, 1, nK[t,j], nV)*beta;
          
          for (k in 1:nK[t,j]) {// sampling events
            
            // lp for species that have been observed at some point
            increment_log_prob(lp_exists(
              segment(X[t,j,k], 1, N), 
              t_logit_psi,
              sub_row(logit_theta, k, 1, N),
              segment(isUnobs[t,j], 1, N),
              lil_lp
            ));
        
            for (s in (N+1):nS) { // padded species
              increment_log_prob(lp_never_obs(logit_psi[j,s], logit_theta[k,s], Omega[t]));
            } // end second species loop
          } // end sample loop
        } // if nK
      } // end site loop
    } // if nJ
  } // end year loop
  
} // end model


// ========================
// = Generated Quantities =
// ========================
// generated quantities {
//
// }

