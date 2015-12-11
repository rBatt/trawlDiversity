
# =========
# = Setup =
# =========
set.seed(1337)
nT <- 2E2
nX <- 2 # don't change
nU <- 1 # don't change
e_sd <- 1

A <- matrix(c(0.1, 10), nrow=1, ncol=nX)
B <- matrix(c(0.5, 0.2, -0.8, 0.8), nrow=nX, ncol=nX)
# B <- matrix(c(0.5, 0, 0.8, 0), nrow=nX, ncol=nX)
C <- matrix(c(0,0), nrow=nU, ncol=nX)

X <- matrix(NA, nrow=nT, ncol=2, dimnames=list(NULL, c("zoop","phyto")))
U <- matrix(rlnorm(nT*nU), nrow=nT, ncol=nU, dimnames=list(NULL, c("nuts")))
E <- matrix(rnorm(nT*nX) * e_sd, nrow=nT, ncol=nX)


# ============
# = Simulate =
# ============
X[1,] <- c(3,10)

for(i in 2:nT){
	X[i,] <- A + matrix(X[i-1,],nrow=1)%*%B + matrix(U[i-1,],nrow=1)%*%C + E[i,]
}


# =================
# = Plot Simulate =
# =================
par(mfrow=c(2,1))
plot(X[,1], type="l", ylim=c(min(X,0), max(X)))
lines(X[,2], col="blue")
plot(X)


# ==============
# = Fit Models =
# ==============
# ---- MAR ----
ar(X) # can't handle covariates :(

# ---- VAR ----
VAR(X, exogen=U) # seems to miss the constant term (intercept, in the Stan, A in Ives)

# ---- Stan VAR ----
library(rstan)
stan_data <- list(T = nT, M = nX, P = 1, Y = X)

fitted.model <- stan(file = 'pkgBuild/test/VAR_p.stan')
model_w_data <- stan(fit = fitted.model, data=stan_data, iter=100, chains=4)

print(model_w_data)

traceplot(model_w_data, inc_warmup=TRUE, pars=c("A"))