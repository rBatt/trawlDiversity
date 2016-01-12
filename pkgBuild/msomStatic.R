
# ========
# = Load =
# ========
library(devtools)
library("trawlData")
library("trawlDiversity")
library("rstan")


# =========================
# = Timing and Efficiency =
# =========================
stan_out <- rm_out[[1]]
(timing <- cbind(get_elapsed_time(stan_out), rowSums(get_elapsed_time(stan_out))))
(max_time <- max(timing))
(neff_mu <- mean(summary(stan_out)$summary[,"n_eff"]))
(neff_sd <- sd(summary(stan_out)$summary[,"n_eff"]))
max_time/neff_mu



# ==================================
# = Printing the Fitted Parameters =
# ==================================

inspect_params <- c(
	"alpha_mu","alpha_sd","beta_mu","beta_sd",
	"Omega"
)

print(stan_out, inspect_params)


# ===============
# = Diagnostics =
# ===============
# traceplot of chains
rstan::traceplot(stan_out, inspect_params, inc_warmup=F)

# historgram of tree depth -- make sure not hugging max
hist_treedepth <- function(fit) { 
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE) 
  hist(sapply(sampler_params, function(x) c(x[,'treedepth__']))[,1], breaks=0:20, main="", xlab="Treedepth") 
  abline(v=10, col=2, lty=1) 
}
# sapply(stan_out@sim$samples, function(x) attr(x, 'args')$control$max_treedepth)
hist_treedepth(stan_out)

# lp
rstan::traceplot(stan_out, "lp__", window=c(1,50), inc_warmup=T)


# =====================
# = Parameter Summary =
# =====================
# what is the distribution of omega?
par(mfrow=c(2,1), mar=c(2.5,2.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=9)
png("~/Desktop/omega_full.png", width=3.5, height=3.5, units="in", res=150)
Omega <- rstan::extract(stan_out, "Omega")[[1]]
plot(density(Omega, from=0, to=1))
omega_prior_q <- seq(0,1,length.out=length(Omega))
omega_prior <- dbeta(omega_prior_q, 2, 2)
lines(omega_prior_q, omega_prior, col="blue")
dev.off()

sims <- rstan::extract(stan_out)
psi_mean <- plogis(apply(sims$logit_psi, 2:4, mean))
theta_mean <- plogis(apply(sims$logit_theta, 2:4, mean))

png("~/Desktop/psi_theta_problem_full.png", width=3.5, height=6, units='in', res=150)
par(mfrow=c(2,1), mar=c(2.5,2.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=9)
hist(psi_mean)
hist(theta_mean)
dev.off()


psi <- data.table(reshape2::melt(psi_mean, varnames=c("year","site","spp"), value.name="psi"), key=c("year","site","spp"))


# ---- Psi and Theta Response Curves ----
png("~/Desktop/psi_theta_responseCurves_full.png", width=3.5, height=6, units="in", res=150)
par(mfrow=c(2,1), mar=c(2.5,2.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=9)
# ---- Species Response Curves via Psi and bt ----
alpha <- apply(sims$alpha, 2:3, mean)
U <- staticData$U[1,,][order(staticData$U[1,,][,2]),]
psi_resp <- plogis(U%*%alpha)

plot(U[,2], psi_resp[,1], ylim=range(psi_resp), type='l', col="gray", xlab="bottom temperature")
for(i in 2:ncol(psi_resp)){
	lines(U[,2], psi_resp[,i], col="gray")
}
par(new=T)
plot(density(staticData$U[,,"bt"]), xaxt="n",yaxt="n", ylab="",xlab="", main="",type="l", col="blue")



# ---- Detectability Response Curve via Theta and doy ----
beta <- apply(sims$beta, 2:3, mean)
V <- apply(staticData$V, 3, function(x)seq(min(x), max(x), length.out=100))
theta_resp <- plogis(V%*%beta)

plot(V[,2], theta_resp[,1], ylim=range(theta_resp), type='l', col="gray", xlab="day of year (scaled)")
for(i in 2:ncol(theta_resp)){
	lines(V[,2], theta_resp[,i], col="gray")
}
par(new=T)
plot(density(staticData$V[,,"doy"]), xaxt="n",yaxt="n", ylab="",xlab="", main="",type="l", col="blue")
dev.off()

# ---- save ----
save.image("trawl/trawlDiversity/pkgBuild/test/msomStatic_medium.RData")


