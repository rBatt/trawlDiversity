
# ========
# = Load =
# ========
library(devtools)
library("trawlData")
library("trawlDiversity")


# ======================
# = Subset Data to Use =
# ======================
# smallest data set
ebs.a2 <- ebs.a[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]

# largest data set
# ebs.a2 <- ebs.agg2[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]

# medium data set
# set.seed(1337)
# ind <- mpick(ebs.agg2, p=c(stratum=10, year=5), weight=TRUE, limit=60)
# logic <- expression(
# 	spp%in%spp[ind]
# 	& stratum%in%stratum[ind]
# 	& year%in%year[ind]
# )
# ebs.a1 <- ebs.agg2[eval(logic)][pick(spp, 25, w=TRUE)]
# ebs.a2 <- ebs.a1[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]


# ======================================
# = Aggregate and Transform Covariates =
# ======================================
# rename columns for shorthand
setnames(ebs.a2, c("btemp"), c("bt"))
ebs.a2[,yr:=scale(as.integer(year))]

# aggregate and transform (^2) btemp
mk_cov_rv_pow(ebs.a2, "bt", across="K", by=c("stratum","year"), pow=2)

# scale and aggregate doy
doy.mu <- ebs.a2[,mean(doy, na.rm=TRUE)]
doy.sd <- ebs.a2[,sd(doy, na.rm=TRUE)]
ebs.a2[,doy:=(doy-doy.mu)/doy.sd]
mk_cov_rv(ebs.a2, "doy", across="K", by=c("stratum","year"))


# ======================
# = Cast Data for Stan =
# ======================
# ---- Get Basic Structure of MSOM Data Input ----
staticData <- msomData(Data=ebs.a2, n0=3, cov.vars=c(bt="bt", bt2="bt2",doy="doy",yr="yr"), u.form=~bt+bt2, v.form=~doy, valueName="abund", cov.by=c("year","stratum"), u_rv=c("bt","bt2"), v_rv=c("doy"))

# ---- Add a counter for nJ (number of sites in each year) ----
staticData$nJ <- apply(staticData$nK, 1, function(x)sum(x>0))

# ---- Aggregate abundance across samples ----
staticData$X <- apply(staticData$X, c(1,2,4), function(x)sum(x))


# =====================
# = Fit Model in Stan =
# =====================
library(rstan)
model_file <- "trawl/trawlDiversity/inst/stan/msomStatic.stan"

ebs_msom <- stan(
	file=model_file, 
	data=staticData, 
	control=list(stepsize=0.01, adapt_delta=0.95, max_treedepth=15),
	chains=4, iter=50, refresh=1, seed=1337, cores=4, verbose=F
)


# =========================
# = Timing and Efficiency =
# =========================
(timing <- cbind(get_elapsed_time(ebs_msom), rowSums(get_elapsed_time(ebs_msom))))
(max_time <- max(timing))
(neff_mu <- mean(summary(ebs_msom)$summary[,"n_eff"]))
(neff_sd <- sd(summary(ebs_msom)$summary[,"n_eff"]))
max_time/neff_mu



# ==================================
# = Printing the Fitted Parameters =
# ==================================

inspect_params <- c(
	"alpha_mu","alpha_sd","beta_mu","beta_sd",
	"Omega"
)

print(ebs_msom, inspect_params)


# ===============
# = Diagnostics =
# ===============
# traceplot of chains
rstan::traceplot(ebs_msom, inspect_params, inc_warmup=F)

# historgram of tree depth -- make sure not hugging max
hist_treedepth <- function(fit) { 
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE) 
  hist(sapply(sampler_params, function(x) c(x[,'treedepth__']))[,1], breaks=0:20, main="", xlab="Treedepth") 
  abline(v=10, col=2, lty=1) 
}
# sapply(ebs_msom@sim$samples, function(x) attr(x, 'args')$control$max_treedepth)
hist_treedepth(ebs_msom)

# lp
rstan::traceplot(ebs_msom, "lp__", window=c(1,50), inc_warmup=T)


# =====================
# = Parameter Summary =
# =====================
# what is the distribution of omega?
par(mfrow=c(2,1), mar=c(2.5,2.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=9)
png("~/Desktop/omega_full.png", width=3.5, height=3.5, units="in", res=150)
Omega <- rstan::extract(ebs_msom, "Omega")[[1]]
plot(density(Omega, from=0, to=1))
omega_prior_q <- seq(0,1,length.out=length(Omega))
omega_prior <- dbeta(omega_prior_q, 2, 2)
lines(omega_prior_q, omega_prior, col="blue")
dev.off()

sims <- rstan::extract(ebs_msom)
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


