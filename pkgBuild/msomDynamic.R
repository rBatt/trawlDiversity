
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
# ebs.a2 <- ebs.a[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]

# medium data set
# set.seed(1337)
# ind <- mpick(ebs.agg2, p=c(stratum=20, year=30), weight=TRUE, limit=60)
# logic <- expression(
# 	spp%in%spp[ind]
# 	& stratum%in%stratum[ind]
# 	& year%in%year[ind]
# )
# ebs.a1 <- ebs.agg2[eval(logic)][pick(spp, 100, w=FALSE)]
# ebs.a2 <- ebs.a1[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, stemp=stemp, doy=yday(datetime))]

# largest data set
# ebs.a2 <- ebs.agg2[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]

ebs.a1 <- trim_msom("ebs", gridSize=1, plot=TRUE)
ebs.a2 <- ebs.a1[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]


# ======================================
# = Aggregate and Transform Covariates =
# ======================================
# rename columns for shorthand
setnames(ebs.a2, c("btemp"), c("bt"))
# setnames(ebs.a2, c("stemp"), c("st"))
ebs.a2[,yr:=scale(as.integer(year))]

# aggregate and transform (^2) btemp
mk_cov_rv_pow(ebs.a2, "bt", across="K", by=c("stratum","year"), pow=2)
# mk_cov_rv_pow(ebs.a2, "st", across="K", by=c("stratum","year"), pow=2)

# scale and aggregate doy
doy.mu <- ebs.a2[,mean(doy, na.rm=TRUE)]
doy.sd <- ebs.a2[,sd(doy, na.rm=TRUE)]
ebs.a2[,doy:=(doy-doy.mu)/doy.sd]
mk_cov_rv(ebs.a2, "doy", across="K", by=c("stratum","year"))


# ======================
# = Cast Data for Stan =
# ======================
# ---- Get Basic Structure of MSOM Data Input ----
# dynData <- msomData(Data=ebs.a2, n0=2, cov.vars=c(bt="bt",bt2="bt2",yr="yr", st="st", st2="st2", doy="doy"), u.form=~bt+bt2+st+yr, v.form=~year+doy, valueName="abund", cov.by=c("year","stratum"))
setkey(ebs.a2, year, stratum, K, spp)
# dynData <- msomData(Data=ebs.a2, n0=2, cov.vars=c(bt="bt",bt2="bt2",yr="yr", doy="doy"), u.form=~bt+bt2, v.form=~doy+year, valueName="abund", cov.by=c("year","stratum"), u_rv=c("bt","bt2"),  v_rv=c("doy"))
dynData <- msomData(Data=ebs.a2, n0=2, cov.vars=c(bt="bt",bt2="bt2",yr="yr", doy="doy"), u.form=~bt+bt2, v.form=~doy+year, valueName="abund", cov.by=c("year","stratum"),  v_rv=c("doy"))

# ---- Add a counter for nJ (number of sites in each year) ----
dynData$nJ <- apply(dynData$nK, 1, function(x)sum(x>0))

# ---- Aggregate abundance across samples ----
dynData$X <- apply(dynData$X, c(1,2,4), function(x)sum(x))


# =====================
# = Fit Model in Stan =
# =====================
library(rstan)
model_file <- "trawl/trawlDiversity/inst/stan/msomDynamic.stan"

ebs_msom <- stan(
	file=model_file, 
	data=dynData, 
	control=list(stepsize=0.01, adapt_delta=0.95, max_treedepth=15),
	chains=1, iter=50, seed=1337, cores=1, verbose=F
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
	"alpha_mu","alpha_sd","beta_mu[1]","beta_sd[1]","phi_mu","phi_sd","gamma_mu","gamma_sd",
	paste0("alpha[1,",sample(1:dynData$nS,3),"]"),
	paste0("alpha[2,",sample(1:dynData$nS,3),"]"),
	paste0("alpha[3,",sample(1:dynData$nS,3),"]"),
	paste0("beta[1,",sample(1:dynData$nS,3),"]"),
	paste0("beta[2,",sample(1:dynData$nS,3),"]"),
	paste0("phi[",sample(1:dynData$nS,10),"]"),
	paste0("gamma[",sample(1:dynData$nS,10),"]"),
	"logit_psi_t1_raw[2,1]",
	"Omega"
)

print(ebs_msom, inspect_params)


# ===============
# = Diagnostics =
# ===============
# traceplot of chains
rstan::traceplot(ebs_msom, inspect_params, inc_warmup=T)

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
png("~/Desktop/omega_full.png", width=3.5, height=3.5, units="in", res=150)
par(mfrow=c(1,1), mar=c(2.5,2.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=9)
Omega <- rstan::extract(ebs_msom, "Omega")[[1]]
plot(density(Omega, from=0, to=1))
omega_prior_q <- seq(0,1,length.out=length(Omega))
omega_prior <- dbeta(omega_prior_q, 2, 2)
lines(omega_prior_q, omega_prior, col="blue")
dev.off()

sims <- rstan::extract(ebs_msom)
psi_mean <- plogis(apply(sims$logit_psi, 2:4, mean))
theta_mean <- plogis(apply(sims$logit_theta, 2:4, mean))
phi_mean <- plogis(apply(sims$phi, 2, mean))
gamma_mean <- plogis(apply(sims$gamma, 2, mean))

png("~/Desktop/psi_theta_problem_full.png", width=5, height=5, units='in', res=150)
par(mfrow=c(2,2), mar=c(2.5,2.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=9)
hist(psi_mean)
hist(theta_mean)
hist(phi_mean)
hist(gamma_mean)
dev.off()


psi <- data.table(reshape2::melt(psi_mean, varnames=c("year","site","spp"), value.name="psi"), key=c("year","site","spp"))


# ---- Psi and Theta Response Curves ----
png("~/Desktop/psi_theta_responseCurves_full.png", width=3.5, height=6, units="in", res=150)
par(mfrow=c(2,1), mar=c(2.5,2.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=9)
# ---- Species Response Curves via Psi and bt ----
alpha <- apply(sims$alpha, 2:3, mean)
U <- dynData$U[1,,][order(dynData$U[1,,][,2]),]
psi_resp <- plogis(U%*%alpha)

plot(U[,2], psi_resp[,1], ylim=range(psi_resp), type='l', col="gray", xlab="bottom temperature")
for(i in 2:ncol(psi_resp)){
	lines(U[,2], psi_resp[,i], col="gray")
}
par(new=T)
plot(density(dynData$U[,,"bt"]), xaxt="n",yaxt="n", ylab="",xlab="", main="",type="l", col="blue")



# ---- Detectability Response Curve via Theta and doy ----
beta <- apply(sims$beta, 2:3, mean)
V <- apply(dynData$V, 3, function(x)seq(min(x), max(x), length.out=100))
theta_resp <- plogis(V%*%beta)

plot(V[,2], theta_resp[,1], ylim=range(theta_resp), type='l', col="gray", xlab="day of year (scaled)")
for(i in 2:ncol(theta_resp)){
	lines(V[,2], theta_resp[,i], col="gray")
}
par(new=T)
plot(density(dynData$V[,,"doy"]), xaxt="n",yaxt="n", ylab="",xlab="", main="",type="l", col="blue")
dev.off()


# ---- alpha boxplots ----
alpha_123 <- apply(sims$alpha, 2:3, mean)[,1:dynData$nS]

# ---- Boxplots of posterior ALPHA ----
png("~/Desktop/alpha_est_boxplots.png", width=5, height=5, res=200, units="in")
par(mfrow=c(3,1), mar=c(2,2,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=10)
alpha_123_df <- reshape2::melt(sims$alpha[,,1:dynData$nS])
for(i in 1:3){ # 3 alphas
	boxplot(value~Var3+Var2, data=alpha_123_df[alpha_123_df[,"Var2"]==i,], ylab=paste0("alpha[",i,"] posterior"), xlab="speices ID")
}
dev.off()


# ---- beta boxplots ----
beta_123 <- apply(sims$beta, 2:3, mean)[,1:dynData$nS]

# ---- Boxplots of posterior Beta ----
png("~/Desktop/beta_est_boxplots.png", width=5, height=5, res=200, units="in")
par(mfrow=c(3,1), mar=c(2,2,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=10)
beta_123_df <- reshape2::melt(sims$beta[,,1:dynData$nS])
boxplot(value~Var3+Var2, data=beta_123_df[beta_123_df[,"Var2"]==1,], ylab=paste0("beta[",1,"] posterior"), xlab="speices ID")
boxplot(value~Var3+Var2, data=beta_123_df[beta_123_df[,"Var2"]==5,], ylab=paste0("beta[",5,"] posterior"), xlab="speices ID")
boxplot(value~Var3+Var2, data=beta_123_df[beta_123_df[,"Var2"]==31,], ylab=paste0("beta[",31,"] posterior"), xlab="speices ID")
dev.off()


# ---- Boxplots of posterior Phi ----
png("~/Desktop/phi_est_boxplots_blueDotTrue.png", width=5, height=5, res=200, units="in")
par(mfrow=c(2,1), mar=c(2,2,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=10)
phi_123_df <- reshape2::melt(sims$phi[,1:dynData$nS])
boxplot(value~Var2, data=phi_123_df, ylab=paste0("phi posterior"), xlab="speices ID", cex=0.5)
gamma_123_df <- reshape2::melt(sims$gamma[,1:dynData$nS])
boxplot(value~Var2, data=gamma_123_df, ylab=paste0("gamma posterior"), xlab="speices ID", cex=0.5)
dev.off()




# ---- save ----
save.image("trawl/trawlDiversity/pkgBuild/test/msomDynamic_medium.RData")


ebs_msom_sso <- as.shinystan(ebs_msom)
deploy_shinystan(ebs_msom_sso)
