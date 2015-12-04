
library("simOcc")
library("trawlDiversity")
library("trawlData")
library("rbLib")
library("rstan")

set.seed(1337)
sim_out <- sim_occ(ns=8, grid.w=3, grid.h=5, grid.t=3, n0s=0, h.slope=5, detect.mus=c(0,1), alpha_mu=c(0.5, 0.75, 0), alpha_sd=c(0.5, 0.01, 0), format.msom="jags")
big.out.obs <- sim_out[["big.out.obs"]]
dims <- attr(big.out.obs[[1]], "dims")
grid.w <- dims["grid.w"]
grid.h <- dims["grid.h"]
grid.t <- dims["grid.t"]
ns <- dims["ns"]
n.obs.reps <- length(big.out.obs)

par(mfrow=c(2,1))
plot_simResp(big.out.obs[[1]])
hist(values(attr(big.out.obs[[1]], "grid.X")))

# ========================================
# = Analyze simData with Stan msomStatic =
# ========================================
spp2msom_2dt <- function(big.simDat, simCov){
	test.spp <- data.table(reshape2:::melt.list(big.simDat))

	test.spp[,c("year","n.obs.reps"):=list(as.numeric(gsub("year([0-9]{1,2})\\.[0-9]{1,2}","\\1",L1)), as.numeric(gsub("year([0-9]{1,2})\\.([0-9]{1,2})","\\2",L1)))]
	test.spp[,L1:=NULL]
	setnames(test.spp, "value", "abund")
	setkey(test.spp, year,stratum,K)

	test.cov1 <- data.table(reshape2:::melt.list(simCov))
	setnames(test.cov1, c("value","L1"), c("bt","year"))
	test.cov1[,stratum:=test.spp[,unique(stratum)]]
	setcolorder(test.cov1, c('year','stratum','bt'))
	setkey(test.cov1, year, stratum)

	test <- test.spp[test.cov1]
	
	setkey(test, year, stratum, K, spp)
	
	return(test)
}
sim_dt <- spp2msom_2dt(sim_out[["big.simDat"]], simCov=sim_out[["formatted"]]$simCov)

# aggregate and transform (^2) btemp
# mk_cov_rv_pow(sim_dt, "bt", across="K", by=c("stratum","year"), pow=2)
mk_cov_rv(sim_dt, "bt", across="K", by=c("stratum","year"))

# scale and aggregate doy
sim_dt[,yr:=as.numeric(year)]
sim_dt[,obs_rep:=as.numeric(n.obs.reps)]
# mk_cov_rv(sim_dt, "yr", across="K", by=c("stratum","year"))

# sim_dt <- sim_dt[n.obs.reps==1]
stopifnot(sim_dt[,all(n.obs.reps==1)])

staticData <- msomData(Data=sim_dt, n0=2, cov.vars=c(bt="bt",yr="yr"), u.form=~bt, v.form=~yr, valueName="abund", cov.by=c("year","stratum"))

# ---- Get Basic Structure of MSOM Data Input ----
# staticData <- msomData(Data=sim_dt, n0=2, cov.vars=c(bt="bt",yr="yr"), u.form=~bt, v.form=~yr, valueName="abund", cov.by=c("year","stratum"))
#
# staticData_sd <- msomData(Data=sim_dt, n0=1, cov.vars=c(bt_sd="bt_sd",yr_sd="yr_sd"), u.form=~bt_sd-1, v.form=~yr_sd-1, valueName="abund", cov.by=c("year","stratum"))[c("U","V")]
#
#
# # ---- Do the Splitting Function ----
# for(UV in c("U","V")){
# 	for(type in c("constant","mu","sd")){
# 		switch(type,
# 			constant={staticData[[paste0(UV,"_c")]] <- getCovType(staticData, staticData_sd, UV,type)},
# 			{staticData[[paste(UV,type,sep="_")]] <- getCovType(staticData, staticData_sd, UV,type)}
# 		)
# 	}
# }
#
# # ---- Add Sizes of Constant and RV U/V Arrays ----
# staticData$nU_rv <- dim(staticData$U_mu)[3]
# staticData$nV_rv <- dim(staticData$V_mu)[3]
# staticData$nU_c <- dim(staticData$U_c)[3]
# staticData$nV_c <- dim(staticData$V_c)[3]
#
# stopifnot(staticData$nV == staticData$nV_c + staticData$nV_rv)
# stopifnot(staticData$nU == staticData$nU_c + staticData$nU_rv)

# ---- Add a counter for nJ (number of sites in each year) ----
staticData$nJ <- apply(staticData$nK, 1, function(x)sum(x>0))

# ---- Aggregate abundance across samples ----
staticData$X <- apply(staticData$X, c(1,2,4), function(x)sum(x))


# ---- Check Data Format and Simpler Reg Results ----
# png("~/Desktop/issue103_diagnosis.png", width=4, height=7, res=150, units="in")
par(mfrow=c(4,1), ps=10, mar=c(3,3,0.1,0.1), mgp=c(1.5,0.1,0), tcl=-0.1)
sim_dt[,plot(stratum, bt, xlab="stratum (from sim_dt, not inferred)", ylab="btemp (from sim_dt)")] # useful in simulated examples
sim_dt[,plot(bt, abund, ylab="presence/ absence from sim_dt", xlab="btemp from sim_dt")]
plot(staticData$U[1,,"bt"], ylab="btemp from U", xlab="stratum (as inferred from column order in U)\nShould be same as first plot, if not, there will be mismatch")

plot(staticData$U[,,"bt"], pmin(staticData$X[,,1],1), pch=21, ylab="presence/ absence from X", xlab="btemp from U\nIf not same as 2nd plot, X is scrambled relative to U")
for(i in 2:ns){
	points(staticData$U[,,"bt"], pmin(staticData$X[,,i],1), pch=21)
}
# dev.off()


naive <- glm(c(pmin(staticData$X[,,], 1)) ~ rep(c(staticData$U[,,2]), dim(staticData$X)[3]), family="binomial")
summary(naive)

plot(rep(c(staticData$U[,,2]), dim(staticData$X)[3]), c(pmin(staticData$X[,,], 1)))
plot(fitted(naive), type="l")

sim_dt2 <- sim_dt
sim_dt2 <- sim_dt2[,list(abund=max(abund), bt=mean(bt)) ,by=c("stratum","spp","year")]
summary(glmer(abund~bt + (1|spp) + (spp-1|spp), data=sim_dt2, family="binomial"))
summary(glm(abund~bt, data=sim_dt2, family="binomial"))


# =====================
# = Fit Model in Stan =
# =====================
model_file <- "trawl/trawlDiversity/inst/stan/msomStatic.stan"

sim_msom <- rstan::stan(
	file=model_file, 
	data=staticData, 
	control=list(stepsize=0.01, adapt_delta=0.95, max_treedepth=15),
	chains=4, iter=100, seed=1337, cores=4, verbose=F
)

# ==================================
# = Printing the Fitted Parameters =
# ==================================

inspect_params <- c(
	"alpha_mu","alpha_sd","beta_mu","beta_sd",# "alpha",
	"Omega"
)
sims <- rstan::extract(sim_msom)

print(sim_msom, inspect_params)
attr(big.out.obs[[1]], "a3")
apply(sims$alpha, 2:3, mean)[,1:ns]


# ===============
# = Diagnostics =
# ===============
# traceplot of chains
rstan::traceplot(sim_msom, inspect_params, inc_warmup=F)

# historgram of tree depth -- make sure not hugging max
hist_treedepth <- function(fit) { 
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE) 
  hist(sapply(sampler_params, function(x) c(x[,'treedepth__']))[,1], breaks=0:20, main="", xlab="Treedepth") 
  abline(v=10, col=2, lty=1) 
}
# sapply(ebs_msom@sim$samples, function(x) attr(x, 'args')$control$max_treedepth)
hist_treedepth(sim_msom)

# lp
rstan::traceplot(sim_msom, "lp__", window=c(1,50), inc_warmup=T)



# ==============================
# = Get True and Estimated Psi =
# ==============================
# ---- True Psi ----
use.logit.psi <- FALSE
agg.psi <- FALSE
dim.conv1 <- c(grid.w*grid.h, ns, grid.t, n.obs.reps)

get.psiTrue <- function(x, use.logit=FALSE, agg=FALSE){
	psi.true <- attr(x, "psi")
	psi.true <- array(psi.true, dim=dim.conv1)
	if(use.logit){
		psi.true <- logit(pmax(pmin(psi.true,1-1E-3),1E-3))
	}
	if(agg){
		psi.true <- apply(psi.true, c(1,2,3), mean)
	}
	return(psi.true)
}
man.psi.true <- attributes(big.out.obs[[1]])$psi
psi.true <- get.psiTrue(big.out.obs[[1]], use.logit.psi, agg.psi)
# psi.true <- aperm(apply(psi.true, c(1,2,3), mean), c(3,1,2))
psi.true <- aperm(psi.true[,,,1], c(3,1,2))

# ---- Estimated Psi ----
psi_mean <- plogis(apply(sims$logit_psi, 2:4, mean))

plot(psi.true, psi_mean[,,1:ns])
# plot(psi.true[,,5], psi_mean[,,5])
# txtplot(psi.true, psi_mean[,,1:ns])

Omega <- rstan::extract(sim_msom, "Omega")[[1]]
plot(density(Omega, from=0, to=1))
omega_prior_q <- seq(0,1,length.out=length(Omega))
omega_prior <- dbeta(omega_prior_q, 2, 2)
lines(omega_prior_q, omega_prior, col="blue")


# txtdensity(psi_mean)
# txtdensity(theta_mean)

# ================================
# = Get True and Estimated Theta =
# ================================

# Options
use.logit.p <- FALSE
agg.p <- FALSE

# True p
get.pTrue <- function(x, use.logit=FALSE, agg=FALSE){
	mini.get.pTrue <- function(x){
		op <- attr(x, "obs.params")
		miniP <- t(op$tax.chance)*op$base.chance
	}
	p.true <- lapply(x, mini.get.pTrue)
	p.true <- array(unlist(p.true,F,F),dim=c(ns,grid.t,n.obs.reps))
	
	return(p.true)
	
}
p.true <- aperm(get.pTrue(big.out.obs, use.logit.p, agg.p)[,,1], dim=c(2,1))
theta_mean <- plogis(apply(sims$logit_theta, 2:4, mean))
plot(p.true, apply(theta_mean[,,1:ns], c(1,3), mean))

# =============================
# = Estimated response curves =
# =============================
par(mfrow=c(2,1), mar=c(2.5,2.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=9)
# ---- Species Response Curves via Psi and bt ----
alpha <- apply(sims$alpha, 2:3, mean)
U <- apply(staticData$U, 3, function(x)seq(min(x), max(x), length.out=100))
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

plot(V[,2], theta_resp[,1], ylim=range(theta_resp), type='l', col="gray", xlab="year")
for(i in 2:ncol(theta_resp)){
	lines(V[,2], theta_resp[,i], col="gray")
}
par(new=T)
plot(density(staticData$V[,,"yr"], from=min(staticData$V[,,"yr"]), to=max(staticData$V[,,"yr"])), xaxt="n",yaxt="n", ylab="",xlab="", main="",type="l", col="blue")

# =====================================
# = True Values in Estimate Intervals =
# =====================================
Rhat_all <- summary(sim_msom)[[1]][,"Rhat"]
Rhat_psi <- aperm(psi_mean, 3:1) # just for array structure
Rhat_psi[] <- Rhat_all[grepl("logit_psi", names(Rhat_all))]
Rhat_psi <- aperm(Rhat_psi, 3:1)
Rhat_psi_good <- (Rhat_psi > 0.99 & Rhat_psi < 1.01)[,,1:ns]

quant_prob <- 0.025
psi_lower <- plogis(apply(sims$logit_psi, 2:4, quantile, probs=quant_prob))
psi_upper <- plogis(apply(sims$logit_psi, 2:4, quantile, probs=1-quant_prob))

sum(psi.true > psi_lower[,,1:ns] & psi.true < psi_upper[,,1:ns])/(length(psi.true))


o_psi <- order(psi.true)
plot(psi_lower[o_psi], ylim=range(c(psi_upper, psi_lower)), type="l", col="blue")
lines(psi_upper[o_psi], type="l", col="red")
lines(psi.true[o_psi], col="black")

# sum(psi.true[Rhat_psi_good] > psi_lower[,,1:ns][Rhat_psi_good] & psi.true[Rhat_psi_good] < psi_upper[,,1:ns][Rhat_psi_good])/(sum(Rhat_psi_good))

# =================
# = Compare Alpha =
# =================
alpha_123 <- apply(sims$alpha, 2:3, mean)[,1:ns]
alpha1_true <- attr(big.out.obs[[1]], "u.a0")
alpha2_true <- attr(big.out.obs[[1]], "a3")
# alpha3_true <- attr(big.out.obs[[1]], "a4")
# alpha_true <- list(alpha1_true, alpha2_true, alpha3_true)
alpha_true <- list(alpha1_true, alpha2_true)

par(mfrow=c(2,1), mar=c(2,2,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=10)
plot(alpha1_true, alpha_123[1,])
plot(alpha2_true, alpha_123[2,])
# plot(alpha3_true, alpha_123[3,])


# ---- Boxplots of posterior, with true in blue, ALPHA ----
par(mfrow=c(2,1), mar=c(2,2,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=10)
alpha_123_df <- reshape2::melt(sims$alpha[,,1:ns])
for(i in 1:2){ # 3 alphas
	boxplot(value~Var3+Var2, data=alpha_123_df[alpha_123_df[,"Var2"]==i,])
	points(alpha_true[[i]], col="blue", pch=19)
}


