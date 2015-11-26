library(devtools)
library("trawlData")
load_all("trawl/trawlDiversity")

# ======================
# = Subset Data to Use =
# ======================
ebs.a2 <- ebs.a[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]

# ebs.a1 <- ebs.agg2[pick(spp, 50, w=T)][pick(year,3)][pick(stratum, 30)]
# ebs.a2 <- ebs.a1[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]

# aggregate process covariates among samples within site-year
ebs.a2[,btemp:=meanna(btemp), by=c("stratum","year")]

# hacked bad word around for dropping
ebs.a2[,cantFill:=all(is.na(btemp)),by=c("year","stratum")]
ebs.a2 <- ebs.a2[!(cantFill)]

# scale day of year
doy.mu <- ebs.a2[,mean(doy, na.rm=TRUE)]
doy.sd <- ebs.a2[,sd(doy, na.rm=TRUE)]
ebs.a2[,doy:=(doy-doy.mu)/doy.sd]


# ======================
# = Cast Data for Stan =
# ======================
staticData <- msomData(Data=ebs.a2, n0=10, cov.vars=c(bt="btemp",doy="doy",yr="year"), u.form=~bt+I(bt^2), v.form=~doy+I(doy^2)+yr, valueName="abund")

# drop K-dimension duplicates of process covariates
staticData$U <- staticData$U[,,1,]

# add a counter for nJ (number of sites in each year)
staticData$nJ <- apply(staticData$nK, 1, function(x)sum(x>0))


# =====================
# = Fit Model in Stan =
# =====================
library(rstan)
model_file <- "trawl/trawlDiversity/inst/stan/msomStatic.stan"
# ebs_msom <- stan(
# 	file=model_file,
# 	data=c("X","U","V","nK","nT","Kmax","Jmax","nU","nV","nS","N"),
# 	# control=list(stepsize=0.05, adapt_delta=0.95),
# 	chains=4, iter=500, refresh=1, seed=1337, cores=4
# )
ebs_msom <- stan(
	file=model_file, 
	data=staticData, 
	control=list(stepsize=0.05, adapt_delta=0.95),
	chains=4, iter=50, refresh=1, seed=1337, cores=4, verbose=FALSE
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
print(ebs_msom, c("alpha[1,1]", "beta[1,1]", "Omega"));

inspect_params <- c(
	"alpha_mu","alpha_sd","beta_mu","beta_sd",
	"alpha[1,1]", "alpha[2,1]", "alpha[3,1]", 
	"beta[1,1]", "beta[2,1]", "beta[3,1]", "beta[4,1]", "beta[5,1]",
	# "logit_psi[1,1,1,1]","logit_psi[1,1,1,5]","logit_psi[1,1,1,9]",
	# "logit_theta[1,1,1,1]","logit_theta[1,1,1,5]","logit_theta[1,1,1,9]",
	"Omega"
)



print(ebs_msom, inspect_params)

traceplot(ebs_msom, inspect_params, inc_warmup=FALSE)

pairs(ebs_msom, pars=c("alpha[1,1]","beta[1,1]"))

hist_treedepth <- function(fit) { 
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE) 
  hist(sapply(sampler_params, function(x) c(x[,'treedepth__']))[,1], breaks=0:20, main="", xlab="Treedepth") 
  abline(v=10, col=2, lty=1) 
}
sapply(ebs_msom@sim$samples, function(x) attr(x, 'args')$control$max_treedepth)
hist_treedepth(ebs_msom)