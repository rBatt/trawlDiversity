library(devtools)
library("trawlData")
load_all("trawl/trawlDiversity")

# ======================
# = Subset Data to Use =
# ======================
ebs.a2 <- ebs.a[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]
# ebs.a2 <- ebs.agg2[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]

# ebs.a1 <- ebs.agg2[pick(year,5)][pick(stratum, 7, w=T)][pick(spp, 10, w=FALSE)]

# set.seed(1337)
# ind <- mpick(ebs.agg2, p=c(stratum=5, year=4), weight=TRUE, limit=60)
# logic <- expression(
# 	spp%in%spp[ind]
# 	& stratum%in%stratum[ind]
# 	& year%in%year[ind]
# )
# ebs.a1 <- ebs.agg2[eval(logic)][pick(spp, 10, w=FALSE)]

# ebs.a2 <- ebs.a1[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]


# aggregate process covariates among samples within site-year
sd2 <- function(x, ...){
	if(sum(!is.na(x))==1){return(0)}else{return(sd(x, ...))}
}
ebs.a2[,bt_sd:=sd2(.SD[,list(btemp=mean(btemp,na.rm=T)),by="K"][,btemp], na.rm=TRUE), by=c("stratum","year")]
ebs.a2[,btemp:=meanna(btemp), by=c("stratum","year")]
ebs.a2[,bt2_sd:=2*btemp*bt_sd]


# hacked bad word around for dropping
ebs.a2[,cantFill:=all(is.na(btemp)),by=c("year","stratum")]
ebs.a2 <- ebs.a2[!(cantFill)]

# scale day of year
doy.mu <- ebs.a2[,mean(doy, na.rm=TRUE)]
doy.sd <- ebs.a2[,sd(doy, na.rm=TRUE)]
ebs.a2[,doy:=(doy-doy.mu)/doy.sd]

# aggregate doy
ebs.a2[,doy_sd:=sd2(.SD[,list(doy=mean(doy,na.rm=T)),by="K"][,doy], na.rm=TRUE), by=c("stratum","year")]
ebs.a2[,doy:=meanna(doy), by=c("stratum","year")]
ebs.a2[,doy2_sd:=2*doy*doy_sd]



# ======================
# = Cast Data for Stan =
# ======================
# cov.vars <- c(bt="btemp",doy="doy",yr="year")
staticData <- msomData(Data=ebs.a2, n0=10, cov.vars=c(bt="btemp",doy="doy",yr="year"), u.form=~bt+I(bt^2), v.form=~doy+I(doy^2)+yr, valueName="abund")

staticData_sd <- msomData(Data=ebs.a2, n0=1, cov.vars=c(bt_sd="bt_sd",doy_sd="doy_sd",bt2_sd="bt2_sd",doy2_sd="doy2_sd",yr="year"), u.form=~bt_sd+bt2_sd-1, v.form=~doy_sd+doy2_sd+yr, valueName="abund")[c("U","V")]

# drop K-dimension of process covariates
# add in sd of covariates (sd among K dimension)
staticData$U <- staticData$U[,,1,-1]
staticData$U_sd <- staticData_sd$U[,,1,]
names(staticData)[names(staticData)=="U"] <- "U_mu"

# drop K-dimension of detection covariates
# add sd of detection covariates (among samples)
staticData$V <- staticData$V[,,1,-1]
staticData$V_sd <- staticData_sd$V[,,1,-1]
names(staticData)[names(staticData)=="V"] <- "V_mu"

# add a counter for nJ (number of sites in each year)
staticData$nJ <- apply(staticData$nK, 1, function(x)sum(x>0))

# aggregate abundance across samples
staticData$X <- apply(staticData$X, c(1,2,4), function(x)sum(x))

# add binary presence/ absence
# isObs <- staticData$X
# isObs[] <- pmin(1, staticData$X)
# staticData$isObs <- isObs

# =====================
# = Fit Model in Stan =
# =====================
library(rstan)
model_file <- "trawl/trawlDiversity/inst/stan/msomStatic.stan"

ebs_msom <- stan(
	file=model_file, 
	data=staticData, 
	control=list(stepsize=0.01, adapt_delta=0.95),
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
	"Omega[1]", "Omega[3]"
)



print(ebs_msom, inspect_params)

traceplot(ebs_msom, inspect_params, inc_warmup=T)

pairs(ebs_msom, pars=c("alpha[1,1]","beta[1,1]"))

hist_treedepth <- function(fit) { 
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE) 
  hist(sapply(sampler_params, function(x) c(x[,'treedepth__']))[,1], breaks=0:20, main="", xlab="Treedepth") 
  abline(v=10, col=2, lty=1) 
}
sapply(ebs_msom@sim$samples, function(x) attr(x, 'args')$control$max_treedepth)
hist_treedepth(ebs_msom)

Omega <- apply(rstan::extract(ebs_msom, "Omega")[[1]], 2, mean)
plot(density(Omega, from=0, to=1))
omega_prior_q <- seq(0,1,length.out=length(Omega))
omega_prior <- dbeta(omega_prior_q, 2, 2)
lines(omega_prior_q, omega_prior, col="blue")

