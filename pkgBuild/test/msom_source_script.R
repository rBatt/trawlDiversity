
# ========
# = Load =
# ========
library("devtools")
library("trawlData")
library("trawlDiversity")
library("rbLib")
library("rstan")
library("parallel")


# =========
# = Cores =
# =========
# if(Sys.info()["sysname"]=="Windows"){
# 	nC <- floor(detectCores()*0.75)
# 	registerDoParallel(cores=nC)
# }else if(Sys.info()["sysname"]=="Linux"){
# 	registerDoParallel(floor(detectCores()*0.5))
# }else{
# 	registerDoParallel()
# }


# ===========
# = Options =
# ===========
chains = 4
cores = 4
iter = 50


# run on amphiprion: 
# nohup R CMD BATCH -cwd --no-save trawlDiversity/pkgBuild/test/msom_source_script.R &

# ======================
# = Subset Data to Use =
# ======================
# smallest data set
ebs.a1 <- trim_msom("ebs", gridSize=1, plot=FALSE)
ebs.a2 <- ebs.a1[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, stemp=stemp, depth=depth, doy=yday(datetime))]


# ==================
# = Set Covariates =
# ==================
cov.vars_use <- c(bt="bt",bt2="bt2",yr="yr", doy="doy", doy2="doy2", depth="depth", depth2="depth2")


# ======================================
# = Aggregate and Transform Covariates =
# ======================================
# rename columns for shorthand
setnames(ebs.a2, c("btemp"), c("bt"))
setnames(ebs.a2, c("stemp"), c("st"))
ebs.a2[,yr:=scale(as.integer(year))]

# aggregate and transform (^2) btemp stemp and yr
mk_cov_rv_pow(ebs.a2, "bt", across="K", by=c("stratum","year"), pow=2)
mk_cov_rv_pow(ebs.a2, "st", across="K", by=c("stratum","year"), pow=2)
mk_cov_rv_pow(ebs.a2, "yr", across="K", by=c("stratum","year"), pow=2)

# scale and aggregate doy
doy.mu <- ebs.a2[,mean(doy, na.rm=TRUE)]
doy.sd <- ebs.a2[,sd(doy, na.rm=TRUE)]
ebs.a2[,doy:=(doy-doy.mu)/doy.sd]
# mk_cov_rv(ebs.a2, "doy", across="K", by=c("stratum","year"))
mk_cov_rv_pow(ebs.a2, "doy", across="K", by=c("stratum","year"), pow=2)

# scale and aggregate depth
depth.mu <- ebs.a2[,mean(depth, na.rm=TRUE)]
depth.sd <- ebs.a2[,sd(depth, na.rm=TRUE)]
ebs.a2[,depth:=(depth-depth.mu)/depth.sd]
# mk_cov_rv(ebs.a2, "depth", across="K", by=c("stratum","year"))
mk_cov_rv_pow(ebs.a2, "depth", across="K", by=c("stratum","year"), pow=2)


# ======================
# = Cast Data for Stan =
# ======================
# ---- Get Basic Structure of MSOM Data Input ----
setkey(ebs.a2, year, stratum, K, spp)
stanData <- msomData(Data=ebs.a2, n0=50, cov.vars=cov.vars_use, u.form=~bt+bt2+depth+depth^2+yr, v.form=~year, valueName="abund", cov.by=c("year","stratum"))

stanData$nJ <- apply(stanData$nK, 1, function(x)sum(x>0)) # number of sites in each year
stanData$X <- apply(stanData$X, c(1,2,4), function(x)sum(x)) # agg abund across samples



# =============
# = Kill NA's =
# =============
# replace NA mu's with 0
# replace NA sd's with 1E3
na_to_0 <- function(x){x[is.na(x)] <- 0; x}
na_to_big <- function(x, big=1E3){x[is.na(x)] <- big; x}
stanData$U_mu <- na_to_0(stanData$U_mu)
stanData$U_sd <- na_to_big(stanData$U_sd)
stanData$V_mu <- na_to_0(stanData$V_mu)
stanData$V_sd <- na_to_big(stanData$V_sd)
stopifnot(!any(sapply(stanData, function(x)any(is.na(x)))))


# =====================
# = Fit Model in Stan =
# =====================
# model_file <- "trawlDiversity/inst/stan/msomStatic.stan"
model_file <- "trawlDiversity/inst/stan/msomDynamic.stan"

tag <- paste0("start", format.Date(Sys.time(),"%Y-%m-%d_%H-%M-%S"))

file_prefix <- paste0(
	"nT",stanData$nT,"_",
	"Jmax",stanData$Jmax, "_",
	"Kmax",stanData$Kmax, "_",
	"N",stanData$N, "_",
	"nS",stanData$nS
)

file_stan_options <- paste0(
	"chains",chains,"_",
	"iter",iter,"_",
	"cores",cores
)

run_info <- paste(file_prefix, file_stan_options, sep="_")

file_path <- file.path("trawlDiversity","pkgBuild","test")
	
save.image(file.path(file_path, paste0("msomDynamic_ebs_", tag,".RData")))

sessionInfo()

ebs_msom_model <- stan_model(file=model_file, auto_write=TRUE)
stan_control <- list(stepsize=0.01, adapt_delta=0.95, max_treedepth=15)

t_out <- sampling(
		object=ebs_msom_model,
		data=stanData,
		control=stan_control,
		pars= c("Omega", "alpha_mu", "alpha_sd", "beta_mu", "beta_sd", "phi_mu", "phi_sd", "gamma_mu", "gamma_sd", "phi", "gamma", "logit_psi", "logit_theta"),
		chains=2, iter=4, seed=1337, cores=1, verbose=F, refresh=1
)


# stan(
# 	file=model_file,
# 	data=stanData,
# 	control=stan_control,
# 	chains=1, iter=2, seed=1337, cores=1, verbose=F, refresh=1
# )

do_stan <- function(x, object, data, iter){
	sampling(
			object=object, 
			data=stanData, 
			control=stan_control,
			pars= c("Omega", "alpha_mu", "alpha_sd", "beta_mu", "beta_sd", "phi_mu", "phi_sd", "gamma_mu", "gamma_sd", "phi", "gamma", "logit_psi", "logit_theta"),
			chains=1, iter=iter, seed=1337, cores=1, verbose=F, refresh=1, chain_id=x
	)
}

out <- do_stan(1, ebs_msom_model, data=stanData, iter=iter)

out <- mclapply(1:chains, do_stan, data=stanData, object=ebs_msom_model, iter=iter, mc.cores=cores, mc.preschedule=FALSE)
# read1 <- read_stan_csv("ebs_msom_samples1")
file.path(file_path, paste0("msomDynamic_ebs_", tag,"_end",".RData"))
save.image(renameNow(file.path(file_path, paste0("msomDynamic_ebs_", tag,"_end",".RData"))), compress="xz")
