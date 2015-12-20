
# ========
# = Load =
# ========
library("devtools")
library("trawlData")
library("trawlDiversity")
library("rbLib")
library("raster")
library("rstan")
library("R2jags")
library("parallel")


# ===========
# = Options =
# ===========
chains = 4
cores = 4
iter = 50
language = c("JAGS", "Stan")[2]

model_type = c("Dynamic", "Static")[1]

reg = c(
	"ai", "ebs", "gmex", 
	"goa", "neus", "newf", 
	"ngulf", "sa", "sgulf", 
	"shelf", "wcann", "wctri"
)[2]

params_out <- c(
	"Omega", "alpha_mu", "alpha_sd", 
	"beta_mu", "beta_sd", "phi_mu", 
	"phi_sd", "gamma_mu", "gamma_sd", 
	"phi", "gamma", "logit_psi", 
	"logit_theta"
)


# run on amphiprion: 
# nohup R CMD BATCH -cwd --no-save trawlDiversity/pkgBuild/test/msom_source_script.R &


# ======================
# = Subset Data to Use =
# ======================
regX.a1 <- trim_msom(reg, gridSize=1, plot=FALSE)


# smallest data set
set.seed(1337)
ind <- mpick(regX.a1, p=c(stratum=12, year=15), weight=TRUE, limit=60)
logic <- expression(
	spp%in%spp[ind]
	& stratum%in%stratum[ind]
	& year%in%year[ind]
)
regX.a2 <- regX.a1[eval(logic)][pick(spp, 20, w=FALSE)]

regX.a2 <- regX.a2[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, stemp=stemp, depth=depth, doy=yday(datetime))]

# smallest data set
# regX.a2 <- regX.a[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]


# ==================
# = Set Covariates =
# ==================
cov.vars_use <- c(bt="bt",bt2="bt2",yr="yr", doy="doy", doy2="doy2", depth="depth", depth2="depth2")


# ======================================
# = Aggregate and Transform Covariates =
# ======================================
# rename columns for shorthand
setnames(regX.a2, c("btemp"), c("bt"))
setnames(regX.a2, c("stemp"), c("st"))
regX.a2[,yr:=scale(as.integer(year))]

# aggregate and transform (^2) btemp stemp and yr
mk_cov_rv_pow(regX.a2, "bt", across="K", by=c("stratum","year"), pow=2)
mk_cov_rv_pow(regX.a2, "st", across="K", by=c("stratum","year"), pow=2)
mk_cov_rv_pow(regX.a2, "yr", across="K", by=c("stratum","year"), pow=2)

# scale and aggregate doy
doy.mu <- regX.a2[,mean(doy, na.rm=TRUE)]
doy.sd <- regX.a2[,sd(doy, na.rm=TRUE)]
regX.a2[,doy:=(doy-doy.mu)/doy.sd]
mk_cov_rv(regX.a2, "doy", across="K", by=c("stratum","year"))
mk_cov_rv_pow(regX.a2, "doy", across="K", by=c("stratum","year"), pow=2)

# scale and aggregate depth
depth.mu <- regX.a2[,mean(depth, na.rm=TRUE)]
depth.sd <- regX.a2[,sd(depth, na.rm=TRUE)]
regX.a2[,depth:=(depth-depth.mu)/depth.sd]
mk_cov_rv(regX.a2, "depth", across="K", by=c("stratum","year"))
mk_cov_rv_pow(regX.a2, "depth", across="K", by=c("stratum","year"), pow=2)


# ======================
# = Cast Data for Stan =
# ======================
# ---- Get Basic Structure of MSOM Data Input ----
setkey(regX.a2, year, stratum, K, spp)
inputData <- msomData(Data=regX.a2, n0=50, cov.vars=cov.vars_use, u.form=~bt+bt2+depth+depth^2+yr, v.form=~year, valueName="abund", cov.by=c("year","stratum"))

inputData$nJ <- apply(inputData$nK, 1, function(x)sum(x>0)) # number of sites in each year
inputData$X <- apply(inputData$X, c(1,2,4), function(x)sum(x)) # agg abund across samples



# =============
# = Kill NA's =
# =============
# replace NA mu's with 0
# replace NA sd's with 1E3
na_to_0 <- function(x){x[is.na(x)] <- 0; x}
na_to_big <- function(x, big=1E3){x[is.na(x)] <- big; x}
inputData$U_mu <- na_to_0(inputData$U_mu)
inputData$U_sd <- na_to_big(inputData$U_sd)
inputData$V_mu <- na_to_0(inputData$V_mu)
inputData$V_sd <- na_to_big(inputData$V_sd)
stopifnot(!any(sapply(inputData, function(x)any(is.na(x)))))


model_file <- paste0("msom", model_type, ".", tolower(language))
model_path <- file.path("trawlDiversity", "inst", tolower(language), model_file)

tag <- paste0("start", format.Date(Sys.time(),"%Y-%m-%d_%H-%M-%S"))

file_prefix <- paste0(
	"nT",inputData$nT,"_",
	"Jmax",inputData$Jmax, "_",
	"Kmax",inputData$Kmax, "_",
	"N",inputData$N, "_",
	"nS",inputData$nS
)

file_options <- paste0(
	"chains",chains,"_",
	"iter",iter,"_",
	"cores",cores, 
	"_",language
)

run_info <- paste(file_prefix, file_options, sep="_")

save_file <- paste0("msom", model_type, "_", reg, "_", tolower(language), "_", tag, ".RData")
save_path <- file.path("trawlDiversity","pkgBuild","test",save_file)

# save.image(save_path)
sessionInfo()


# =====================
# = Fit Model in Stan =
# =====================

stan_control <- list(stepsize=0.01, adapt_delta=0.95, max_treedepth=15)

out <- stan(
	file=model_path,
	data=inputData,
	control=stan_control,
	chains=chains, iter=iter, seed=1337, cores=cores, verbose=F, refresh=1
)


# ---- alternative approach to doing stan model, maybe ----
# probably requires installing StanHeaders and rstan from source on development GitHub
# msom_model <- stan_model(file=model_path, auto_write=TRUE)
#
# do_stan <- function(x, object, data, iter){
# 	sampling(
# 			object=object,
# 			data=inputData,
# 			control=stan_control,
# 			pars= params_out,
# 			chains=1, iter=iter, seed=1337, cores=1, verbose=F, refresh=1, chain_id=x
# 	)
# }
#
#
# t_out <- sampling(
# 		object=msom_model,
# 		data=inputData,
# 		control=stan_control,
# 		pars= params_out,
# 		chains=2, iter=4, seed=1337, cores=1, verbose=F, refresh=1
# )
#
# out <- do_stan(1, msom_model, data=inputData, iter=iter)
# # read1 <- read_stan_csv("ebs_msom_samples1")
#
# out <- mclapply(1:chains, do_stan, data=inputData, object=msom_model, iter=iter, mc.cores=cores, mc.preschedule=FALSE)



save.image(renameNow(save_path), compress="xz")
