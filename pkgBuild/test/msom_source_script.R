
# ========
# = Load =
# ========
library("devtools")
library("trawlData")
library("trawlDiversity")
library("rbLib")
library("rstan")


# run on amphiprion: 
# nohup R CMD BATCH -cwd --no-save trawlDiversity/pkgBuild/test/msom_source_script.R &

# ======================
# = Subset Data to Use =
# ======================
# smallest data set
# ebs.a2 <- ebs.a[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]
#
# # medium data set
# set.seed(1337)
# ind <- mpick(ebs.agg2, p=c(stratum=10, year=10), weight=TRUE, limit=60)
# logic <- expression(
# 	spp%in%spp[ind]
# 	& stratum%in%stratum[ind]
# 	& year%in%year[ind]
# )
# ebs.a1 <- ebs.agg2[eval(logic)][pick(spp, 20, w=FALSE)]
# ebs.a2 <- ebs.a1[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, stemp=stemp, depth=depth, doy=yday(datetime))]

# # Medium-large data set
# set.seed(1337)
# ebs.a1 <- ebs.agg2[pick(spp, 100, w=TRUE)]
# ebs.a2 <- ebs.a1[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, stemp=stemp, depth=depth, doy=yday(datetime))]

# largest data set
ebs.a2 <- ebs.agg2[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, stemp=stemp, depth=depth, doy=yday(datetime))]


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
model_file <- "trawlDiversity/inst/stan/msomStatic.stan"
# model_file <- "trawlDiversity/inst/stan/msomDynamic.stan"

tag <- paste0("start_", format.Date(Sys.time(),"%Y-%m-%d_%H-%M-%S"))

save.image(paste0("trawlDiversity/pkgBuild/test/msomStatic_fullEBS_preSave_",tag,".RData"))

sessionInfo()

ebs_msom <- stan(
	file=model_file, 
	data=stanData, 
	control=list(stepsize=0.01, adapt_delta=0.95, max_treedepth=15),
	chains=1, iter=200, seed=1337, cores=1, verbose=F, refresh=1
)

save.image(renameNow(paste0("trawlDiversity/pkgBuild/test/msomStatic_fullEBS_",tag,"_end",".RData")), compress="xz")
