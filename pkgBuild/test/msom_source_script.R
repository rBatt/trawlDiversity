
# ========
# = Load =
# ========
library("devtools")
library("trawlData")
library("trawlDiversity")
library("rbLib")
library("rstan")


# ======================
# = Subset Data to Use =
# ======================
# smallest data set
# ebs.a2 <- ebs.a[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]
#
# # medium data set
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


# =============
# = Kill NA's =
# =============
# replace NA's with mean 0 and sd 1E3
na_to_0 <- function(x){x[is.na(x)] <- 0; x}
na_to_big <- function(x, big=1E3){x[is.na(x)] <- big; x}
cov.vars_use_n <- names(cov.vars_use)
cov.vars_use_n_sd <- paste0(cov.vars_use_n, "_sd")
ebs.a2[,c(cov.vars_use_n):=lapply(eval(s2c(cov.vars_use_n)), na_to_0)]
ebs.a2[,c(cov.vars_use_n_sd):=lapply(eval(s2c(cov.vars_use_n_sd)), na_to_big)]
# ebs.a2[stratum=="-171.5 57.5"  & year==1994]


# ======================
# = Cast Data for Stan =
# ======================
# ---- Get Basic Structure of MSOM Data Input ----
stanData <- msomData(Data=ebs.a2, n0=50, cov.vars=cov.vars_use, u.form=~bt+bt2+depth+depth^2+yr, v.form=~doy+doy2+year, valueName="abund", cov.by=c("year","stratum"), v_rv=c("doy","doy2"))

stanData$nJ <- apply(stanData$nK, 1, function(x)sum(x>0)) # number of sites in each year
stanData$X <- apply(stanData$X, c(1,2,4), function(x)sum(x)) # agg abund across samples


# =====================
# = Fit Model in Stan =
# =====================
model_file <- "trawlDiversity/inst/stan/msomStatic.stan"
# model_file <- "trawlDiversity/inst/stan/msomDynamic.stan"

save.image(renameNow("trawlDiversity/pkgBuild/test/msomStatic_fullEBS_preSave.RData"))

ebs_msom <- stan(
	file=model_file, 
	data=stanData, 
	control=list(stepsize=0.01, adapt_delta=0.95, max_treedepth=15),
	chains=4, iter=150, seed=1337, cores=4, verbose=F
)


save.image(renameNow("trawlDiversity/pkgBuild/test/msomStatic_fullEBS.RData"), compress="xz")