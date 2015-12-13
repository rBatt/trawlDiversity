
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
ebs.a2 <- ebs.agg2[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, stemp=stemp, doy=yday(datetime))]


# ======================================
# = Aggregate and Transform Covariates =
# ======================================
# rename columns for shorthand
setnames(ebs.a2, c("btemp"), c("bt"))
setnames(ebs.a2, c("stemp"), c("st"))
ebs.a2[,yr:=scale(as.integer(year))]

# aggregate and transform (^2) btemp
mk_cov_rv_pow(ebs.a2, "bt", across="K", by=c("stratum","year"), pow=2)
mk_cov_rv_pow(ebs.a2, "st", across="K", by=c("stratum","year"), pow=2)

# scale and aggregate doy
doy.mu <- ebs.a2[,mean(doy, na.rm=TRUE)]
doy.sd <- ebs.a2[,sd(doy, na.rm=TRUE)]
ebs.a2[,doy:=(doy-doy.mu)/doy.sd]
mk_cov_rv(ebs.a2, "doy", across="K", by=c("stratum","year"))


# ======================
# = Cast Data for Stan =
# ======================
# ---- Get Basic Structure of MSOM Data Input ----
dynData <- msomData(Data=ebs.a2, n0=2, cov.vars=c(bt="bt",bt2="bt2",yr="yr", doy="doy"), u.form=~bt+bt2, v.form=~doy+year, valueName="abund", cov.by=c("year","stratum"), v_rv=c("doy"))

# ---- Add a counter for nJ (number of sites in each year) ----
dynData$nJ <- apply(dynData$nK, 1, function(x)sum(x>0))

# ---- Aggregate abundance across samples ----
dynData$X <- apply(dynData$X, c(1,2,4), function(x)sum(x))


# =====================
# = Fit Model in Stan =
# =====================
model_file <- "trawl/trawlDiversity/inst/stan/msomDynamic.stan"

ebs_msom <- stan(
	file=model_file, 
	data=dynData, 
	control=list(stepsize=0.01, adapt_delta=0.95, max_treedepth=15),
	chains=4, iter=200, seed=1337, cores=4, verbose=F
)


save.image(renameNow("trawl/trawlDiversity/pkgBuild/test/msomDynamic_fullEBS.RData"))