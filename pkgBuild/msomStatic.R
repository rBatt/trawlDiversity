library(devtools)
library("trawlData")
load_all("trawl/trawlDiversity")

# ======================
# = Subset Data to Use =
# ======================
ebs.a2 <- ebs.a[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]

# ebs.a1 <- ebs.agg2[pick(spp, 50, w=T)][pick(year,3)][pick(stratum, 30)]
# ebs.a2 <- ebs.a1[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]

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
staticData <- msomData(ebs.a2, n0=10, cov.vars=c(bt="btemp",doy="doy",yr="year"), u.form=~bt+I(bt^2), v.form=~doy+I(doy^2)+yr, valueName="abund")


# =====================
# = Fit Model in Stan =
# =====================
library(rstan)
model_file <- "trawl/trawlDiversity/inst/stan/msomStatic.stan"
ebs_msom <- stan(
	file=model_file, 
	data=c("X","U","V","nK","nT","Kmax","Jmax","nU","nV","nS","N"), 
	# control=list(stepsize=0.05, adapt_delta=0.95), 
	chains=4, iter=500, refresh=1, seed=1337, cores=4
)


# ==================================
# = Printing the Fitted Parameters =
# ==================================
print(ebs_msom, c("alpha[1,1]", "beta[1,1]", "Omega"));

print(ebs_msom, c(
	"alpha[1,1]", "alpha[2,1]", "alpha[3,1]", 
	"beta[1,1]", "beta[2,1]", "beta[3,1]", "beta[4,1]", "beta[5,1]",
	"Omega"
))