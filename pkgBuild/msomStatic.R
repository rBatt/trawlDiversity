load_all("trawlData")
load_all("trawl/trawlDiversity")

ebs.a2 <- ebs.a[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]

# Cast Data
# Not quite ready for Stan
# Order of dimensions is optimized for Stan looping:
# last dimension loops most quickly
Yc <- trawlCast(x=ebs.a2, formula=year~stratum~K~spp, valueName="abund", grandNamesOut=c("t","j","k","s"))

nK <- trawlCast(ebs.a2, 
	year~stratum, 
	valueName="K", 
	fixAbsent=FALSE, 
	fun.aggregate=max, 
	valFill=0, 
	grandNamesOut=c("j","t")
) # used to indicate which values in Yc are NA, basically


# Get covariates for U
U.tjk <- ebs.a2[,list(btemp=una(btemp,na.rm=TRUE), doy=una(doy,na.rm=TRUE)),by=c("year","stratum","K")]
U.tjk <- U.tjk[,list(K=K, btemp=fill.mean(btemp), doy=fill.mean(doy)),by=c("year","stratum")]
U <- 


Yc <- trawlCast(x=ebs.a2, formula=year~stratum~K, valueName="doy", grandNamesOut=c("t","j","k","year"), fun.aggregate=meanna)

U_year <- trawlCast(ebs.a2, 
	year~stratum, 
	valueName="K", 
	fixAbsent=FALSE, 
	fun.aggregate=max, 
	valFill=0, 
	grandNamesOut=c("j","t")
) # used to indicate which values in Yc are NA, basically


# test if it'd work in stan
Jmax <- length(dimnames(Yc)$j)
Kmax <- length(dimnames(Yc)$k)
nS <- length(dimnames(Yc)$s)
nT <- length(dimnames(Yc)$t)
for(t in 1:Tmax){
	for(j in 1:Jmax){
		t.K <- nK[t,j]
		if(t.K){
			for(k in 1:t.K){
				t.dat <- unname(Yc[t,j,k,])
				print(t.dat)
			}
		}
	}
}
# if you don't see any NA's, should all be good!


# Convert NA's to 0's, so they Stan doesn't get mad.
# But don't worry, these converts will be skipped by loop,
# thanks to indexing provided by nK. 
# Reminder: nK tells us how many reps, 
# or if 0, that a J-T combo doesn't exist
Y <- Yc
Y[is.na(Y)] <- 0






