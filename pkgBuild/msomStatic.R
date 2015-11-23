library(devtools)
load_all("trawlData")
load_all("trawl/trawlDiversity")

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

# Cast Data
# Not quite ready for Stan
# Order of dimensions is optimized for Stan looping:
# last dimension loops most quickly
n0 <- 10
gno <- c("year","stratum","K","spp") # need to get this by formula parsing
formula <- year~stratum~K~spp
cov.vars <- c(bt="btemp",doy="doy",yr="year") # order does not matter


Xc <- trawlCast(x=ebs.a2, formula=formula, valueName="abund", grandNamesOut=gno)
nK <- trawlCast(ebs.a2, 
	year~stratum, 
	valueName=tail(gno,2)[1], #"K", 
	fixAbsent=FALSE, 
	fun.aggregate=max, 
	valFill=0, 
	grandNamesOut=head(gno, length(gno)-2)# c("j","t")
) # used to indicate which values in Xc are NA, basically


# Get covariates for U
# Choose names of covariates
# The names of the vector are the names you'll get in the end
# The character elements of the vector should correspond to columns
# cov.vars <- c(bt="btemp",doy="doy",yr="year") # order does not matter
# cov.vars <- c(btemp="bt",doy="doy",year="yr") # if wanted opposite convention

# Grouping for covariates
# Correspond to column names
cov.by <- c("year","stratum","K") # the order matters! most specific last

# Aggregate covariates
# Define expression to agg by 
# applying una() function
# Sets names, too
una.cov <- expression({
	structure(lapply(eval(s2c(cov.vars)), una, na.rm=TRUE),.Names=names(cov.vars))
})
cov.tjk <- ebs.a2[,eval(una.cov), keyby=cov.by] # also setting key

# Fill in NA covariates
# Define expression & functions to fill using mean
# It is implied that to get the mean, you look 1 level higher
# than the most specific favor in cov.by
# This is why the order of cov.by matters, a lot
# PRETTX PRESUMPTUOUS / FRAGILE CODE
is.ci <- function(x)is.numeric(x) | is.integer(x)
fm2 <- function(x){
	if(!any(is.na(x))){
		return(x)
	}
	if(!is.ci(x)){
		warning("Covariate contains NA's, but fill.mean needs numeric or integer")
		return(x)
	}else{
		cl <- class(x)
		as(fill.mean(x), cl)
	}	
}
fillMean.cov <- expression({
	structure(
		c(
			# eval(s2c(cov.by[length(cov.by)])),
			lapply(eval(s2c(names(cov.vars))), fm2)
		),
		.Names=names(cov.vars)
	)
})

cov.tjk[,c(names(cov.vars)):=eval(fillMean.cov),by=c(cov.by[-length(cov.by)])]



# Check
stopifnot(!any(is.na(cov.tjk))) # can't have any NA's with my current simple approach

# Set up template for expanding covariates
template <- unique(data.table(reshape2:::melt(Xc), key=c(cov.by)))[,eval(s2c(cov.by))]
template[,c(names(template)):=lapply(.SD, as.character)]
cov.tjk[,c(cov.by):=lapply(.SD[,eval(s2c(cov.by))], as.character)]

# Fill out (expand) covariate data.table
cov.f <- merge(template, cov.tjk, all=TRUE, by=cov.by) # filled cov

# cov.tjk.long <- data.table:::melt.data.table(cov.tjk, id.vars=c("stratum","K","year"), measure.vars=c("btemp","doy","t"))
# cov.a <- trawlCast(cov.f,
# 	year~stratum,
# 	valueName="K",
# 	fixAbsent=FALSE,
# 	fun.aggregate=max,
# 	valFill=0,
# 	grandNamesOut=c("t","j")
# )

getUV <- function(form){
	UV.m <- model.matrix(form, model.frame(form,data=cov.f, na.action=na.pass))
	ncUV <- ncol(UV.m)
	nrUV <- nrow(UV.m)
	namesUV <- gsub("I\\(|\\(|\\)|\\^", "", colnames(UV.m))
	colnames(UV.m) <- namesUV

	UV.dt0 <- data.table(cov.f[,eval(s2c(cov.by))],UV.m)
	UV.dt <- data.table:::melt.data.table(UV.dt0, id.var=cov.by)
	# U.form <- paste(c(tail(cov.by,1),"variable",rev(head(cov.by,-1))),collapse="~")
	UV.form <- paste(c(cov.by,"variable"),collapse="~")
	
	reshape2::acast(UV.dt, formula=UV.form, value.var="value")
}

u.form <- ~bt+I(bt^2)
U <- getUV(u.form)
nU <- tail(dim(U),1)

v.form <- ~doy+I(doy^2)+yr
V <- getUV(v.form)
nV <- tail(dim(V),1)


# Xc <- trawlCast(x=ebs.a2, formula=year~stratum~K, valueName="doy", grandNamesOut=c("t","j","k","year"), fun.aggregate=meanna)




# if you don't see any NA's, should all be good! (in old version, with t.K loop)


# Add never-observed species to array
add_neverObs <- function(x, n0){
	fillA <- do.call(`[`, c(list(x), rep(TRUE, length(dim(x))-1), 1))
	fillA[!is.na(fillA)] <- 0
	X0 <- replicate(n0, fillA)
	
	outDim <- c(head(dim(x),-1), tail(dim(x),1)+n0)
	outA <- array(c(x, X0), dim=outDim)
	
	return(outA)
	
}
X <- add_neverObs(Xc, n0=n0)


# test if it'd work in stan
nT <- length(dimnames(Xc)$year)
Jmax <- length(dimnames(Xc)$stratum)
Kmax <- length(dimnames(Xc)$K)
N <- length(dimnames(Xc)$spp)
nS <- N + n0

# for(t in 1:nT){
# 	for(j in 1:Jmax){
# 		t.K <- nK[t,j]
# 		if(t.K){
# 			t.dat <- unname(Xc[t,j,,])
# 			print(t.dat)
# 			# for(k in 1:t.K){
# #
# # 			}
# 		}
# 	}
# }


# Convert NA's to 0's, so they Stan doesn't get mad.
# But don't worry, these converts will be skipped by loop,
# thanks to indexing provided by nK. 
# Reminder: nK tells us how many reps, 
# or if 0, that a J-T combo doesn't exist
X[is.na(X)] <- 0

U[is.na(U)] <- 0
V[is.na(V)] <- 0


# =====================
# = Fit Model in Stan =
# =====================
library(rstan)
model_file <- "trawl/trawlDiversity/inst/stan/msomStatic.stan"
ebs_msom <- stan(
	file=model_file, 
	data=c("X","U","V","nK","nT","Kmax","Jmax","nU","nV","nS","N"), 
	# control=list(stepsize=0.05, adapt_delta=0.95), 
	chains=1, iter=500, refresh=1, seed=1337, cores=1
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