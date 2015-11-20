load_all("trawlData")
load_all("trawl/trawlDiversity")

ebs.a2 <- ebs.a[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]

# Cast Data
# Not quite ready for Stan
# Order of dimensions is optimized for Stan looping:
# last dimension loops most quickly
gno <- c("year","stratum","K","spp") # need to get this by formula parsing
Yc <- trawlCast(x=ebs.a2, formula=year~stratum~K~spp, valueName="abund", grandNamesOut=gno)
nK <- trawlCast(ebs.a2, 
	year~stratum, 
	valueName=tail(gno,2)[1], #"K", 
	fixAbsent=FALSE, 
	fun.aggregate=max, 
	valFill=0, 
	grandNamesOut=head(gno, length(gno)-2)# c("j","t")
) # used to indicate which values in Yc are NA, basically


# Get covariates for U
# Choose names of covariates
# The names of the vector are the names you'll get in the end
# The character elements of the vector should correspond to columns
cov.vars <- c(bt="btemp",doy="doy",yr="year") # order does not matter
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
# PRETTY PRESUMPTUOUS / FRAGILE CODE
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
template <- unique(data.table(reshape2:::melt(Yc), key=c(cov.by)))[,eval(s2c(cov.by))]
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


# Yc <- trawlCast(x=ebs.a2, formula=year~stratum~K, valueName="doy", grandNamesOut=c("t","j","k","year"), fun.aggregate=meanna)



# test if it'd work in stan
nT <- length(dimnames(Yc)$year)
Jmax <- length(dimnames(Yc)$stratum)
Kmax <- length(dimnames(Yc)$K)
nS <- length(dimnames(Yc)$spp)

for(t in 1:Tmax){
	for(j in 1:Jmax){
		t.K <- nK[t,j]
		if(t.K){
			t.dat <- unname(Yc[t,j,,])
			print(t.dat)
			# for(k in 1:t.K){
#
# 			}
		}
	}
}
# if you don't see any NA's, should all be good! (in old version, with t.K loop)


# Convert NA's to 0's, so they Stan doesn't get mad.
# But don't worry, these converts will be skipped by loop,
# thanks to indexing provided by nK. 
# Reminder: nK tells us how many reps, 
# or if 0, that a J-T combo doesn't exist
Y <- Yc
Y[is.na(Y)] <- 0

U[is.na(U)] <- 0
V[is.na(V)] <- 0


# items to use: Y, U, V, nK ...
# nT, Jmax, Kmax, nS ...
# nV, nU ...
# those should be passed to Stan as data