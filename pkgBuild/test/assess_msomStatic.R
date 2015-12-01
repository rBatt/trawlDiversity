
# =================
# = Load Packages =
# =================
# Data structure
library(raster)
library(data.table)

# Graphing
library(fields)

# Statistics
library(igraph)
library(R2jags)

# Computing
library(parallel)
library(doParallel)
library(foreach)

# Other
library(rbLib) # library(devtools); install_github("rBatt/rbLib")
library(devtools)
library("trawlData")
library("trawlDiversity")

# ===============================
# = Guess appropriate directory =
# ===============================
if(Sys.info()["sysname"]=="Linux"){
	setwd("~/Documents/School&Work/pinskyPost")
}else{
	setwd("~/Documents/School&Work/pinskyPost")
}


# ======================
# = Load Sim Functions =
# ======================
sim.location <- "~/Documents/School&Work/pinskyPost/trawl/Scripts/SimFunctions"
invisible(sapply(paste(sim.location, list.files(sim.location), sep="/"), source, .GlobalEnv))

data.location <- "./trawl/Scripts/DataFunctions"
invisible(sapply(paste(data.location, list.files(data.location), sep="/"), source, .GlobalEnv))

stat.location <- "./trawl/Scripts/StatFunctions"
invisible(sapply(paste(stat.location, list.files(stat.location), sep="/"), source, .GlobalEnv))


# =====================================
# = Set cores for parallel processing =
# =====================================
if(Sys.info()["sysname"]=="Windows"){
	nC <- floor(detectCores()*0.75)
	registerDoParallel(cores=nC)
}else if(Sys.info()["sysname"]=="Linux"){
	# registerDoParallel(cores=min(c(25,floor(detectCores()*0.75))))
	registerDoParallel(floor(detectCores()*0.8))
	# registerDoParallel(floor(detectCores()*0.90))
}else{
	registerDoParallel()
}


# ================
# = Grid Options =
# ================
# Grid Size
grid.w <- 9 # Width # 6
grid.h <- 9 # Height # 11
grid.t <- 4 # Time


# ===================
# = Species Options =
# ===================
ns <- 30 # Number of Species
n0s <- 20


# ======================
# = Simulation Options =
# ======================
n.obs.reps <- 1 # number of time to observe the same true process (each observation is analyzed separately)
n.ss <- 4 # number of substrata (for observation)
n.ss.mu <- trunc((n.ss*grid.w*grid.h)*(50/100)) #max(trunc((n.ss*grid.w*grid.h)/3*2), grid.w*grid.h) # total substrata observed
base.chance <- 1 #plogis(rnorm(ns)) #rbeta(ns,2,2) #runif(n=ns, 0.2, 0.8) # baseline detectability (before ID chance)

# Create chance to be identified if caught
# Can also be used to represent a generic
#  time-varying chance of being detected
obs.chance <- function(dim2=ns, dim1=grid.t, dim3=n.obs.reps, rand.gen=rnorm, chances, ...){
	# If chances is supplied, then each species will have the same chance in each year, and each replicate
	# species may differ tho
	
	rand.gen <- match.fun(rand.gen)
	dots <- list(...)
	# stopifnot(all.same(sapply(dots, length)))
	
	if(missing(chances)){
		
		if(length(dots)>0){
			chances <- matrix(mapply(rand.gen, MoreArgs=list(n=dim2), ...), nrow=dim2)
		}else{
			chances <- matrix(rand.gen(n=dim2,...))
		}
		
		v.rep <- function(...){rep(c(...),each=max(floor(dim1/ncol(chances)),1))}
		
		t.lvls <- unlist(apply(chances, 1, function(x)list(v.rep(x))),F,F)
		# t.noID <- aperm(simplify2array(lapply(t.lvls, roll.recycle, dim3, dim1)), c(2, 3, 1))
		# t.noID0 <- lapply(t.lvls, roll.recycle, dim3, dim1)
		# t.noID <- aperm(array(simplify2array(t.noID0), dim=c(dim3,dim1,dim2)), c(2, 3, 1))
		# t.noID
		
	}else{
		stopifnot(dim2%%length(chances)==0)
		chances <- rep(chances, each=dim2%/%length(chances))
		
		t.lvls <- lapply(chances, c)
		# t.noID <- aperm(simplify2array(lapply(t.lvls, roll.recycle, dim1, dim3)), c(1, 3, 2))
		# t.noID0 <- lapply(t.lvls, roll.recycle, dim3, dim1)
		# t.noID <- aperm(array(simplify2array(t.noID0), dim=c(dim3,dim1,dim2)), c(2, 3, 1))
		# t.noID
	}
	
	t.noID0 <- lapply(t.lvls, roll.recycle, dim3, dim1)
	t.noID <- aperm(array(simplify2array(t.noID0), dim=c(dim3,dim1,dim2)), c(2, 3, 1))
	t.noID	
}

# obs.chance 

# the main use of the function will be when chances is not supplied, 
# and the function performs the random number generation
# In a given year, all species will be drawn from the same distribution
# That distribution can change among years by supplying vectors to ... . Between replicates, 
# the order of which years correspond to which distribution changes

# dim2 is number of species; chances/ rng's always differ across dim2 (unless chances supplied with 1 unique)
# dim1 is time steps; chances will only differ across dim1 if length of arguments in ... are > 1
# dim3 is replicates; chances will be same across dim3, but dim1 shift by 1 index between each dim3
# rand.gen is a function for random number generation; its first argument must be \code{n}. Used to generate chances when chances is not specified.
# chances are predetermined probabilities that can optionally be provided for each species. length(chances) must be a multiple of dim2. If chances is supplied, rand.gen will not be used.
# ... arguments to be passed to rand.gen

# ===========================
# = Examples for obs.chance =
# ===========================

# species (dim2) are drawn from same distribution, but each is different
# no changes between years (dim1), thus no changes between replicates(dim3)
# plogis(obs.chance(dim2=6, dim3=2, dim1=7))

t.noID.mus <- c(-0.5, 0.5)
t.noID.sd <- 1
t.noID <- plogis(obs.chance(dim2=ns, dim1=grid.t, dim3=n.obs.reps, mean=t.noID.mus, sd=t.noID.sd))


# =================================
# = Do Simulation of True Process =
# =================================
# Simulate environment
# env <- sim.env(grid.w=grid.w, grid.h=grid.h, grid.t=grid.t, X.slope=0.75*(12/grid.t))
X.slope <- 0
env <- sim.env(grid.w=grid.w, grid.h=grid.h, grid.t=grid.t, X.slope=X.slope)

# Simulate Species
out <- sim.spp.proc(env, ns=ns, dynamic=FALSE)

# name output attributes for easy access
spp.bio <- attr(out, "spp.bio")
grid.X <- attr(out, "grid.X")
S.dens.X <- attr(out, "spp.densX")
dims <- attr(out, "dims")

# Get S, a list of bricks of length grid.t specifying species presence
S <- getS(out)


# =========================================
# = Observation Level and Format for MSOM =
# =========================================
# The loop is for re-observing the same true process multiple times
for(i in 1:n.obs.reps){
	if(i==1){
		out.obs <- obs.spp(out, n.ss, n.ss.mu, base.chance=rep(1,ns), t.noID[,,i])
		formatted <- spp2msom(out.obs)
		new.simDat <- formatted$simDat 
		simCov <- formatted$simCov 
		simCov.NA <- formatted$simCov.NA 
		simCov.precs <- formatted$simCov.precs 
		simCov.precs.bad <- formatted$simCov.precs.bad 
		sim.cov.names <- formatted$sim.cov.names
		
		big.out.obs <- list(out.obs)
		
		names(new.simDat) <- paste(names(new.simDat),i, sep=".")
		big.simDat <- new.simDat
	}else{
		big.out.obs[[i]] <- obs.spp(out, n.ss, n.ss.mu, base.chance=rep(1,ns), t.noID[,,i])
		new.simDat <- spp2msom(big.out.obs[[i]])$simDat
		names(new.simDat) <- paste(names(new.simDat),i, sep=".")
		big.simDat <- c(big.simDat, new.simDat)
	}
}


# ========================================
# = Analyze simData with Stan msomStatic =
# ========================================
spp2msom_2dt <- function(big.simDat, simCov){
	test.spp <- data.table(reshape2:::melt.list(big.simDat))

	test.spp[,c("year","n.obs.reps"):=list(as.numeric(gsub("year([0-9]{1,2})\\.[0-9]{1,2}","\\1",L1)), as.numeric(gsub("year([0-9]{1,2})\\.([0-9]{1,2})","\\2",L1)))]
	test.spp[,L1:=NULL]
	setnames(test.spp, "value", "abund")
	setkey(test.spp, year,stratum,K)

	test.cov1 <- data.table(reshape2:::melt.list(simCov))
	setnames(test.cov1, c("value","L1"), c("bt","year"))
	test.cov1[,stratum:=test.spp[,unique(stratum)]]
	setcolorder(test.cov1, c('year','stratum','bt'))
	setkey(test.cov1, year, stratum)

	test <- test.spp[test.cov1]
	
	return(test)
}
sim_dt <- spp2msom_2dt(big.simDat, simCov)

# aggregate and transform (^2) btemp
mk_cov_rv_pow(sim_dt, "bt", across="K", by=c("stratum","year"), pow=2)

# scale and aggregate doy
sim_dt[,yr:=as.numeric(year)]
sim_dt[,obs_rep:=as.numeric(n.obs.reps)]
mk_cov_rv(sim_dt, "yr", across="K", by=c("stratum","year"))

sim_dt <- sim_dt[n.obs.reps==1]

# ---- Get Basic Structure of MSOM Data Input ----
staticData <- msomData(Data=sim_dt, n0=n0s, cov.vars=c(bt="bt", bt2="bt2",yr="yr"), u.form=~bt+bt2, v.form=~yr, valueName="abund", cov.by=c("year","stratum"))

staticData_sd <- msomData(Data=sim_dt, n0=1, cov.vars=c(bt_sd="bt_sd",yr_sd="yr_sd",bt2_sd="bt2_sd"), u.form=~bt_sd+bt2_sd-1, v.form=~yr_sd-1, valueName="abund", cov.by=c("year","stratum"))[c("U","V")]


# ---- Split Covariates into Constants and Random Variables ----
getCovType <- function(UV, type=c("constant","mu","sd")){
	
	# Dimension Sizes, Number, and Names
	dims <- dim(staticData[[UV]])
	nD <- length(dims)
	dn <- dimnames(staticData[[UV]])[[nD]]
	dn_sd <- dimnames(staticData_sd[[UV]])[[nD]]
	
	# Get Names of Desired Covariates
	covNames <- switch(type,
		constant = dn[!dn%in%gsub("_sd","",dn_sd)],
		mu = dn[dn%in%gsub("_sd","",dn_sd)],
		sd = dn_sd
	)
	
	# Get the Desired Covariates
	covOut <- switch(type,
		constant = staticData[[UV]][,,covNames],
		mu = staticData[[UV]][,,covNames],
		sd = staticData_sd[[UV]][,,covNames]
	)
	
	# Convert to Array, Get Dimension Information
	covOut <- as.array(covOut)
	cO_dims <- dim(covOut)
	cO_nD <- length(cO_dims)
	
	# Make Sure covOut Dimensions match Full Cov Dimensions
	if(cO_nD<nD){
		# if here, I'm guessing it's because length(covNames) is 1
		# and in that case, I want an array of the orignal major dims to be returned
		stopifnot(length(covNames)==1)
		
		cO_newDim <- c(cO_dims,rep(1,nD-cO_nD)) 
		covOut <- array(covOut, dim=cO_newDim, dimnames=c(dimnames(covOut),covNames))
	}
	
	# Return
	return(covOut)
}

# ---- Do the Splitting Function ----
for(UV in c("U","V")){
	for(type in c("constant","mu","sd")){
		switch(type,
			constant={staticData[[paste0(UV,"_c")]] <- getCovType(UV,type)},
			{staticData[[paste(UV,type,sep="_")]] <- getCovType(UV,type)}
		)
	}
}

# ---- Add Sizes of Constant and RV U/V Arrays ----
staticData$nU_rv <- dim(staticData$U_mu)[3]
staticData$nV_rv <- dim(staticData$V_mu)[3]
staticData$nU_c <- dim(staticData$U_c)[3]
staticData$nV_c <- dim(staticData$V_c)[3]

stopifnot(staticData$nV == staticData$nV_c + staticData$nV_rv)
stopifnot(staticData$nU == staticData$nU_c + staticData$nU_rv)

# ---- Add a counter for nJ (number of sites in each year) ----
staticData$nJ <- apply(staticData$nK, 1, function(x)sum(x>0))

# ---- Aggregate abundance across samples ----
staticData$X <- apply(staticData$X, c(1,2,4), function(x)sum(x))


# =====================
# = Fit Model in Stan =
# =====================
library(rstan)
model_file <- "trawl/trawlDiversity/inst/stan/msomStatic.stan"

sim_msom <- stan(
	file=model_file, 
	data=staticData, 
	control=list(stepsize=0.01, adapt_delta=0.95, max_treedepth=15),
	chains=4, iter=150, refresh=1, seed=1337, cores=4, verbose=F
)


# ==============================
# = Get True and Estimated Psi =
# ==============================
# ---- True Psi ----
use.logit.psi <- FALSE
agg.psi <- FALSE
dim.conv1 <- c(grid.w*grid.h, ns, grid.t, n.obs.reps)

get.psiTrue <- function(x, use.logit=FALSE, agg=FALSE){
	psi.true <- attr(x, "psi")
	psi.true <- array(psi.true, dim=dim.conv1)
	if(use.logit){
		psi.true <- logit(pmax(pmin(psi.true,1-1E-3),1E-3))
	}
	if(agg){
		psi.true <- apply(psi.true, c(1,2,3), mean)
	}
	return(psi.true)
}
man.psi.true <- attributes(big.out.obs[[1]])$psi
psi.true <- get.psiTrue(big.out.obs[[1]], use.logit.psi, agg.psi)
# psi.true <- aperm(apply(psi.true, c(1,2,3), mean), c(3,1,2))
psi.true <- aperm(psi.true[,,,1], c(3,1,2))

# ---- Estimated Psi ----
sims <- rstan::extract(sim_msom)
psi_mean <- plogis(apply(sims$logit_psi, 2:4, mean))

plot(psi.true, psi_mean[,,1:30])

# ================================
# = Get True and Estimated Theta =
# ================================



