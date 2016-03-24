#' Run an MSOM Model
#' 
#' Run a multispecies occupancy model on trawl data. The model is multi-year, and there are two basic versions; first one is a 'static' model for which no temporal process is explicitly modelled (the years are 'stacked'). The second is a 'dynamic' model that includes persistence and colonization parameters as processes that facilitate the state transition between years. Each of these models can be run in either JAGS or Stan. Both models accept covariates for the detection and presence processes, and both sets of covariates can be modelled as random variables, constants (varying, but known precisely), or a mixture of the two.
#' 
#' @param reg character, region name
#' @param regX.a1 optional data set formatted according to output from \code{\link{trim_msom}}; if missing, gets a data set using region name and defaults on trimming and aggregating
#' @param params_out character vector specifying category of parameters to report in output; if custom, points to \code{custom_params} as well
#' @param custom_params possible characters specifying custom parameters, potentially indexed; see Details 
#' @param model_type type of model. Currently both are multiyear. 'Dynamic' includes terms for persistence and colonization (species-specific), whereas 'Static' simply 'stacks' the years.
#' @param n0 integer indicating the number of species to use in data augmentation; defaults to 50, is passed to \code{\link{msomData}}
#' @param chains integer number of chains, default is 4
#' @param cores integer number of cores for parallel processes; default is half of avaialble cores
#' @param iter number of iterations; the default depends on \code{language} and on /code{test}; testing in Stan is 50, non-testing is 200. Testing in JAGS is 500, non-testing is 5000.
#' @param thin Integer, thinning rate. Default is set to yield 200 draws from the posterior, or if that is not possible, thinning rate is 1 and draws will be \code{iter}/2.
#' @param language character indicating the language to be used --- JAGS or Stan
#' @param test Logical, whether to do this run as a 'test' run. Default is FALSE. If TRUE, fewer iterations are run (unless overridden by non-default), and the data set is subsetted 
#' @param test_sub a named list with elements 'stratum', 'year', and 'spp'. Each element should be an integer indicating the number of levels to select for each of those dimensions. Used for subsetting when \code{test} is TRUE.
#' @param seed integer random number seed
#' @param pre_save Logical; if TRUE (default) saves a workspace image before running the model
#' @param save_dir Character string indicating the location of the directory to save the intermediate image; default is current directory
#' @param model_dir Character string indicating the location of the model file; default is selected automatically based on \code{language} and \code{model_type}, and looks to models that come with this package
#' @param compiled_model Only used when language="Stan"; a Stan model compiled using \code{\link{stan_model}}. Useful for preventing repeated model recompilation and/or reapeated loading of DLLs, which can cause an error if done enough.
#' @param ... arguments passed to \code{rstan::sampling}
#' 
#' @details
#' Both \code{params_out} and \code{custom_params} must find matches in the output of \code{\link{msom_params}}. For all parameters, use 'params'; main-effect parameters specified via 'params_main'; random-effect parameters via 'params_random'. Latent stochastic nodes/ parameters via 'params_latent'. Additional flexibility offered by specifying 'custom', which will add manually specified parameters from \code{custom_params}.
#' 
#' @return
#' Returns a named list of length two. The first element, 'out', contains a JAGS or Stan model object. The second element is a named character vector containing potential file names or prefixes to file names. These names contain collapsed information related to the model run settings, time, file paths, etc.
#' 
#' @import trawlData
#' 
#' @export
run_msom <- function(reg = c("ai", "ebs", "gmex", "goa", "neus", "newf", "ngulf", "sa", "sgulf", "shelf", "wcann", "wctri"), regX.a1, params_out=c("params","params_main","params_random","params_latent","custom"), custom_params=NULL, model_type=c("Dynamic", "Static"), n0=50, chains=4, cores=parallel::detectCores()/2, iter, thin=max(1, floor((iter/2)/200)), language=c("JAGS", "Stan"), test=FALSE, test_sub=list(stratum=4, year=3, spp=10), seed=1337, pre_save=FALSE, save_dir=".", model_dir=file.path(system.file(package="trawlDiversity"), tolower(language)), compiled_model=NULL, ...){
	
	model_type <- match.arg(model_type)
	language <- match.arg(language)
	reg <- match.arg(reg)
	params_out <- match.arg(params_out, several.ok=TRUE)
	
	
	requireNamespace("parallel", quietly=TRUE)
	requireNamespace("trawlData", quietly=TRUE)
	requireNamespace("rbLib", quietly=TRUE)
	if(language=="Stan"){requireNamespace("rstan", quietly=TRUE)}
	if(language=="JAGS"){requireNamespace("R2jags", quietly=TRUE)}
	if(language=="JAGS"){library("R2jags")}
		
	if(missing(iter)){
		if(test){
			iter <- ifelse(language=="Stan", 50, 500)
		}else{
			iter <- ifelse(language=="Stan", 200, 5E3)
		}
	}

	
	
	# check files/ paths
	if(!file.exists(save_dir)){
		stop("save directory (", save_dir, ") does not exist")
	}
	if(!file.exists(model_dir) & ((is.null(compiled_model) & language=="Stan") | language=="JAGS")){
		stop("model directory (", model_dir, ") does not exist")
	}


	# ======================
	# = Subset Data to Use =
	# ======================
	if(missing(regX.a1)){
		regX.a1 <- trim_msom(reg, gridSize=1, plot=FALSE)
	}
	

	if(test){
		set.seed(seed)
		ind <- mpick(regX.a1, p=c(stratum=test_sub$stratum, year=test_sub$year), weight=TRUE, limit=60)
		logic <- expression(
			spp%in%spp[ind]
			& stratum%in%stratum[ind]
			& year%in%year[ind]
		)
		regX.a1 <- regX.a1[eval(logic)][pick(spp, test_sub$spp, w=FALSE)]
	}
	
	regX.a2 <- regX.a1[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, depth=depth, doy=yday(datetime))]


	# ==================
	# = Set Covariates =
	# ==================
	cov.vars_use <- c(bt="bt",bt2="bt2",yr="yr", yr2="yr2", doy="doy", doy2="doy2", depth="depth", depth2="depth2")


	# ======================================
	# = Aggregate and Transform Covariates =
	# ======================================
	# rename columns for shorthand
	setnames(regX.a2, c("btemp"), c("bt"))
	# setnames(regX.a2, c("stemp"), c("st"))
	# regX.a2[,yr:=scale(as.integer(year))] # fails with 1 year; also, mu and sd not weighted to unique years, so kinda weird
	regX.a2[,yr:=as.numeric(year)]
	regX.a2[,year:=as.character(year)]

	# aggregate and transform (^2) btemp stemp
	btemp.mu <- regX.a2[,mean(bt, na.rm=TRUE)]
	btemp.sd <- regX.a2[,sd(bt, na.rm=TRUE)]
	if(btemp.sd == 0){btemp.sd <- 1}
	regX.a2[,bt:=(bt-btemp.mu)/btemp.sd]
	mk_cov_rv_pow(regX.a2, "bt", across="K", by=c("stratum","year"), pow=2)
	# mk_cov_rv_pow(regX.a2, "st", across="K", by=c("stratum","year"), pow=2)
	
	# make smallest yr 0, then aggregate
	regX.a2[,yr:=(yr-min(yr, na.rm=TRUE))]
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
	# inputData <- msomData(Data=regX.a2, n0=n0, cov.vars=cov.vars_use, u.form=~bt+bt2, v.form=~year+doy, valueName="abund", cov.by=c("year","stratum"), u_rv=c("bt","bt2"), v_rv=c("doy"))
	inputData <- msomData(Data=regX.a2, n0=n0, cov.vars=cov.vars_use, u.form=~bt+bt2+depth+depth2, v.form=~1, valueName="abund", cov.by=c("year","stratum"))

	inputData$nJ <- as.array(apply(inputData$nK, 1, function(x)sum(x>0))) # number of sites in each year
	inputData$X <- apply(inputData$X, c(1,2,4), function(x)sum(x)) # agg abund across samples
	
	
	# ==========================================================
	# = Check parameter names against model_type and grid size =
	# ==========================================================
	if("custom"%in%params_out){
		if(length(params_out)>1){
			params_out <- params_out[params_out!="custom"]
		}
		if(any(grepl("\\[", custom_params))){
			mps <- do.call(msom_params, c(language=language, model_type=model_type, inputData[c("nT","Jmax","nU","nV","nS")], noIndex=FALSE))[-1]
			mps_keep <- unlist(mps[params_out], use.names=FALSE)
			mps_custom <- unlist(mps, use.names=FALSE)
			custom_params <- match.arg(custom_params, choices=mps_custom, several.ok=TRUE)
			mps_keep <- unique(ifInd_strip_noInd(unique(c(mps_keep,custom_params))))
		}else{
			mps <- do.call(msom_params, c(language=language, model_type=model_type, inputData[c("nT","Jmax","nU","nV","nS")], noIndex=TRUE))[-1]
			mps_keep <- unlist(mps[params_out], use.names=FALSE)
			mps_custom <- unlist(mps, use.names=FALSE)
			mps_custom <- mps_custom[mps_custom%in%custom_params]
			custom_params <- match.arg(custom_params, choices=mps_custom, several.ok=TRUE)
			mps_keep <- unique(c(mps_keep,custom_params))
		}
	}else{
		mps_keep <- unlist(do.call(msom_params, c(language=language, model_type=model_type, inputData[c("nT","Jmax","nU","nV","nS")], noIndex=TRUE))[params_out], use.names=FALSE)
	}
	
	
	# =============
	# = Kill NA's =
	# =============
	# replace NA mu's with 0
	# replace NA sd's with 1E3
	na_to_0 <- function(x){x[is.na(x)] <- 0; x}
	na_to_big <- function(x, big=1E3){x[is.na(x)] <- big; x}
	
	inputData$U_mu <- na_to_0(inputData$U_mu)
	inputData$V_mu <- na_to_0(inputData$V_mu)
	inputData$U_sd <- na_to_big(inputData$U_sd)
	inputData$V_sd <- na_to_big(inputData$V_sd)
	
	stopifnot(!any(sapply(inputData, function(x)any(is.na(x)))))
	
	if(language=="JAGS"){
		# change 0 sd's to tiny sd's, so inverse is finite
		zero_to_small <- function(x, small=1E-3){x[x==0] <- small; x}
		inputData$U_sd <- zero_to_small(inputData$U_sd)
		inputData$V_sd <- zero_to_small(inputData$V_sd)
		
		# change names from sd to tau
		names(inputData)[names(inputData)%in%c("U_sd","V_sd")] <- c("U_tau", "V_tau")
		
		# convert sd to precision (tau)
		inputData$U_tau <- 1/inputData$U_tau^2
		inputData$V_tau <- 1/inputData$V_tau^2
		
		stopifnot(!any(sapply(inputData, function(x)any(!is.finite(x)))))
		# stopifnot(!any(sapply(inputData, function(x)any(x==0))))
	}
	
	
	# ---- Set run info before removing inputData elements ----
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
	
	
	# ---- Remove inputData elements to make compatible w/ JAGS ----
	# Ensure valid initial values for JAGS
	if(language=="JAGS"){
		make.inits <- function(){
			Z <- inputData$X
			Z[] <- pmin(1, inputData$X)
			Omega <- runif(1, inputData$N/(inputData$nS), 1)
			w <- c(rep(1, inputData$N), rbinom(inputData$nS-inputData$N, size=1, prob=Omega))
			list(Z = Z, Omega = Omega, w = w)	
		}
		inits <- list()
		for(i in 1:chains){
			inits[[i]] <- make.inits()
		}
	}
	
	# Remove Unused list elements
	if(inputData$nU_rv == 0 & inputData$nV_rv == 0){
		inputData$U_c <- NULL
		inputData$U_mu <- NULL
		inputData$U_sd <- NULL
		inputData$V_c <- NULL
		inputData$V_mu <- NULL
		inputData$V_sd <- NULL
		
		inputData$nU_rv <- NULL
		inputData$nV_rv <- NULL
		inputData$nU_c <- NULL
		inputData$nV_c <- NULL
		
		if(language=="JAGS"){
			inputData$U_tau <- NULL
			inputData$V_tau <- NULL
		}
		
		model_type <- paste0(model_type, "_norv")
		
		if(inputData$nT == 1){
			drop_year_dim <- TRUE
			model_type <- paste0(model_type, "_1yr")
		}else{
			drop_year_dim <- FALSE
		}
		
	}else{
		inputData$U <- NULL
		inputData$V <- NULL
		drop_year_dim <- FALSE
	}

	
	if(drop_year_dim){
		inputData$X <- as.matrix(inputData$X[1,,])
		inputData$U <- as.matrix(inputData$U[1,,])
		inputData$V <- as.matrix(inputData$V[1,,])
		inputData$isUnobs <- as.matrix(inputData$isUnobs[1,,])
		inputData$nK <- inputData$nK[1,]
		
		if(language=="JAGS"){
			for(i in 1:chains){
				inits[[i]]$Z <- as.matrix(inits[[i]]$Z[1,,])
			}
		}

		inputData$nT <- NULL
		inputData$nJ <- NULL
	}
	
	iD_N <- inputData$N # grab this here so that it can be added back in before return (useful value)
	if(language=="JAGS"){
		inputData$nJ <- NULL
		inputData$N <- NULL
		inputData$isUnobs <- NULL # Stan models could be modified to make unncessary for them, too
		inputData$Kmax <- NULL # Stan models could be modified to make unncessary for them, too
	}
	
	
 # ---- Model names for save and model read, etc ----
	model_file <- paste0("msom", model_type, ".", tolower(language))
	model_path <- file.path(model_dir, model_file)
	# model_path = "~/Documents/School&Work/pinskyPost/trawl/trawlDiversity/inst/stan/msomStatic_norv.stan"
	# model_path = "~/Documents/School&Work/pinskyPost/trawl/trawlDiversity/inst/jags/msomStatic_norv.jags"

	tag <- paste0("start", format.Date(Sys.time(),"%Y-%m-%d"))

	save_file <- paste0("msom", model_type, "_", reg, "_", tolower(language), "_", iter/1E3, "kIter_", n0, "nZ_", tag, ".RData")
	save_path <- file.path(save_dir,save_file)

	if(pre_save){
		# save.image(save_path)
		save(list=ls(all.names=TRUE), file=save_path)
	}


	# =====================
	# = Fit Model in Stan =
	# =====================
	if(language=="Stan"){
		if(is.null(compiled_model)){
			compiled_model <- rstan::stan_model(model_path)
		}
		stan_control <- list(stepsize=0.01, adapt_delta=0.99, max_treedepth=15)
		out <- rstan::sampling(
			object=compiled_model,
			data=inputData,
			pars = mps_keep,
			control=stan_control, 
			chains=chains, iter=iter, seed=seed, cores=cores, verbose=F, refresh=100, ...
		)
	}else if(language=="JAGS"){
		
		# mps_keep <- gsub("phi", "Phi", mps_keep)
		
		# out <- R2jags::jags(
		# 	data=inputData,
		# 	inits=inits,
		# 	parameters.to.save=mps_keep,
		# 	model.file=model_path,
		# 	jags.seed=seed,
		# 	n.chains=chains,
		# 	n.iter=iter,
		# 	n.thin=thin
		# 	# ,working.directory=paste0(getwd(),"/","trawl/Scripts/Analysis/JAGS")
		# )
	
	for(i in 1:length(inputData)){
		assign(names(inputData)[i], inputData[[i]])
	}
	
	out <- R2jags::jags.parallel(
		data=names(inputData),
		inits=inits[1],
		parameters.to.save=mps_keep,
		model.file=model_path,
		jags.seed=seed,
		n.chains=chains,
		n.iter=iter,
		n.thin=thin, 
		export_obj_names=c("inits", "mps_keep", "model_path", "seed", "chains", "iter", "thin"),
		envir = environment()
		# export_obj_names=c("iter", "thin")
		# ,working.directory=paste0(getwd(),"/","trawl/Scripts/Analysis/JAGS")
	)
		
		
		
	}
	
	# ---- Add N back in to inputData (was removed in the case of JAGS) ----
	inputData$N <- iD_N
	
	# ---- Add other relevant info to inputData ----
	inputData$scaling <- c(
		"doy.mu" = doy.mu,
		"doy.sd" = doy.sd,
		"depth.mu" = depth.mu,
		"depth.sd" = depth.sd,
		"btemp.mu" = btemp.mu,
		"btemp.sd" = btemp.sd
	)
	
	
	# ---- Return ----
	return(list(out=out, inputData=inputData, info=c(reg=reg, language=language, run_info=run_info, model_file=model_file, model_path=model_path, tag=tag, save_file=save_file, save_path=save_path)))
	
	
	# ---- Helpful Reinstall Code (not actually part of function) ----
	# setwd("~/Documents/School&Work/pinskyPost/trawl")
	# remove.packages('trawlDiversity')
	# library(devtools)
	# install('trawlDiversity', upgrade_dependencies=FALSE)
	
}