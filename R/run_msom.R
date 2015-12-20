#' Run an MSOM Model
#' 
#' Run a multispecies occupancy model on trawl data. The model is multi-year, and there are two basic versions; first one is a 'static' model for which no temporal process is explicitly modelled (the years are 'stacked'). The second is a 'dynamic' model that includes persistence and colonization parameters as processes that facilitate the state transition between years. Each of these models can be run in either JAGS or Stan. Both models accept covariates for the detection and presence processes, and both sets of covariates can be modelled as random variables, constants (varying, but known precisely), or a mixture of the two.
#' 
#' @param reg
#' 


run_msom <- function(reg = c("ai", "ebs", "gmex", "goa", "neus", "newf", "ngulf", "sa", "sgulf", "shelf", "wcann", "wctri"), params_out=c("params","params_main","params_random","params_latent","custom"), custom_params=NULL, model_type=c("Dynamic", "Static"), chains=4, cores=parallel::detectCores()/2, iter=50, language=c("JAGS", "Stan"), test=FALSE, test_sub=list(stratum=12, year=15, spp=20), seed=1337){
	
	model_type <- match.arg(model_type)
	language <- match.arg(language)
	reg <- match.arg(reg)
	params_out <- match.arg(params_out, several.ok=TRUE)
	
	
	requireNamespace("parallel", quietly=TRUE)
	requireNamespace("trawlData", quietly=TRUE)
	requireNamespace("rbLib", quietly=TRUE)
	if(language=="Stan"){requireNamespace("rstan", quietly=TRUE)}
	if(language=="JAGS"){requireNamespace("R2jags", quietly=TRUE)}


	# ======================
	# = Subset Data to Use =
	# ======================
	regX.a1 <- trim_msom(reg, gridSize=1, plot=FALSE)

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
	
	regX.a2 <- regX.a1[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, stemp=stemp, depth=depth, doy=yday(datetime))]


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
	
	
	# Ensure valid initial values for JAGS
	if(language=="JAGS"){
		make.inits <- function(){
			Z <- inputData$X
			Z[] <- 1
			list(Z = Z)	
		}
		z.inits <- list()
		for(i in 1:chains){
			z.inits[[i]] <- make.inits()
		}
	}
	
	# Remove Unused list elements
	inputData$U <- NULL
	inputData$V <- NULL
	inputData$nJ <- NULL
	inputData$isUnobs <- NULL
	inputData$Kmax <- NULL
	if(language=="JAGS"){
		inputData$N <- NULL
	}
	
	
	
	model_file <- paste0("msom", model_type, ".", tolower(language))
	# model_path <- file.path(system.file(package="trawlDiversity"), "inst", tolower(language), model_file)
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


	# =====================
	# = Fit Model in Stan =
	# =====================
	if(language=="Stan"){
		stan_control <- list(stepsize=0.01, adapt_delta=0.95, max_treedepth=15)
		out <- stan(
			file=model_path,
			data=inputData,
			pars = mps_keep,
			control=stan_control,
			chains=chains, iter=iter, seed=seed, cores=cores, verbose=F, refresh=1
		)
	}else if(language=="JAGS"){
		
		# mps_keep <- gsub("phi", "Phi", mps_keep)
		
		out <- jags(
			data=inputData,
			inits=z.inits,
			parameters.to.save=mps_keep,
			model.file=model_path,
			jags.seed=seed,
			n.chains=chains,
			n.iter=iter,
			n.thin=thin
			# ,working.directory=paste0(getwd(),"/","trawl/Scripts/Analysis/JAGS")
		)
		
		
		
	}
	
	
	
}