#' Process Output from msomStatic (annual)
#' 
#' Processes a list structured as run_msom output (level 2) within a year (level 1) to be used to summarize diversity. Note that this function then only processes 1 region at a time, so it expects a list whose length is equal to the number of years for a region. Currently, intended to work with Stan model.
#' 
#' @param reg_out A list with length equal to number of years in a region, with each element containing output from run_msom
#' @param save_mem Save memory be deleting intermediate steps as function progresses; default TRUE (only reason for FALSE is for debugging)
#' @param obs_yrs Optional vector of integers of calendar years used to process richness from observed data; used to subset years of MSOM analysis
#' 
#' @details Right now only intended for use with specific structuring of the output, so that it matches the output expected from running each year separately using the Stan version of the msomStatic model.
#' 
#' @export
process_msomStatic <- function(reg_out, save_mem=TRUE, obs_yrs){
	
	# ====================
	# = Organize Read-in =
	# ====================
	out <- lapply(reg_out, function(x)x$out)
	empty_ind <- sapply(out, is.null)
	out <- out[!empty_ind]
	inputData <- lapply(reg_out, function(x)x$inputData)[!empty_ind]
	info <- lapply(reg_out, function(x)x$info)[!empty_ind]
	
	regs <- sapply(info, function(x)x["reg"])
	stopifnot(lu(regs)==1)
	reg <- unique(regs)
	
	langs <- unlist(sapply(info, function(x)x["language"]))
	stopifnot(lu(langs)==1)
	lang <- unique(langs)
	
	
	# =====================
	# = Get Full Data Set =
	# =====================
	# ---- Makes [rd] ----
	info_yrs <- sapply(info, function(x)as.integer(x['year']))
	sub_reg <- reg
	# rd <- data_all[year %in% info_yrs & sub_reg==(reg)] # data_all is an object associated with the trawlDiversity package!
	# rd_yr <- rd[,sort(unique(year))]
	
	
	# ================================================
	# = Subset MSOM Output to Years in Full Data Set =
	# ================================================
	if(!missing(obs_yrs)){
		yr_match_msomSub <- info_yrs %in% obs_yrs
		out <- out[yr_match_msomSub]
		inputData <- inputData[yr_match_msomSub]
		info <- info[yr_match_msomSub]
		info_yrs <- info_yrs[yr_match_msomSub]
	
		# stopifnot(all(info_yrs==rd_yr))
	}

	
	
	# ======================================
	# = Memory Save by Deleting from 'out' =
	# ======================================
	if(save_mem){
		# rm(list="reg_out")
		if(lang == "JAGS"){
			for(i in 1:length(out)){
				# out[[i]]$BUGSoutput[c("sims.list","sims.matrix")] <- NULL
				out[[i]]$BUGSoutput[c("sims.matrix")] <- NULL
			}
		}
	}
	
	

	
	
	# ================================================
	# = Get Colonization Patterns from Observations =
	# ================================================
	# ---- Makes [colonization] ----
	# colonization <- get_colonizers(d=rd)
	
	
	# ==============
	# = Covariates =
	# ==============
	# strat2lld <- function(x){
	# 	s <- strsplit(x, split=" ")
	# 	lon <- sapply(s, function(x)x[1])
	# 	lat <- sapply(s, function(x)x[2])
	# 	depth_interval <- sapply(s, function(x)x[3])
	# 	data.table(lon=as.numeric(lon), lat=as.numeric(lat), depth_interval=as.numeric(depth_interval))
	# }
	#
	# # ---- Bottom Temperature ----
	# # ---- Makes [bt] ----
	# get_bt <- function(id){
	# 	(id$U[,"bt"]*id$scaling["btemp.sd"]) + id$scaling["btemp.mu"]
	# }
	# bt0 <- lapply(inputData, get_bt)
	# bt2dt <- function(x,y)data.table(stratum=names(x), bt=x, year=y)
	# bt <- rbindlist(mapply(bt2dt, bt0, info_yrs, SIMPLIFY=FALSE))
	#
	# bt[,c("lon","lat","depth_interval"):=strat2lld(stratum)]
	# bt[,bt_col:=zCol(256, bt)]
	#
	# # ---- Depth ----
	# # ---- Adds to [bt] ----
	# get_depth <- function(id){
	# 	(id$U[,"depth"]*id$scaling["depth.sd"]) + id$scaling["depth.mu"]
	# }
	# depth0 <- lapply(inputData, get_depth)
	# depth2dt <- function(x,y)data.table(stratum=names(x), depth=x, year=y)
	# depth <- rbindlist(mapply(depth2dt, depth0, info_yrs, SIMPLIFY=FALSE))
	#
	# depth[,c("lon","lat","depth_interval"):=strat2lld(stratum)]
	# depth[,depth_col:=zCol(256, depth)]
	#
	# bt <- merge(bt, depth, by=c("stratum","lon","year","lat","depth_interval"), all=TRUE)
	
	
	# ====================================================
	# = Get Posterior Samples (iterations) of Parameters =
	# ====================================================
	# ---- Hyperparameters and Omega ----
	# ---- Makes [param_iters] ----
	if(lang == "Stan"){
		pars_trace <- c("Omega","alpha_mu[1]", "alpha_mu[2]", "alpha_mu[3]", "alpha_mu[4]", "alpha_mu[5]", "beta_mu[1]")
	}else{
		pars_trace <- c("Omega","alpha_mu[1]", "alpha_mu[2]", "alpha_mu[3]", "alpha_mu[4]", "alpha_mu[5]", "beta_mu")
	}
	
	param_iters <- list()
	for(i in 1:length(out)){
		param_iters[[i]] <- data.table(reg=reg, year=info_yrs[i], get_iters(out[[i]], pars_trace, lang))
	}
	param_iters <- rbindlist(param_iters)
	
	# ---- Species-specific Alpha and Beta Parameters ----
	# ---- Makes [ab] ----
	# Only for observed species (i.e., parameters that I can tie to a Latin name)
	ab_all <- mapply(get_ab, inputData, out, info_yrs, SIMPLIFY=FALSE)
	ab <- rbindlist(ab_all)
	
	# ---- Unscale Alpha ----
	# ---- Makes [alpha_unscale] ----
	unscale_ab <- function(abDat){
		ab_un <- unscale(
			abDat[ab_ind==1, value],
			abDat[ab_ind==2, value],
			abDat[ab_ind==3, value],
			abDat[ab_ind==4, value],
			abDat[ab_ind==5, value],
			abDat[ab_ind==1,btemp.mu],
			abDat[ab_ind==1,depth.mu],
			abDat[ab_ind==1,btemp.sd],
			abDat[ab_ind==1,depth.sd]	
		)
		return(as.data.table(ab_un))
	}
	scaling <- rbindlist(lapply(inputData, function(x)as.list((x$scaling))))
	scaling[,year:=as.numeric(sapply(info, function(x)x["year"]))]
	ab2 <- merge(ab, scaling, by="year")
	alpha_unscale <- ab2[par=="alpha", unscale_ab(.SD),by=c("spp","spp_id","year")]
	ab <- ab[,list(value_mu=mean(value), value_sd=sd(value)),by=c("parameter","par","ab_ind","spp_id","spp","year")]
	rm(list="ab2")
	
	
	
	# =========================================
	# = Get Estimates of Region-Wide Richness =
	# =========================================
	# ---- Makes [processed] ----
	
	# ---- First, psi and theta (Stan), then X_obs (Stan and JAGS) ----
	if(lang == "Stan"){
		# ---- Species Detection Probs ----
		# first dimension is iterations
		# second is year (which is just 1)
		# third is site (stratum; but detection same across)
		# fourth is species
		theta_dist <- lapply(out, function(x)plogis(extract(x, pars="logit_theta")[[1]][,1,,]))
		theta_dist_reg <- lapply(theta_dist, function(x)x[,1,])
	
		# ---- Species Presence Probs ----
		# same dimensions as theta, except in this model presence changes across strata
		psi_dist <- lapply(out, function(x)plogis(extract(x, pars="logit_psi")[[1]][,1,,]))
		
		# ---- Observed Presences (Stan) ----
		X_obs <- lapply(inputData, function(x)x$X[1,,])
		
	}else if(lang == "JAGS"){
		# ---- Observed Presences (JAGS) ----
		# (theta and psi aren't needed if JAGS)
		X_obs <- lapply(inputData, function(x)x$X)
	}
	
	# ---- Richness in the Region ----
	Omega_iter <- param_iters[,list(year,Omega)] #lapply(out, get_iters, pars="Omega", lang="JAGS")
	Omega_mean <- Omega_iter[,mean(Omega), by="year"][,V1] #sapply(Omega_iter, function(x)x[,mean(Omega)])
	# naive_rich <- sapply(inputData, function(x)x$N)
	rich_pureModel <- Omega_mean * sapply(inputData, function(x)x$nS)
	
	if(lang=="Stan"){
		reg_pres_dist <- lapply(psi_dist, function(x)apply(x, c(1,3), function(x)(1-prod(1-x))))
		reg_pres <- lapply(reg_pres_dist, function(x)colMeans(x)) # see note below when I do ~this for JAGS --- I think this is kinda wrong, because I'm averaging over posterior iterations for each species, rather than getting the region-wide richness for each iteration then taking the average of richness over iterations
		reg_rich <- sapply(reg_pres, sum)
		
	}else if(lang == "JAGS"){
		z_spp <- function(x){gsub("Z\\[[0-9]+\\,([0-9]+)\\]$", "\\1", x)}
		z_j <- function(x){gsub("Z\\[([0-9]+)\\,[0-9]+\\]$", "\\1", x)}
		reg_rich <- rep(NA, length(out))
		for(i in 1:length(out)){
			t_Z_big <- get_iters(out[[i]], pars="Z", lang="JAGS")
			t_Z_big[,iter:=(1:nrow(t_Z_big))]
				
			t_Z_big_long <- data.table:::melt.data.table(t_Z_big, id.vars=c("iter","chain"))
			t_Z_big_long[, c("sppID","jID"):=list(z_spp((variable)), z_j((variable))), by=c("variable")]
			if(save_mem){rm(list="t_Z_big")}
			
			mu_site_Z <- t_Z_big_long[,j={ list(mu_site_Z = max(value)) }, by=c("iter","chain","sppID")]
			reg_rich_iter <- mu_site_Z[, j={list(reg_rich = sum(mu_site_Z))}, by=c("iter","chain")]
			reg_rich[i] <- reg_rich_iter[,mean(reg_rich)]
			if(save_mem){rm(list="t_Z_big_long")}
				
			 # this is different than how I did it for Stan -- for Stan I calculated took the mean (across iterations) probability of a species being present somewhere in the region, then summed up across species to get richness. Here I get the pres/abs for each species in the region for each iteration, then I sum across species but w/in an iteration to get the posterior of region wide richness (that summing is done on this line), and then in the next line I take the mean of the posterior distribution of region-wide richness
		}
	}
	
	# create processed object
	processed <- data.table(reg = reg, year=info_yrs, Omega=Omega_mean, reg_rich=reg_rich)
	# processed <- merge(processed, colonization$n_cep, by="year", all=TRUE)
	# processed <- merge(processed, bt[,list(bt_ann=mean(bt)), by="year"], by="year", all=TRUE)
	

	return(list(param_iters=param_iters, processed=processed, ab=ab, alpha_unscale=alpha_unscale))
	
}
	
	
	
	