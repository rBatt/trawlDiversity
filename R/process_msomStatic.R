#' Process Output from msomStatic (annual)
#' 
#' Processes a list structured as run_msom output (level 2) within a year (level 1) to be used to summarize diversity. Note that this function then only processes 1 region at a time, so it expects a list whose length is equal to the number of years for a region. Currently, intended to work with Stan model.
#' 
#' @param reg_out A list with length equal to number of years in a region, with each element containing output from run_msom
#' @param save_mem Save memory be deleting intermediate steps as function progresses; default TRUE (only reason for FALSE is for debugging)
#' 
#' @details Right now only intended for use with specific structuring of the output, so that it matches the output expected from running each year separately using the Stan version of the msomStatic model.
#' 
#' @export
process_msomStatic <- function(reg_out, save_mem=TRUE){
	
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
		
	if(save_mem){
		# rm(list="reg_out")
		if(lang == "JAGS"){
			for(i in 1:length(out)){
				# out[[i]]$BUGSoutput[c("sims.list","sims.matrix")] <- NULL
				out[[i]]$BUGSoutput[c("sims.matrix")] <- NULL
			}
		}
	}
	
	
	# =====================
	# = Get Full Data Set =
	# =====================
	# ---- Makes [rd] ----
	rd <- data_all[reg==(reg)] # data_all is an object associated with the trawlDiversity package!
	rd_yr <- rd[,sort(unique(year))]
	
	
	# ================================================
	# = Get Colonization Patterns from Observations =
	# ================================================
	# ---- Makes [colonization] ----
	colonization <- get_colonizers(d=rd)
	
	
	# ==============
	# = Covariates =
	# ==============
	# ---- Bottom Temperature ----
	# ---- Makes [bt] ----
	bt0 <- lapply(inputData, function(x)x$U[,"bt"])	
	bt2dt <- function(x,y)data.table(stratum=names(x), bt=x, year=y)
	bt <- rbindlist(mapply(bt2dt, bt0, rd_yr, SIMPLIFY=FALSE))
	strat2lld <- function(x){
		s <- strsplit(x, split=" ")
		lon <- sapply(s, function(x)x[1])
		lat <- sapply(s, function(x)x[2])
		depth_interval <- sapply(s, function(x)x[3])
		data.table(lon=as.numeric(lon), lat=as.numeric(lat), depth_interval=as.numeric(depth_interval))
	} 
	bt[,c("lon","lat","depth_interval"):=strat2lld(stratum)]
	bt[,bt_col:=zCol(256, bt)]
	
	
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
		param_iters[[i]] <- data.table(reg=reg, year=rd_yr[i], get_iters(out[[i]], pars_trace, lang))
	}
	param_iters <- rbindlist(param_iters)
	
	# ---- Species-specific Alpha and Beta Parameters ----
	# ---- Makes [ab] ----
	# Only for observed species (i.e., parameters that I can tie to a Latin name)
	ab_all <- mapply(get_ab, inputData, out, rd_yr, SIMPLIFY=FALSE)
	ab <- rbindlist(ab_all)
	
	
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
	naive_rich <- sapply(inputData, function(x)x$N)
	rich_pureModel <- Omega_mean * sapply(inputData, function(x)x$nS)
	
	if(lang=="Stan"){
		reg_pres_dist <- lapply(psi_dist, function(x)apply(x, c(1,3), function(x)(1-prod(1-x))))
		reg_pres <- lapply(reg_pres_dist, function(x)colMeans(x)) # see note below when I do ~this for JAGS --- I think this is kinda wrong, because I'm averaging over posterior iterations for each species, rather than getting the region-wide richness for each iteration then taking the average of richness over iterations
		reg_rich <- sapply(reg_pres, sum)
		
	}else if(lang == "JAGS"){
		Z_big <- lapply(out, get_iters, pars="Z", lang="JAGS")
		nr <- nrow(Z_big[[1]])
		for(i in 1:length(Z_big)){
			Z_big[[i]] <- Z_big[[i]][,iter:=(1:nr)]
		}
		
		Z_big_long <- lapply(Z_big, data.table:::melt.data.table, id.vars=c("iter","chain"))
		z_spp <- function(x){gsub("Z\\[[0-9]+\\,([0-9]+)\\]$", "\\1", x)}
		z_j <- function(x){gsub("Z\\[([0-9]+)\\,[0-9]+\\]$", "\\1", x)}
		for(i in 1:length(Z_big_long)){
			td <- Z_big_long[[i]]
			td[, c("sppID","jID"):=list(z_spp((variable)), z_j((variable))), by=c("variable")]
			Z_big_long[[i]] <- td
		}
		
		reg_rich <- rep(NA, length(Z_big_long))
		for(i in 1:length(Z_big_long)){
			mu_site_Z <- Z_big_long[[i]][,j={ list(mu_site_Z = max(value)) }, by=c("iter","chain","sppID")]
			reg_rich_iter <- mu_site_Z[, j={list(reg_rich = sum(mu_site_Z))}, by=c("iter","chain")] # this is different than how I did it for Stan -- for Stan I calculated took the mean (across iterations) probability of a species being present somewhere in the region, then summed up across species to get richness. Here I get the pres/abs for each species in the region for each iteration, then I sum across species but w/in an iteration to get the posterior of region wide richness (that summing is done on this line), and then in the next line I take the mean of the posterior distribution of region-wide richness
			reg_rich[i] <- reg_rich_iter[,mean(reg_rich)]
		}

		if(save_mem){rm(list="Z_big_long")}
	}
	
	# get unobserved but present species
	unobs_rich <- reg_rich - naive_rich
	frac_unobs_rich <- unobs_rich/reg_rich
	
	# create processed object
	processed <- data.table(reg = reg, year=rd[,sort(una(year))], Omega=Omega_mean, reg_rich=reg_rich, naive_rich=naive_rich, unobs_rich=unobs_rich)
	processed <- merge(processed, get_colonizers(rd)$n_cep, by="year", all=TRUE)
	processed <- merge(processed, bt[,list(bt_ann=mean(bt)), by="year"], by="year", all=TRUE)
	

	return(list(rd=rd, colonization=colonization, bt=bt, param_iters=param_iters, processed=processed, ab=ab))
	
}
	
	
	
	