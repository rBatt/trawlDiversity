#' Process Output from msomStatic (annual)
#' 
#' Processes a list structured as run_msom output (level 2) within a year (level 1) to be used to summarize diversity. Note that this function then only processes 1 region at a time, so it expects a list whose length is equal to the number of years for a region. Currently, intended to work with Stan model.
#' 
#' @param rm_out A list with length equal to number of years in a region, with each element containing output from run_msom
#' @param reg Region name (character)
#' 
#' @details Right now only intended for use with specific structuring of the output, so that it matches the output expected from running each year separately using the Stan version of the msomStatic model.
#' 
#' @export
process_msomStatic <- function(rm_out, reg){
	
	library(rstan)
	library(trawlDiversity)
	library(rbLib)
	
	
	# load("./trawlDiversity/pkgBuild/results/msomStatic_norv_ebs_stan_start2016-02-03_14-57-36_r1-1.RData")
	# load("./trawlDiversity/pkgBuild/results/msomStatic_norv_wctri_stan_start2016-02-03_21-59-14_r1-4.RData")
	# load("./trawlDiversity/pkgBuild/results/msomStatic_norv_wcann_stan_start2016-02-04_01-34-58_r5-5.RData")
	# load("./trawlDiversity/pkgBuild/results/msomStatic_norv_sa_stan_start2016-02-04_07-40-01_r6-7.RData")
	# load("./trawlDiversity/pkgBuild/results/msomStatic_norv_neus_stan_start2016-02-05_13-26-12_r8-8.RData")
	
	# load("./trawlDiversity/pkgBuild/results/msomStatic_norv_ai_stan_start2016-02-09_04-44-20_r1-2.RData")
	# load("./trawlDiversity/pkgBuild/results/msomStatic_norv_shelf_stan_start2016-02-09_00-14-00_r10-9.RData")
	
	# load("./trawlDiversity/pkgBuild/results/msomStatic_norv_ai_stan_start2016-02-11_19-09-15_r2.RData")
	load("./trawlDiversity/pkgBuild/results/msomStatic_norv_goa_stan_start2016-02-11_15-06-49_r3.RData")
	# load("./trawlDiversity/pkgBuild/results/msomStatic_norv_gmex_stan_start2016-02-12_05-34-10_r6.RData")
	
	orig_rm_out <- rm_out
	rm_out <- orig_rm_out[[3]]
	reg <- "GOA"
	
	inputData <- lapply(rm_out, function(x)x$inputData)
	out <- lapply(rm_out, function(x)x$out)
	
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
	X_obs <- lapply(inputData, function(x)x$X[1,,])
	
	# ---- Work on Species Availability Probs ----
	# K <- lapply(inputData, function(x)x$nK[1,])
	# Omega_dist <- lapply(out, function(x)extract(x, pars="Omega")[[1]])
	# lp_unavail_dist <- lapply(Omega_dist, function(x)dbinom(0, size=1, prob=x, log=TRUE))
	# lp_avail_dist1 <- lapply(Omega_dist, function(x)dbinom(1, size=1, prob=x, log=TRUE))
	#
	# dim_match_K <- function(ref, K){
	# 	# a dimension-matching helper function
	# 	# intended for matching K to either theta or psi (as ref)
	# 	# created for the case in which K has a length equal to the second dimension of ref
	#
	# 	# i.e., the second dimension of ref (psi or theta) is usually the 'stratum' index
	# 	# and K is the number of samples in a stratum
	# 	# the other dimensions of ref are things like iteration index (from posterior sampling) and species
	#
	# 	# ref is 2-3 dim
	# 	# K is vector
	# 	stopifnot(is.null(dim(K)))
	# 	stopifnot(any(dim(ref)==length(K)))
	# 	how_k <- which(dim(ref)==length(K))
	# 	K2 <- ref
	# 	if(how_k == 1){ # untested
	# 		K2[] <- K
	# 	} else if(how_k == 2){ # only case currently using, so only one tested
	# 		K1 <- rep(K, each=nrow(ref))
	# 		K2[] <- K1
	# 	} else if(how_k == 3){ # untested
	# 		K1 <- rep(K, each=prod(dim(ref)[2:3]))
	# 		K2[] <- K1
	# 	}
	#
	# 	return(K2)
	#
	# }
	#
	# lp_observed_func <- function(x, K, psi, theta){
	# 	# based on Carpenter's species-site-occupancy.pdf document
	# 	K2 <- dim_match_K(theta, K)
	# 	log(psi) + pbinom(q=x, size=K2, prob=theta, log=TRUE)
	#
	# }
	# lp_unobserved_func <- function(K, psi, theta){
	# 	# based on Carpenter's species-site-occupancy.pdf document
	# 	log(exp(lp_observed_func(0, K, psi, theta) + exp(log(1 - psi))))
	# }
	#
	# lp_avail_dist2 <- mapply(lp_unobserved_func, K, psi_dist, theta_dist, SIMPLIFY=FALSE)
	# match_avail_dist_12 <- function(x,y){
	# 	# helper function to match array dimensions
	# 	# used when x is of length equal to the first dimension of array y
	# 	# simply repeats the values in x across the other y dimensions
	# 	lx <- length(x)
	# 	stopifnot(dim(y)[1]==lx)
	# 	x2 <- y
	# 	x2[] <- x
	#
	# 	return(x2)
	# }
	#
	# lp_avail_dist <- mapply(function(x,y) match_avail_dist_12(x,y) + y, lp_avail_dist1, lp_avail_dist2, SIMPLIFY=FALSE)
	#
	#
	# lse <- function(x,y){
	# 	# like the log_sum_exp function in Stan
	# 	log(exp(x) + exp(y))
	# }
	# get_p_avail <- function(x,y){
	# 	# based on Carpenter's species-site-occupancy.pdf document, in the bottom of the loop for E-N_2
	# 	# get Pr_available given x=lp_available and y=lp_unavailable
	# 	y2 <- match_avail_dist_12(y, x) # note the ordering of the reference and to-be-changed arrays is reversed from previous use of this function
	# 	exp(x - lse(y2, x))
	# }
	# PR_available_dist_pureModel <- mapply(get_p_avail, lp_avail_dist, lp_unavail_dist, SIMPLIFY=FALSE)
	#
	# update_pr_avail_withObs <- function(pure, obs){
	# 	# PR_available_dist_pureModel is the model estimate of each species existing in a given stratum
	# 	# but for strata where we observed the species, we know the probability to be 1
	# 	# This function changes values in PR_available_dist_pureModel
	# 	# to 1 for all indices that correspond to confirmed observations
	#
	# 	d_pure <- dim(pure)
	# 	d_obs <- dim(obs)
	# 	stopifnot(length(d_pure)==3 & length(d_obs)==2)
	# 	stopifnot(d_obs[1]==d_pure[2] & d_obs[2] == d_pure[3])
	# 	obs2 <- aperm(pure, c(2,3,1))
	# 	obs2[] <- obs
	# 	obs2 <- aperm(obs2, c(3,1,2))
	# 	obs_ind <- obs2 > 0
	# 	pure[obs_ind] <- 1
	#
	# 	return(pure)
	# }
	#
	# PR_available_dist <- mapply(update_pr_avail_withObs, PR_available_dist_pureModel, X_obs, SIMPLIFY=FALSE) # the problem with this is that it doesn't really make sense to estimate the species availability on a per-stratum basis. The w parameter is defined per species; only the probability of the species being present varies among strata or years.
	
	# seen <- PR_available_dist[[1]][1,,] == 1
	# PR_available_dist_pureModel[[1]][1,,][seen]
	# PR_available_dist_pureModel[[1]][1,,][!seen]
	
	# ---- Richness in Each Stratum ----
		
	
	
	# ---- Richness in the Region ----
	Omega_mean <- sapply(out, function(x)mean(extract(x, pars="Omega")[[1]]))
	naive_rich <- sapply(inputData, function(x)x$N)
	rich_pureModel <- Omega_mean * naive_rich
	
	reg_pres_dist <- lapply(psi_dist, function(x)apply(x, c(1,3), function(x)(1-prod(1-x))))
	reg_pres <- lapply(reg_pres_dist, function(x)colMeans(x))
	reg_rich <- sapply(reg_pres, sum)
	
	# psi_dist_upObs <- mapply(update_pr_avail_withObs, psi_dist, X_obs, SIMPLIFY=FALSE)
	# reg_pres_dist_upObs <- lapply(psi_dist_upObs, function(x)apply(x, c(1,3), function(x)(1-prod(1-x))))
	# reg_pres_upObs <- lapply(reg_pres_dist_upObs, function(x)colMeans(x))
	# reg_rich_upObs <- sapply(reg_pres_upObs, sum)
	
	
	# ---- Covariates ----
	bt0 <- lapply(inputData, function(x)x$U[1,,"bt"])
	bt_ann <- sapply(bt0, mean)
	
	yr <- 1:length(bt0)
	bt2dt <- function(x,y)data.table(stratum=names(x), bt=x, yr=y)
	bt <- rbindlist(mapply(bt2dt, bt0, yr, SIMPLIFY=FALSE))
	
	strat2lld <- function(x){
		s <- strsplit(x, split=" ")
		lon <- sapply(s, function(x)x[1])
		lat <- sapply(s, function(x)x[2])
		depth_interval <- sapply(s, function(x)x[3])
		data.table(lon=as.numeric(lon), lat=as.numeric(lat), depth_interval=as.numeric(depth_interval))
	} 
	bt[,c("lon","lat","depth_interval"):=strat2lld(stratum)]
	
	zCol <- function(nCols, Z){
		cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(nCols)
		colVec_ind <- cut(Z, breaks=nCols)
		colVec <- cols[colVec_ind]
	}
	bt[,bt_col:=zCol(256, bt)]
	

	
	
	mytrace <- function(x, pars){
		sims <- lapply(x@sim$samples, function(x)x[pars])
		for(i in 1:length(sims)){
			sims[[i]]$chain <- rep(i, length(sims[[i]][[1]]))
		}
		sims <- rbindlist(sims)
		cn <- colnames(sims)
		for(h in 1:(ncol(sims)-1)){
			ylim <- sims[,range(eval(s2c(cn[h]))[[1]])]
			for(i in 1:sims[,lu(chain)]){
				
				if(i == 1){
					sims[chain==i, plot(eval(s2c(cn[h]))[[1]], ylim=ylim, type="l", col=i, xlab="iter",ylab=cn[h])]
				}else{
					sims[chain==i, lines(eval(s2c(cn[h]))[[1]], col=i)]
				}
			}
		}
	}
	
	pars_trace <- c("Omega","alpha_mu[1]", "alpha_mu[2]", "alpha_mu[3]", "alpha_mu[4]", "beta_mu[1]")
	
	
	
	


	

	
	
	# ---- Figures ----
	fig1_name <- paste0("richness_bt_timeSeries_", reg, ".png")
	png(file.path("trawlDiversity/pkgBuild/figures",fig1_name), width=3.5, height=6, units="in", res=200)
	par(mfrow=c(3,1), mar=c(1.75,1.5,0.25,0.25), oma=c(0.1,0.1,0.75,0.1), mgp=c(0.75,0.1,0), tcl=-0.1, ps=8, cex=1)
	plot(naive_rich, type="o", ylab="Naive Region Richness", xlab="Year")
	plot(reg_rich, type="o", xlab="Year", ylab="MSOM Region Richness")
	plot(bt_ann, type="o", xlab="Year", ylab="Annual Mean Bottom Temperature")
	mtext(reg, side=3, line=0, outer=TRUE, font=2)
	dev.off()
	
	fig2_name <- paste0("richness_bt_scatter_", reg, ".png")
	png(file.path("trawlDiversity/pkgBuild/figures",fig2_name), width=3.5, height=5, units="in", res=200)
	par(mfrow=c(2,1), mar=c(1.75,1.5,0.25,0.25), oma=c(0.1,0.1,0.75,0.1), mgp=c(0.75,0.1,0), tcl=-0.1, ps=8, cex=1)
	plot(bt_ann, naive_rich, type="p", ylab="Naive Region Richness", xlab="Annual Mean Bottom Temperature")
	plot(bt_ann, reg_rich, type="p", ylab="MSOM Region Richness", xlab="Annual Mean Bottom Temperature")
	mtext(reg, side=3, line=0, outer=TRUE, font=2)
	dev.off()
	
	fig3_name <- paste0("btempMap_", reg, ".png")
	png(file.path("trawlDiversity/pkgBuild/figures",fig3_name), width=8, height=3, units="in", res=200)
	par(mfrow=auto.mfrow(bt[,lu(yr)]), oma=c(0.1,0.1, 1,0.1), mar=c(1,1,0.1,0.1), mgp=c(0.75,0.1,0), tcl=-0.15, cex=1, ps=8)
	bt[,j={
		plot(lon, lat, type="n")
		map(add=TRUE)
		points(lon, lat, col=bt_col, pch=20)
		mtext(unique(yr), side=3, adj=0.1, line=-0.75, font=2)
	}, by="yr"]
	mtext(paste(reg, "Bottom Temperature"), outer=TRUE, side=3, line=0.25, font=2)
	dev.off()
	
	fig4_name <- paste0("traceplot_", reg, ".png")
	png(file.path("trawlDiversity/pkgBuild/figures",fig4_name), width=12, height=6, units="in", res=200)
	par(mfrow=c(length(pars_trace), length(out)), oma=c(1,1, 1,0.1), mar=c(0.5,0.5,0.1,0.1), mgp=c(0.5,0.1,0), tcl=-0.1, cex=1, ps=6)
	for(h in 1:length(pars_trace)){
		for(i in 1:length(out)){
			mytrace(out[[i]], pars=pars_trace[h])
			if(i == 1){
				mtext(pars_trace[h], side=2, line=0.75)
			}
		}
	}
	dev.off()
	
}
	
	
	
	