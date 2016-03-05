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
	
	library("rstan")
	library("trawlDiversity")
	library("rbLib")
	library("R2jags")
	library("maps")
	
	
	# load("trawlDiversity/pkgBuild/results/msomStatic_norv_1yr_shelf_jags_start2016-03-02_23-14-33_r9.RData")
	# reg <- "shelf"
	# lang <- "JAGS"
	
	load("trawlDiversity/pkgBuild/results/msomStatic_norv_1yr_goa_jags_start2016-03-03_04-46-31_r3.RData")
	reg <- "goa"
	lang <- "JAGS"
	
	reg_results_ind <- which(sapply(rm_out, function(x)!is.null(x)))
	stopifnot(length(reg_results_ind) == 1)
	reg_out <- rm_out[[reg_results_ind]]
	
	inputData <- lapply(reg_out, function(x)x$inputData)
	out <- lapply(reg_out, function(x)x$out)
	
	
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
		
		X_obs <- lapply(inputData, function(x)x$X[1,,])
	}else{
		X_obs <- lapply(inputData, function(x)x$X)
	}


	
	
	
	# ---- Richness in the Region ----
	Omega_iter <- lapply(out, get_iters, pars="Omega", lang="JAGS")
	Omega_mean <- sapply(Omega_iter, function(x)x[,mean(Omega)])
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
		z_spp <- function(x){
			gsub("Z\\[[0-9]+\\,([0-9]+)\\]$", "\\1", x)
		}
		z_j <- function(x){
			gsub("Z\\[([0-9]+)\\,[0-9]+\\]$", "\\1", x)
		}
		
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

		
	}

# plot(reg_rich, type="o")
	
	# psi_dist_upObs <- mapply(update_pr_avail_withObs, psi_dist, X_obs, SIMPLIFY=FALSE)
	# reg_pres_dist_upObs <- lapply(psi_dist_upObs, function(x)apply(x, c(1,3), function(x)(1-prod(1-x))))
	# reg_pres_upObs <- lapply(reg_pres_dist_upObs, function(x)colMeans(x))
	# reg_rich_upObs <- sapply(reg_pres_upObs, sum)
	
	
	# ---- Covariates ----
	# bt0 <- lapply(inputData, function(x)x$U[1,,"bt"])
	bt0 <- lapply(inputData, function(x)x$U[,"bt"])
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
	

	
	
	mytrace <- function(x, pars, lang, ...){
		
		sims <- get_iters(x, pars, lang=lang)
		
		cn <- colnames(sims)
		for(h in 1:(ncol(sims)-1)){
			ylim <- sims[,range(eval(s2c(cn[h]))[[1]])]
			n_chains <- sims[,lu(chain)]
			cols <- adjustcolor(1:n_chains, 0.5)
			for(i in 1:n_chains){
				
				if(i == 1){
					sims[chain==i, plot(eval(s2c(cn[h]))[[1]], ylim=ylim, type="l", col=cols[i], xlab="",ylab=cn[h], ...)]
				}else{
					sims[chain==i, lines(eval(s2c(cn[h]))[[1]], col=cols[i])]
				}
			}
		}
	}
	
	
	if(lang == "Stan"){
		pars_trace <- c("Omega","alpha_mu[1]", "alpha_mu[2]", "alpha_mu[3]", "alpha_mu[4]", "beta_mu[1]")
	}else{
		pars_trace <- c("Omega","alpha_mu[1]", "alpha_mu[2]", "alpha_mu[3]", "alpha_mu[4]", "beta_mu")
	}
	
	
	
	


	

	
	
	# ---- Figures ----
	fig1_name <- paste0("richness_bt_timeSeries_", reg, ".png")
	# png(file.path("trawlDiversity/pkgBuild/figures",fig1_name), width=3.5, height=6, units="in", res=200)
	dev.new(width=3.5, height=6)
	par(mfrow=c(3,1), mar=c(1.75,1.5,0.25,0.25), oma=c(0.1,0.1,0.75,0.1), mgp=c(0.75,0.1,0), tcl=-0.1, ps=8, cex=1)
	plot(naive_rich, type="o", ylab="Naive Region Richness", xlab="Year")
	plot(reg_rich, type="o", xlab="Year", ylab="MSOM Region Richness")
	plot(bt_ann, type="o", xlab="Year", ylab="Annual Mean Bottom Temperature")
	mtext(reg, side=3, line=0, outer=TRUE, font=2)
	# dev.off()
	
	fig2_name <- paste0("richness_bt_scatter_", reg, ".png")
	# png(file.path("trawlDiversity/pkgBuild/figures",fig2_name), width=3.5, height=5, units="in", res=200)
	dev.new(width=3.5, height=5)
	par(mfrow=c(2,1), mar=c(1.75,1.5,0.25,0.25), oma=c(0.1,0.1,0.75,0.1), mgp=c(0.75,0.1,0), tcl=-0.1, ps=8, cex=1)
	plot(bt_ann, naive_rich, type="p", ylab="Naive Region Richness", xlab="Annual Mean Bottom Temperature")
	plot(bt_ann, reg_rich, type="p", ylab="MSOM Region Richness", xlab="Annual Mean Bottom Temperature")
	mtext(reg, side=3, line=0, outer=TRUE, font=2)
	# dev.off()
	
	fig3_name <- paste0("btempMap_", reg, ".png")
	# png(file.path("trawlDiversity/pkgBuild/figures",fig3_name), width=8, height=3, units="in", res=200)
	dev.new(width=8, height=3)
	par(mfrow=auto.mfrow(bt[,lu(yr)]), oma=c(0.1,0.1, 1,0.1), mar=c(1,1,0.1,0.1), mgp=c(0.75,0.1,0), tcl=-0.15, cex=1, ps=8)
	bt[,j={
		plot(lon, lat, type="n")
		map(add=TRUE)
		points(lon, lat, col=bt_col, pch=20)
		mtext(unique(yr), side=3, adj=0.1, line=-0.75, font=2)
	}, by="yr"]
	mtext(paste(reg, "Bottom Temperature"), outer=TRUE, side=3, line=0.25, font=2)
	# dev.off()
	
	fig4_name <- paste0("traceplot_", reg, ".png")
	# png(file.path("trawlDiversity/pkgBuild/figures",fig4_name), width=12, height=6, units="in", res=200)
	dev.new(width=12, height=6)
	par(mfrow=c(length(pars_trace), length(out)), oma=c(1,1, 1,0.1), mar=c(0.5,0.5,0.1,0.1), mgp=c(0.25,0.1,0), tcl=-0.1, cex=1, ps=6)
	for(h in 1:length(pars_trace)){
		for(i in 1:length(out)){
			mytrace(out[[i]], pars=pars_trace[h], lang=lang, xaxt='n')
			if(i == 1){
				mtext(pars_trace[h], side=2, line=0.75)
			}
		}
	}
	# dev.off()
	
}
	
	
	
	