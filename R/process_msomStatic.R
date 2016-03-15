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
	
	# save_mem=F
	#
	# library("rstan")
	# library("trawlDiversity")
	# library("rbLib")
	# library("R2jags")
	# library("maps")
	#
	# reg_file <- c(
	# 	"msomStatic_norv_1yr_ebs_jags_start2016-03-10_16-44-44_r1.RData", # ebs 6k iter, converged
	# 	# "msomStatic_norv_1yr_ai_jags_start2016-03-11_16-01-45_r2.RData", #ai 6k iter, didn't converge
	# 	"msomStatic_norv_1yr_ai_jags_start2016-03-06_19-10-49_r2.RData", # ai 30k iter, didn't converge
	# 	# "msomStatic_norv_1yr_goa_jags_start2016-03-11_17-55-00_r3.RData", # goa 6k iter, didn't converge
	# 	"msomStatic_norv_1yr_goa_jags_start2016-03-07_06-21-58_r3.RData", # goa 30k iter, didn't converge
	# 	# "msomStatic_norv_1yr_wctri_jags_start2016-03-11_21-14-14_r4.RData", # wctri 6k iter, didn't converge (also too few 0 spp)
	# 	"msomStatic_norv_1yr_wctri_jags_start2016-03-08_10-53-18_r4.RData", # wctri 30k iter, didn't converge (better than most tho); also, needs more than 100 0 spp
	# 	# "msomStatic_norv_1yr_wcann_jags_start2016-03-11_23-18-31_r5.RData", # wcann 6k iter, didn't converge, too few 0 sp (50)
	# 	"msomStatic_norv_1yr_wcann_jags_start2016-03-07_13-25-01_r5.RData", # wcann 30k iter, didn't converge, too few 0 spp
	# 	"msomStatic_norv_1yr_gmex_jags_start2016-03-12_02-53-26_r6.RData", # gmex, 30k iter, didn't converge, too few 0 spp
	# 	# "msomStatic_norv_1yr_sa_jags_start2016-03-12_04-26-53_r7.RData", # sa, 6k iter, almost converged, but too few n0 spp
	# 	"msomStatic_norv_1yr_sa_jags_start2016-03-06_16-15-08_r7.RData", # sa, 30k iter, very nearly converged, but too few n0 spp
	# 	"msomStatic_norv_1yr_neus_jags_start2016-03-12_06-01-10_r8.RData", # neus, 6k iter, mostly converged, but too few n0 spp
	# 	# "msomStatic_norv_1yr_shelf_jags_start2016-03-12_12-00-07_r9.RData", # shelf, 6k iter, very nearly converged
	# 	"msomStatic_norv_1yr_shelf_jags_start2016-03-06_14-03-15_r9.RData", # shelf, 30k iter, better than 6k, but not perfect
	# 	"msomStatic_norv_1yr_newf_jags_start2016-03-12_13-30-24_r10.RData" # newf, 6k iter, really good mixing, mild lack of stationary
	#
	# )[10]
	#
	# load(paste0("trawlDiversity/pkgBuild/results/", reg_file))
	#
	# reg_results_ind <- which(sapply(rm_out, function(x)!is.null(x)))
	# stopifnot(length(reg_results_ind) == 1)
	# reg_out <- rm_out[[reg_results_ind]]
	#
	# if(save_mem){
	# 	rm(list="rm_out")
	# }
	
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
	ab_all <- mapply(get_ab, inputData, out, SIMPLIFY=FALSE)
	for(i in 1:length(ab_all)){
		tyr <- rd_yr[i]
		ab_all[[i]] <- ab_all[[i]][,year:=tyr] 
	}
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
	

	
	
	# # ---- Figures ----
# 	fig1_name <- paste0("richness_bt_timeSeries_", reg, ".png")
# 	# png(file.path("trawlDiversity/pkgBuild/figures",fig1_name), width=3.5, height=6, units="in", res=200)
# 	dev.new(width=3.5, height=6)
# 	par(mfrow=c(3,1), mar=c(1.75,1.5,0.25,0.25), oma=c(0.1,0.1,0.75,0.1), mgp=c(0.75,0.1,0), tcl=-0.1, ps=8, cex=1)
# 	plot(naive_rich, type="o", ylab="Naive Region Richness", xlab="Year")
# 	plot(reg_rich, type="o", xlab="Year", ylab="MSOM Region Richness")
# 	plot(bt_ann, type="o", xlab="Year", ylab="Annual Mean Bottom Temperature")
# 	mtext(reg, side=3, line=0, outer=TRUE, font=2)
# 	# dev.off()
#
# 	fig2_name <- paste0("richness_bt_scatter_", reg, ".png")
# 	# png(file.path("trawlDiversity/pkgBuild/figures",fig2_name), width=3.5, height=5, units="in", res=200)
# 	dev.new(width=3.5, height=5)
# 	par(mfrow=c(2,1), mar=c(1.75,1.5,0.25,0.25), oma=c(0.1,0.1,0.75,0.1), mgp=c(0.75,0.1,0), tcl=-0.1, ps=8, cex=1)
# 	plot(bt_ann, naive_rich, type="p", ylab="Naive Region Richness", xlab="Annual Mean Bottom Temperature")
# 	plot(bt_ann, reg_rich, type="p", ylab="MSOM Region Richness", xlab="Annual Mean Bottom Temperature")
# 	mtext(reg, side=3, line=0, outer=TRUE, font=2)
# 	# dev.off()
#
# 	fig3_name <- paste0("btempMap_", reg, ".png")
# 	# png(file.path("trawlDiversity/pkgBuild/figures",fig3_name), width=8, height=3, units="in", res=200)
# 	dev.new(width=8, height=3)
# 	par(mfrow=auto.mfrow(bt[,lu(year)]), oma=c(0.1,0.1, 1,0.1), mar=c(1,1,0.1,0.1), mgp=c(0.75,0.1,0), tcl=-0.15, cex=1, ps=8)
# 	bt[,j={
# 		plot(lon, lat, type="n")
# 		map(add=TRUE)
# 		points(lon, lat, col=bt_col, pch=20)
# 		mtext(unique(year), side=3, adj=0.1, line=-0.75, font=2)
# 	}, by="yr"]
# 	mtext(paste(reg, "Bottom Temperature"), outer=TRUE, side=3, line=0.25, font=2)
# 	# dev.off()
#
# 	fig4_name <- paste0("traceplot_", reg, ".png")
# 	# png(file.path("trawlDiversity/pkgBuild/figures",fig4_name), width=12, height=6, units="in", res=200)
# 	dev.new(width=12, height=6)
# 	par(mfrow=c(length(pars_trace), length(out)), oma=c(1,1, 1,0.1), mar=c(0.5,0.5,0.1,0.1), mgp=c(0.25,0.1,0), tcl=-0.1, cex=1, ps=6)
# 	for(h in 1:length(pars_trace)){
# 		for(i in 1:length(out)){
# 			mytrace(out[[i]], pars=pars_trace[h], lang=lang, xaxt='n')
# 			if(i == 1){
# 				mtext(pars_trace[h], side=2, line=0.75)
# 			}
# 		}
# 	}
# 	# dev.off()
#
#
# 	# ---- look for correlation in posterior distribution of parameteres ----
# 	# i.e., do parameter values covary
# 	dev.new()
# 	pairs(get_iters(out[[1]], pars_trace, "JAGS"))
#
#
#
#
# 	dev.new()
# 	processed[,plot(unobs_rich[-length(unobs_rich)], n_col[-1], xlab="Unobserved species present last year", ylab="Species colonizing this year")]
# 	abline(a=0, b=1)
#
#
# 	# ---- Number of Colonizations per Stratum ----
# 	dev.new(width=6, height=3)
# 	par(mfrow=c(1,2), mar=c(1.5,1.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=8, cex=1)
# 	colonization$n_spp_col_weighted_tot[,plot_space(lon, lat, n_spp_col_weighted, pch=19)]
# 	map(add=TRUE, fill=TRUE, col="white")
#
# 	colonization$n_spp_col_weighted_tot[,plot_space(lon, lat, n_spp_col_weighted, TRUE, pch=19)]
# 	map(add=TRUE, fill=TRUE, col="white")
#
#
# 	# ---- Location of Colonizations ----
# 	dev.new()
# 	par(mfrow=auto.mfrow(ab[,lu(spp)]), mar=c(1,1,1,0.1), ps=6, mgp=c(0.6,0.1,0), tcl=-0.1, cex=1)
# 	yc <- ab[,zCol(lu(year), una(year))]
# 	names(yc) <- ab[,una(year)]
# 	ab[par=="beta",j={
# 		xyd <- .SD[,density(value)[c("x","y")],by="year"]
# 		xlim <- xyd[,range(x)]
# 		ylim <- xyd[,range(y)]
#
# 		# plot(1,1, type='n', xlim=xlim, ylim=ylim)
# 		# xyd[,lines(x,y,col=yc[as.character(unique(year))]),by="year"]
# 		ts <- una(spp)
# 		tcom <- rd[spp==ts,una(common)]
# 		si <- sppImg(ts, common=tcom)
# 		if(!is.null(si)){
# 			par(new=T)
# 		}
# 		plot(year, value, xlim=range(rd_yr), col=adjustcolor('gray', 0.01), cex=0.5, pch=21, bg=adjustcolor('white',0.01))
# 		mu <- .SD[,list(mu=mean(value)),by="year"]
# 		mu[,lines(year, mu, lwd=2, col='gray')]
# 		mu[,lines(year, mu, lwd=1, col='white')]
# 		if(is.null(si)){
# 			mtext(paste(ts, tcom, sep="\n"), side=3)
# 		}
#
# 	},by=c("spp")]

	return(list(rd=rd, colonization=colonization, bt=bt, param_iters=param_iters, processed=processed, ab=ab))
	
}
	
	
	
	