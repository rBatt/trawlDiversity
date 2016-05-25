
library(trawlDiversity)


setwd("~/Documents/School&Work/pinskyPost/trawl/")

load("trawlDiversity/pkgBuild/results/processedMsom/p.RData")

# =============
# = Functions =
# =============

# ---- Time Metrics for Colonization/ Extinction ----
# add more columns pertaining to how long spp has been around
streak_length <- function(x, streak_val=1, fill_val=0){
	rl <- rle(x)
	rlv <- rl$values
	streaks <- sapply(rl$lengths[rlv==streak_val & is.finite(rlv)], seq_len)	
	
	finite_x <- is.finite(x)
	x[x!=streak_val | !finite_x] <- fill_val
	x[x==streak_val & finite_x] <- Reduce(c, streaks)
	
	return(x)
}

event_distance <- function(x, positions=seq_along(x), event_value=1, keep_sign=FALSE){
	# distance of each event to nearest non-event
	
	event_index <- which(x == event_value)
	nonEvent_index <- seq_along(x)[-event_index]
	
	if(length(event_index) == 0 | length(nonEvent_index) == 0){
		return(as(rep(NA, length(x)), class(x)))
	}
	
	event_position <- positions[event_index]
	nonEvent_position <- positions[nonEvent_index]
	
	# dists <- abs(outer(event_position, nonEvent_position, "-"))
	# event_dists <- rep(0, length(x))
	# event_dists[event_index] <- apply(dists, 1, min)
	
	dists_s <- outer(nonEvent_position, event_position, "-")
	dists <- abs(dists_s)
	event_dists <- rep(0, length(x))
	signed_dists <- apply(dists_s, 2, function(x)x[which.min(abs(x))])
	if(keep_sign){
		event_dists[event_index] <- signed_dists
	}else{
		event_dists[event_index] <- abs(signed_dists)
	}
	
	return(as(event_dists, class(x)))
}

# ================
# = Data Objects =
# ================
# ---- Generic ----
regs <- sapply(p, function(x)x$processed[,una(reg)])

# ---- Processed ----
processed_dt <- lapply(p, function(x)x$processed)
processed_dt <- rbindlist(processed_dt)[reg!="wcann"]
setkey(processed_dt, reg, year)

# ---- Beta Diversity ----
beta_div_dt <- list()
for(i in 1:length(p)){
	t_reg <- p[[i]]$processed[,una(reg)]
	t_beta_div_dt_obs <- p[[i]]$beta_div_obs
	t_beta_div_dt_msom <- p[[i]]$beta_div
	beta_div_dt[[i]] <- data.table(reg=t_reg, merge(t_beta_div_dt_msom, t_beta_div_dt_obs, by=c("year","method")))
}
beta_div_dt <- rbindlist(beta_div_dt)[reg!="wcann"]
setkey(beta_div_dt, reg, year)

# ---- Detection ----
detect_dt <- list()
for(i in 1:length(p)){
	t_reg <- p[[i]]$processed[,una(reg)]
	t_detect <- p[[i]]$ab[par=="beta", list(year, spp, value_mu, value_sd)]
	detect_dt[[i]] <- data.table(reg=t_reg, t_detect)
}
detect_dt <- rbindlist(detect_dt)[reg!="wcann"]
setkey(detect_dt, reg, year, spp)

# ---- Colonization/ Extinction ----
# colonization/ extinction
ce_dt <- list()
for(i in 1:length(p)){
	t_reg <- p[[i]]$processed[,una(reg)]
	ce_dt[[i]] <- data.table(reg=t_reg, p[[i]]$colonization$col_ext_dt)
}
ce_dt <- rbindlist(ce_dt)[reg!="wcann"]
setkey(ce_dt, reg, year, spp)

# ---- Add Detection to Colonization/ Extinction ----
detect_ce_dt <- merge(detect_dt, ce_dt, all=TRUE)
detect_ce_dt[,present:=as.integer(!is.na(value_mu))]
setnames(detect_ce_dt, c("value_mu", "value_sd"), c("detect_mu", "detect_sd"))

# ---- Add 'Year-Of' Timing for Extinction Flag ----
detect_ce_dt[,now_ext:=(c(0,as.integer(diff(present)==-1))),by=c("reg","spp")]

# ---- Add Time Metrics to Detection + Col/Ext ----
detect_ce_dt[,cumm_yrs_pres:=(cumsum(present)), by=c("reg", "spp")] # cummulative years
detect_ce_dt[, consec_yrs:=streak_length(present), by=c("reg", "spp")] # years in a row
detect_ce_dt[,yrs_since_1st:=cumsum(c(0,1)[(cumsum(present)>=1)+1L]), by=c("reg", "spp")] # years since first
detect_ce_dt[,ext_dist:=event_distance(x=present, positions=year, event_value=1), by=c("reg", "spp")] # years (actual years, not sampling events) away from an absence
detect_ce_dt[,ext_dist_sign:=event_distance(x=present, positions=year, event_value=1, keep_sign=TRUE), by=c("reg", "spp")]
detect_ce_dt[,ext_dist_samp:=event_distance(x=present, event_value=1), by=c("reg", "spp")]


# ====================
# = Richness Figures =
# ====================

# ---- Time Series Naive Richness ----
dev.new()
par(mfrow=c(3,3))
processed_dt[,j={plot(year, naive_rich, type="o", main=una(reg))},by=c("reg")]

# ---- Time Series MSOM Richness ----
dev.new()
par(mfrow=c(3,3))
processed_dt[,j={plot(year, reg_rich, type="o", main=una(reg))},by=c("reg")]

# ---- Scatter Richness Naive vs MSOM ----
dev.new()
par(mfrow=c(3,3))
processed_dt[,j={plot(naive_rich, reg_rich, main=una(reg)); abline(a=0, b=1)},by=c("reg")]


# ==========================
# = Beta Diversity Figures =
# ==========================

# ---- Spatial Variance Time Series Naive ----
dev.new()
par(mfrow=c(3,3))
beta_div_dt[,j={plot(year, beta_div_obs, type="o", main=una(reg))}, by=c("reg")]

# ---- Spatial Variance Time Series MSOM ----
dev.new()
par(mfrow=c(3,3))
beta_div_dt[,j={plot(year, beta_div_mu, type="o", main=una(reg))}, by=c("reg")]

# ---- Scatter Beta Diversity: Obs vs MSOM ----
dev.new()
par(mfrow=c(3,3))
beta_div_dt[,j={plot(beta_div_obs, beta_div_mu, main=una(reg)); abline(a=0, b=1)},by=c("reg")]


# =========================
# = Detectability Figures =
# =========================

# ---- Time Series of Community Average Detectability ----
dev.new()
par(mfrow=c(3,3))
detect_dt[,j={
	comm_detect <- .SD[,list(detctability=plogis(mean(value_mu))),by=c("year")]
	plot(comm_detect, type="o", main=una(reg))
}, by=c("reg")]

# ---- Time Series of Species Detectability ----
dev.new()
par(mfrow=c(3,3))
for(r in 1:length(regs)){
	t_reg <- regs[r]
	t_dt <- detect_dt[reg==t_reg]
	
	xlim <- range(t_dt[,year])
	ylim <- range(plogis(t_dt[,value_mu]))
	
	us <- una(t_dt[,spp])
	t_col <- adjustcolor("black", alpha=0.25)
	for(s in 1:lu(us)){
		xy <- t_dt[spp==us[s], list(year, detectability=plogis(value_mu))]
		if(s==1){
			plot(xy, xlim=xlim, ylim=ylim, type="l", lwd=0.5, main=t_reg, col=t_col)
		}else{
			lines(xy, lwd=0.5, col=t_col)
		}
	}
	
}

# ---- Density Plots of Species-Specific Trends in Detectability ----
dev.new()
par(mfrow=c(3,3))
detect_ce_dt[,j={
	t_slopes <- .SD[,j={as.numeric(timeSlope(cumm_yrs_pres, detect_mu))},by=c("spp")]
	ce_spp <- .SD[col==1 | ext==1, una(spp)]
	
	dens_ce <- t_slopes[spp%in%ce_spp, density(V1, na.rm=TRUE)]
	dens_pres <- t_slopes[!spp%in%ce_spp, density(V1, na.rm=TRUE)]
	
	ylim <- range(c(dens_pres$y, dens_ce$y))
	xlim <- range(c(dens_pres$x, dens_ce$x))
	
	plot(dens_ce[c("x","y")], col="red", xlim=xlim, ylim=ylim, type="l", xlab="Detectability Trend", ylab="Density", main=una(reg))
	lines(dens_pres[c("x","y")], col="blue")
	abline(v=0, lty="dashed")
	
	if(reg[1]=="ai"){
		legend("topright",legend=c("all yrs", "some yrs"), text.col=c("blue","red"), bty="n", title="spp presence", title.col="black")
	}
	
	NULL
	
}, by=c("reg")]


# =====================
# = Proportion Strata =
# =====================
propStrat <- list()
for(i in 1:length(p)){
	t_reg <- p[[i]]$processed[,una(reg)]
	t_rd <- p[[i]]$rd
	
	ssy_tbl <- apply(t_rd[,table(spp, stratum, year)>0], c(1,3), sum)
	n_strat <- t_rd[,colSums(table(stratum,year)>0)]	
	n_strat <- t(matrix(n_strat, nrow=length(n_strat), ncol=nrow(ssy_tbl)))
	prop_strat_tbl <- ssy_tbl/n_strat

	prop_strat_dt_wide <- data.table(spp=rownames(prop_strat_tbl), as.data.table(prop_strat_tbl))
	prop_strat_dt <- data.table:::melt.data.table(prop_strat_dt_wide, id.vars="spp", variable.name="year", value.name="propStrata")
	prop_strat_dt[,year:=as.integer(as.character(year))]
	
	propStrat[[i]] <- data.table(reg=t_reg, prop_strat_dt)
}
propStrat <- rbindlist(propStrat)[reg!="wcann"]
setkey(propStrat, reg, year, spp)


# ---- Time series of Community Average of Proportion Strata ----
dev.new()
par(mfrow=c(3,3))
propStrat[,j={
	comm_propStrat <- .SD[,list(propStrat_mu=mean(propStrata)),keyby=c("year")]
	plot(comm_propStrat, type="o", main=una(reg))
}, by=c("reg")]

# ---- Time Series of Community Average Deviation from 50% Strata Occupied ----
dev.new()
par(mfrow=c(3,3))
propStrat[,j={
	comm_propStrat <- .SD[,list(propStrat_mu_dev50=mean(abs(propStrata - 0.5))),keyby=c("year")]
	plot(comm_propStrat, type="o", main=una(reg))
}, by=c("reg")]

# ---- Time Series of Proportion of Strata Occupied by Each Species ----
dev.new()
par(mfrow=c(3,3))
for(r in 1:length(regs)){
	t_reg <- regs[r]
	t_dt <- propStrat[reg==t_reg]
	
	xlim <- range(t_dt[,year])
	ylim <- range(t_dt[,propStrata])
	
	us <- una(t_dt[,spp])
	t_col <- adjustcolor("black", alpha=0.25)
	for(s in 1:lu(us)){
		xy <- t_dt[spp==us[s], list(year, propStrata)]
		if(s==1){
			plot(xy, xlim=xlim, ylim=ylim, type="l", lwd=0.5, main=t_reg, col=t_col)
		}else{
			lines(xy, lwd=0.5, col=t_col)
		}
	}
}

# ---- Density Plots of Species-Specific Trends in Proportion Strata Occupied ----
dev.new()
par(mfrow=c(3,3))
propStrat[,j={
	t_slopes <- .SD[,j={as.numeric(timeSlope(year, propStrata))},by=c("spp")]
	
	ps_cat <- function(x){
		all_below <- all(x < 0.5 | is.na(x))
		if(all_below){return("below")}
		
		all_above <- all(x > 0.5 | is.na(x))
		if(all_above){return("above")}
		
		return("both")
	}
	spp_cats <- .SD[,list(spp_cat=ps_cat(propStrata)), by="spp"]
	above_spp <- spp_cats[spp_cat=="above",una(spp)]
	below_spp <- spp_cats[spp_cat=="below",una(spp)]
	both_spp <- spp_cats[spp_cat=="both",una(spp)]
	
	d2 <- function(x){
		if(all(is.na(x))){
			return(list(x=NA, y=NA))
		}else{
			return(density(x, na.rm=TRUE))
		}
	}
	
	dens_above <- t_slopes[spp%in%above_spp, d2(V1)]
	dens_below <- t_slopes[spp%in%below_spp, d2(V1)]
	dens_both <- t_slopes[spp%in%both_spp, d2(V1)]
	d_list <- list("above"=dens_above, "below"=dens_below, "both"=dens_both)
	
	ylim <- range(sapply(d_list, function(x)range(x$y)))
	xlim <- range(sapply(d_list, function(x)range(x$x)))
	
	plot(d_list[["above"]][c("x","y")], col="red", xlim=xlim, ylim=ylim, type="l", xlab="Trend in Proportion Strata Occupied", ylab="Density", main=una(reg))
	lines(d_list[["below"]][c("x","y")], col="blue")
	lines(d_list[["both"]][c("x","y")], col="black")
	abline(v=0, lty="dashed")
	
	if(reg[1]=="ebs"){
		legend("topright",legend=c("spp above 50%", "spp that cross", "spp below 50%"), text.col=c("red", "black", "blue"), bty="n")
	}
	
	NULL
	
}, by=c("reg")]

# ---- Scatter Plot of Trend vs Mean Strata Occupied ----
dev.new()
par(mfrow=c(3,3))
propStrat[,j={
	t_slopes_means <- .SD[,j={list(slope=as.numeric(timeSlope(year, propStrata)), mean=mean(propStrata, na.rm=TRUE))},by=c("spp")]
	t_slopes_means[,plot(slope, mean, main=una(reg), ylab="Species Mean % Strata", xlab="Species Slope in % Strata")]
	abline(v=0, h=0.5, col="blue")
	NULL
}, by=c("reg")]


# ====================================================================
# = % Strat, Richness, Beta Diversity, & Detectability Relationships =
# ====================================================================

# ---- richness time series with sparklines for colonizer/ leaver proportion strata ----
# data
spp_master <- merge(detect_ce_dt, propStrat, all=TRUE)
comm_master <- merge(beta_div_dt, processed_dt, all=TRUE)


# define characteristics related to colonization-extinction "stretches"
# for any stretches to exist, the full series will have all(c(0,1)%in%una(present))
spp_master[,has_stretches:=all(c(0,1)%in%una(present)), by=c("reg","spp")]


# present <- detect_ce_dt[reg=="ai" & spp =="Lethotremus muticus", present]
# year <- detect_ce_dt[reg=="ai" & spp =="Lethotremus muticus", year]
# now_ext <- detect_ce_dt[reg=="ai" & spp =="Lethotremus muticus", now_ext]
# col <- detect_ce_dt[reg=="ai" & spp =="Lethotremus muticus", col]
# ext_dist_sign <- detect_ce_dt[reg=="ai" & spp =="Lethotremus muticus", ext_dist_sign]

event_stretches <- function(X){
	# must have columns for present, year, now_ext, col, ext_dist_sign
	# returns columns of stretch_id, hybrid_part
	
	present <- X[, present]
	year <- X[, year]
	now_ext <- X[, now_ext]
	col <- X[, col]
	ext_dist_sign <- X[, ext_dist_sign]

	# For convenience, define rle of present
	rls0 <- rle(present)
	
	# Output components, empty
	# stretch_id = integer vector identifying the type and ID of stretch to which each element belongs
	stretch_id <- rep(0, length(present)) # 0 means no stretch
	hybrid_part <- rep(0, length(present))

	# For any stretches to exist, the full series will have all(c(0,1)%in%una(present)) 
	has_stretches <- all(c(0,1)%in%una(present))
	
	if(!has_stretches){
		return(list(stretch_id=stretch_id, hybrid_part=hybrid_part))
	}

	# If it has stretches, the series will be comprised of stretches taking on one or more of the following forms:
		# extinction-only
		# colonization-only
		# hybrid stretch

	# Define start and stop of e-only and c-only stretches
	has_col <- any(col==1)
	has_now_ext <- any(now_ext==1)
	e_only_start <- (year==min(year) & present == 1) # extinction-only
	e_only_end <- (now_ext==1 & year==min(year[now_ext==1])) # extinction-only
	c_only_start <- (col==1 & year==max(year[col==1])) # colonization-only
	c_only_end <- (year==max(year) & present==1) # colonization-only

	# A e-only or c-only stretch exists if (any(start_logic) & any(end_logic))
	has_eo <- any(e_only_start) & any(e_only_end)
	has_co <- any(c_only_start) & any(c_only_end)

	# The existance and number of hybrid stretches in a series is defined by:
	n_hybrid <- sum(head(rls0$values[-1],-1)==1)

	# ID for extinction-only stretch
	if(has_eo){
		ind <- seq_along(present)
		eo_ind <- ind[e_only_start] : ind[e_only_end]
		stretch_id[eo_ind] <- -1
	}

	# ID for colonization-only stretch
	if(has_co){
		ind <- seq_along(present)
		co_ind <- ind[c_only_start] : ind[c_only_end]
		stretch_id[co_ind] <- -2
	}

	# ID for hybrid stretches
	if(n_hybrid>0){
		rls <- rls0
		rls$values[c(1,length(rls$values))] <- 0
		rls$values <- rls$values * cumsum(rls$values)
		after_hybrid <- c(0, diff(rls$values)) < 0
		rls$lengths <- rls$lengths - as.integer(after_hybrid) # removes the now_ext==1 from non-hybrid classification
		rls$lengths <- rls$lengths + as.integer(rls$values > 0) # adds the now_ext==1 to the hybrid classification
		hybrid_id <- inverse.rle(rls)
		stretch_id[hybrid_id!=0] <- hybrid_id[hybrid_id!=0]
	
		# A single hybrid stretch can be further subdivided into two parts:
			# post-colonization: an initial portion of the stretch that is closer to its start than end
			# pre-extinction: a latter portion of the stretch that is closer to its end than start
		# Categorization as post-c or pre-e corresponds to whether the nearest absence is in the past or future:
			# post-c: ext_dist_sign<0
			# pre-e: ext_dist_sign>=0
		hybrid_part[stretch_id>0 & ext_dist_sign<0] <- 1 # post-colonization (first part)
		hybrid_part[stretch_id>0 & ext_dist_sign>=0] <- 2 # pre-extinction (second part)
	}
	out <- list(stretch_id=stretch_id, hybrid_part=hybrid_part)
	
	
	return(out)
	
}

detect_ce_dt[reg=="ai" & spp =="Lethotremus muticus", data.table(present, as.data.table(event_stretches(.SD)))]
detect_ce_dt[reg=="ai",event_stretches(.SD), by="spp"]

detect_ce_dt[,j={
	
	
	cbind(year, present, col, now_ext, ext_dist_sign, stretch_id, hybrid_part)

	
},by=c("reg","spp")]





# setup figure
dev.new()
par(mfrow=c(3,3))

# loop through each region

	# plot richness time series

	# loop through each richness year
	
		# skip first year
		# identify species that colonize or left
		
		# get time-to-event distances for every species
		# subset event distance vectors to region bounded by 0 distance with current at center
		# split subsetted distance into maximum prior and leading into current
		# split subsetted distance into current leading into maximum after
		
		# define viable leaver species as those having prior maximum distance of at least 2
		# define 
		
		# for each colonizer
			
			
		# for each leaver



		



# ==================
# = exploring fast =
# ==================
detect_ce_dt[ext_dist!=0, plot(ext_dist, detect_mu, main="Detectability vs Time to Absence")]
detect_ce_dt[ext_dist!=0, abline(lm(detect_mu ~ ext_dist), col="blue")]
detect_ce_dt[ext_dist!=0, summary(lm(detect_mu ~ ext_dist))] # significant, 0.04 R2 tho

detect_ce_dt[ext_dist!=0, plot(cumm_yrs_pres, detect_mu, main="Detectability vs Cummulative 'Years' Present")]
detect_ce_dt[ext_dist!=0, abline(lm(detect_mu ~ cumm_yrs_pres), col="blue")]
detect_ce_dt[ext_dist!=0, summary(lm(detect_mu ~ cumm_yrs_pres))] # R2 0.04

# ---- winner ----
detect_ce_dt[ext_dist!=0, plot(consec_yrs, detect_mu, main="Detectability vs Consecutive 'Years' Present")]
detect_ce_dt[ext_dist!=0, abline(lm(detect_mu ~ consec_yrs), col="blue")]
detect_ce_dt[ext_dist!=0, summary(lm(detect_mu ~ consec_yrs))] # R2 0.06

detect_ce_dt[ext_dist!=0, plot(yrs_since_1st, detect_mu, main="Detectability vs 'Years' Since First Present")]
detect_ce_dt[ext_dist!=0, abline(lm(detect_mu ~ yrs_since_1st), col="blue")]
detect_ce_dt[ext_dist!=0, summary(lm(detect_mu ~ yrs_since_1st))] # R2 0.03

detect_ce_dt[ext_dist!=0, predict(lm(detect_mu ~ yrs_since_1st), newdata=data.frame(yrs_since_1st=range(yrs_since_1st)))]

pred_range(mod){
	mr <- lapply(mod$model, range)
}






