
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
for(r in 1:length(regs[regs!="wcann"])){
	t_reg <- regs[regs!="wcann"][r]
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
for(r in 1:length(regs[regs!="wcann"])){
	t_reg <- regs[regs!="wcann"][r]
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
spp_master <- merge(detect_ce_dt, propStrat, all=TRUE) # species-specific data
comm_master <- merge(beta_div_dt, processed_dt, all=TRUE) # community-level data


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
	
	# If the current species does not have both presences and absenences, just reutrn 0's
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

# testing function a bit
# detect_ce_dt[reg=="ai" & spp =="Lethotremus muticus", data.table(present, as.data.table(event_stretches(.SD)))]
# detect_ce_dt[reg=="ai" & spp =="Anoplopoma fimbria", data.table(as.data.table(event_stretches(.SD))),by="spp"]
# detect_ce_dt[reg=="ai",event_stretches(.SD), by="spp"]

stretches <- spp_master[,data.table(year, present, propStrata, col, now_ext, consec_yrs, ext_dist_sign, as.data.table(event_stretches(.SD))), keyby=c("reg","spp")]
# stretches[,event_year:=(year+ext_dist_sign)*c(1,NA)[1L+(ext_dist_sign==0|is.na(ext_dist_sign))],by="reg"]
# stretches[,event_year:=(year+ext_dist_sign)*c(1,NA)[1L+(is.na(ext_dist_sign))],by="reg"]
# stretches[,event_year:=(year+ext_dist_sign),by="reg"]

stretches[stretch_id==-1 | hybrid_part==2, stretch_type:="pre_ext"]
stretches[stretch_id==-2 | hybrid_part==1, stretch_type:="post_col"]

stretches[!is.na(stretch_type),event_year:=c(post_col=min(year), pre_ext=max(year))[stretch_type[1]],by=c("reg","spp","stretch_type", "stretch_id", "hybrid_part")]

blah1 <- stretches[!is.na(stretch_type), list(col=paste(col,collapse=" "),year=paste(year,collapse=" "), event_year=paste(event_year, collapse=" "), pres_seq=paste(present, collapse=" "), max_consec=lu(consec_yrs)) ,by=c("reg","spp","stretch_type", "stretch_id", "hybrid_part")]
blah1[max_consec <6 & max_consec>3]
blah1[col==1]

# Stretch ID Key:
#  -1 is extinction only
#  -2 is colonization only
#  1 is hybrid
# Hybrid Part Key:
#  1 is post-colonization
#  2 is pre-extinction


# setup figure
dev.new(width=6.8, height=7)
# par(mfrow=c(3,3), mar=c(1,1,0.5,0.1), oma=c(0.1,0.1,0.1,0.1), cex=1, ps=6)
par(mar=c(1.5,1.5,0.5,0.1), oma=c(0.1,0.1,0.1,0.1), cex=1, ps=6, mgp=c(0.5, 0.1, 0), tcl=-0.1)

u_regs <- stretches[,una(reg)]

# lay_mat <- matrix(0, ncol=8, nrow=5)
# lay_mat[1,] <- which(u_regs=="ebs") # ebs
# lay_mat[2,] <- which(u_regs=="shelf") # shelf
# lay_mat[3,] <- which(u_regs=="neus") # neus
# lay_mat[4,1:4] <- which(u_regs=="sa")
# lay_mat[4,5:8] <- which(u_regs=="gmex")
# lay_mat[5,1:2] <- which(u_regs=="ai")
# lay_mat[5,3:4] <- which(u_regs=="goa")
# lay_mat[5,5:6] <- which(u_regs=="newf")
# lay_mat[5,7:8] <- which(u_regs=="wctri")

lay_mat <- matrix(0, ncol=12, nrow=4)
lay_mat[1,1:8] <- which(u_regs=="ebs") # ebs
lay_mat[1,9:12] <- which(u_regs=="ai")
lay_mat[2,1:8] <- which(u_regs=="shelf") # shelf
lay_mat[2,9:12] <- which(u_regs=="goa")
lay_mat[3,1:8] <- which(u_regs=="neus") # neus
lay_mat[3,9:12] <- which(u_regs=="wctri")
lay_mat[4,1:4] <- which(u_regs=="sa")
lay_mat[4,5:8] <- which(u_regs=="gmex")
lay_mat[4,9:12] <- which(u_regs=="newf")


layout(lay_mat)

# rr <- 

# loop through each region
u_regs <- stretches[,una(reg)]
for(r in 1:length(u_regs)){
	t_stretches <- stretches[reg==u_regs[r]]
	setorder(t_stretches, year)
	# t_stretches[stretch_id==-1 | hybrid_part==2,table(event_year)] # number of spp-years leading into each year's extinctions
	# t_stretches[stretch_id==-2 | hybrid_part==1,table(event_year)] # number of spp-years following each year's colonizations

		# plot richness time series
		plt_stretch <- t_stretches[,list(richness=sum(present)),keyby="year"]
		c1 <- mean(par("cin"))/par("pin")
		ylim <- plt_stretch[,range(richness)+c(-c1[2],c1[2])*diff(par()$usr[3:4])*c(0.1,0.5)]
		xlim <- plt_stretch[,range(year)+c(-c1[1],c1[1])*diff(par()$usr[1:2])*0.15]
		plot(plt_stretch, type='l', main=u_regs[r], ylim=ylim, xlim=xlim) # use msom?
	
		# t_stretches[propStrata!=0,propStrata:=(propStrata-min(propStrata)), by=c("spp")]
		# t_stretches[propStrata!=0,propStrata:=(propStrata-min(propStrata)), by=c("spp","stretch_id","hybrid_part")]
		# t_stretches[,propStrata:=propStrata/max(propStrata), by="spp"]
		# spark_scale <- t_stretches[,sum(present),by="year"][,diff(range(V1))/10]
		# t_stretches[,propStrata:=propStrata*spark_scale]
	
		# loop through each richness year
		u_ev_yrs <- sort(t_stretches[,una(event_year, na.rm=TRUE)])
		for(y in 1:length(u_ev_yrs)){
			t_ev_yr <- u_ev_yrs[y]
			ty_stretches <- t_stretches[event_year==t_ev_yr]
			t_r <- t_stretches[year==t_ev_yr,sum(present)] # or could use msom
	
			# rescale proportion by species
			# ty_stretches[,propStrata:=(propStrata-min(propStrata[propStrata!=0])+sd(propStrata[propStrata!=0])), by="spp"]
			ty_stretches[,propStrata:=(propStrata-min(propStrata)), by="spp"]
			ty_stretches[,propStrata:=propStrata/max(propStrata), by="spp"]
	
			t_plot <- ty_stretches[,j={
				t_plot <- .SD[,list(year,propStrata)]
				# t_plot <- rbind(data.table(year=event_year[1], propStrata=0), t_plot)
				# t_plot[,propStrata:=propStrata+t_r]
				setorder(t_plot, year)
		
				coloD <- c("blue","red")[(stretch_id==-1 | hybrid_part==2)+1L]
				colo <- adjustcolor(coloD, 0.25)
		
				# lines(t_plot, col=colo)
				# pyr <- t_plot[,year-t_ev_yr] # t_plot[,year - min(year)]
				# pyr <- pyr/20#max(abs(pyr))*2.5
				# pyr <- pyr + t_ev_yr #t_plot[,min(year)]
				# t_plot[,pyr:=pyr]
				# t_plot[,segments(x0=pyr[1],x1=tail(pyr,1), y0=t_r-1E-2, y1=t_r-1E-2, lwd=0.25)]
				# t_plot[,lines(pyr, propStrata, col=colo)]
				
				# t_plot[lines(pyr, propStrata, col=colo)]
				# lines(t_plot, col=colo)
		
				t_plot[,list(year=year, propStrata=propStrata, stretch_type=stretch_type, colo=colo, coloD=(coloD))]
		
			},by=c("spp","stretch_id","hybrid_part")]
	
			# mu_plot <- t_plot[,list(propStrata=mean(propStrata), .N),by=c("year","colo","coloD")]
			# mu_plot[,lines(year, propStrata, col=colo),by=c("colo")]
			
			mu_plot <- t_plot[,list(propStrata=mean(propStrata[is.finite(propStrata)]), .N),by=c("year","colo","coloD", "stretch_type")]
			setorder(mu_plot, year, coloD)
			ust <- mu_plot[,unique(stretch_type)]
			for(st in 1:length(ust)){
				x <- mu_plot[stretch_type==ust[st],year]
				y <- mu_plot[stretch_type==ust[st],propStrata]
				if(length(x)>=3 & sum(is.finite(y))>=3){
					if(sum(is.finite(y))>=4){
						y <- fitted(loess(y~x, data.frame(x,y))) # spline(x, y, n=length(x))$y 
					}
					y_align <- c(pre_ext="right", post_col="left")[ust[st]]
					x_align <- c(pre_ext="right", post_col="left")[ust[st]]
					col <- mu_plot[stretch_type==ust[st],unique(coloD)]
					ax_sides <- list(pre_ext=c(1,2),post_col=c(1,4))[[ust[st]]]
					acol <- mu_plot[stretch_type==ust[st],unique(colo)]
					x_cex <- 0.75 #t_stretches[,1/lu(year)*10]
					y_cex <- x_cex
				
					sparklines(x, y, x_pt=t_ev_yr, y_pt=t_r, x_align=x_align, y_align=y_align, col=col, ax_sides=ax_sides, lwd=0.75, awd=0.5, acol=acol, x_cex=x_cex, y_cex=y_cex)
					
				}
			
			}
			
			# mu_plot[,segments(x0=pyr[1],x1=tail(pyr,1), y0=t_r-1E-2, y1=t_r-1E-2, lwd=0.5)]
		}
}


		



# ==================
# = exploring fast =
# ==================
detect_ce_dt[ext_dist!=0, plot(ext_dist, detect_mu, main="Detectability vs Time to Absence")]
detect_ce_dt[ext_dist!=0, abline(lm(detect_mu ~ ext_dist), col="blue")]
detect_ce_dt[ext_dist!=0, summary(lm(detect_mu ~ ext_dist))] # significant, 0.04 R2 tho

detect_ce_dt[ext_dist!=0, plot(cumm_yrs_pres, detect_mu, main="Detectability vs Cumulative 'Years' Present")]
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







