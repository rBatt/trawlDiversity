
library(trawlDiversity)


setwd("~/Documents/School&Work/pinskyPost/trawl/")

load("trawlDiversity/pkgBuild/results/processedMsom/p.RData")


# ====================
# = Richness Figures =
# ====================
processed_dt <- list()
for(i in 1:length(p)){
	processed_dt[[i]] <- p[[i]]$processed
}
processed_dt <- rbindlist(processed_dt)[reg!="wcann"]
setkey(processed_dt, reg, year)

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
beta_div_dt <- list()
for(i in 1:length(p)){
	t_reg <- p[[i]]$processed[,una(reg)]
	t_beta_div_dt_obs <- p[[i]]$beta_div_obs
	t_beta_div_dt_msom <- p[[i]]$beta_div
	beta_div_dt[[i]] <- data.table(reg=t_reg, merge(t_beta_div_dt_msom, t_beta_div_dt_obs, by=c("year","method")))
}
beta_div_dt <- rbindlist(beta_div_dt)[reg!="wcann"]
setkey(beta_div_dt, reg, year)

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
# detection
detect_dt <- list()
for(i in 1:length(p)){
	t_reg <- p[[i]]$processed[,una(reg)]
	t_detect <- p[[i]]$ab[par=="beta", list(year, spp, value_mu, value_sd)]
	detect_dt[[i]] <- data.table(reg=t_reg, t_detect)
}
detect_dt <- rbindlist(detect_dt)[reg!="wcann"]
setkey(detect_dt, reg, year, spp)

# colonization/ extinction
ce_dt <- list()
for(i in 1:length(p)){
	t_reg <- p[[i]]$processed[,una(reg)]
	ce_dt[[i]] <- data.table(reg=t_reg, p[[i]]$colonization$col_ext_dt)
}
ce_dt <- rbindlist(ce_dt)[reg!="wcann"]
setkey(ce_dt, reg, year, spp)

# detection + colonization/extinction
detect_ce_dt <- merge(detect_dt, ce_dt, all=TRUE)
detect_ce_dt[,present:=as.integer(!is.na(value_mu))]
setnames(detect_ce_dt, c("value_mu", "value_sd"), c("detect_mu", "detect_sd"))

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
detect_ce_dt[,cumm_yrs_pres:=(cumsum(present)), by=c("reg", "spp")] # cummulative years
detect_ce_dt[, consec_yrs:=streak_length(present), by=c("reg", "spp")] # years in a row
detect_ce_dt[,yrs_since_1st:=cumsum(c(0,1)[(cumsum(present)>=1)+1L]), by=c("reg", "spp")] # years since first

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
regs <- detect_dt[,una(reg)]
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

#
# detect_dt[,j={
# 	xlim <- range(year)
# 	ylim <- range(plogis(value_mu))
#
# 	us <- una(spp)
# 	for(s in 1:lu(spp)){
# 		xy <- .SD[spp==us[s], list(year, detectability=plogis(value_mu))]
# 		if(s==1){
# 			plot(xy, xlim=xlim, ylim=ylim, type="l", lwd=0.5, main=una(reg))
# 		}else{
# 			lines(xy, lwd=0.5)
# 		}
# 	}
# }, by=c("reg")]

# ---- Density Plots of Species-Specific Trends in Detectability ----
dev.new()
par(mfrow=c(3,3))
detect_ce_dt[,j={
	t_slopes <- .SD[,j={as.numeric(timeSlope(cumm_yrs_pres, detect_mu))},by=c("spp")]
	ce_spp <- .SD[col==1 | ext==1, una(spp)]
	
	# if(length(ce_spp)==0 | all(t_slopes[,spp%in%ce_spp])){
		# print(t_slopes, nrow=Inf)
	# }
	
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






