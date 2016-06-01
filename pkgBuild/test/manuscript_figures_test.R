
library("trawlDiversity")
library("rbLib")
library("vegan")
library("maps")
library("spatstat")
library("fields")


setwd("~/Documents/School&Work/pinskyPost/trawl/")

load("trawlDiversity/pkgBuild/results/processedMsom/p.RData")


# =============
# = Functions =
# =============
# ---- make mapDat ----
make_mapDat <- function(p){
	mapDat <- list()
	for(r in 1:length(p)){
		pt1 <- trawlAgg(
			p[[r]]$rd,
			bioFun=meanna, envFun=meanna,
			envCols=c("btemp","depth","stemp","lon","lat"),
			bio_lvl="spp",time_lvl="year",space_lvl="stratum",
			metaCols=c("reg"),meta.action="unique1"	
		)
		pt2 <- pt1[,j={
			avgRich <- .SD[,lu(spp),by="time_lvl"][,meanna(V1)]
			sdRich <- .SD[,lu(spp),by="time_lvl"][,sd(V1, na.rm=TRUE)]
			avgBtemp <- .SD[,meanna(btemp),by="time_lvl"][,meanna(V1)]
			sdBtemp <- .SD[,meanna(btemp),by="time_lvl"][,sd(V1, na.rm=TRUE)]
			data.table(avgRich=avgRich, sdRich=sdRich, avgBtemp=avgBtemp, sdBtemp=sdBtemp)
		},by=c("reg","stratum")]
		to_merge <- c(p[[r]]$colonization[c("n_spp_col_weighted_tot","n_spp_ext_weighted_tot","n_cep")])
		mapDat[[r]] <- merge(to_merge[["n_spp_col_weighted_tot"]], to_merge[["n_spp_ext_weighted_tot"]], by=c("stratum","lon","lat","depth"),all=TRUE)
		mapDat[[r]] <- merge(mapDat[[r]], pt2, by="stratum",all=TRUE)
	}
	mapDat <- rbindlist(mapDat)[reg!="wcann"]
	# mapDim <- mapDat[,list(r_lon=diff(range(lon)),r_lat=diff(range(lat))),by="reg"]
	# mapDim[,c("lon_scale","lat_scale"):=list(r_lon/min(r_lon), r_lat/min(r_lat))]
	# mapDim[,ll_ratio:=r_lon/r_lat]
	
	return(mapDat)
}


# ---- Packed Map Layout ----
trawl_layout <- function(){
	lay_grid <- matrix(1:84, nrow=6, ncol=14)
	squares <- list(
		ebs = c(1,26),
		ai = c(6,42),
		goa = c(3,41),
		wctri = c(45,48),
		gmex = c(53,66),
		sa = c(70,84), # should really be c(70,77), ... just made it bigger to fill in gap and b/c it's a square
		neus = c(67,81),
		shelf = c(31,44),
		newf = c(49,64)
	)
	squares_ind <- lapply(squares, function(x)arrayInd(which(lay_grid%in%x),.dim=dim(lay_grid)))
	map_layout <- array(NA, dim(lay_grid))
	for(i in 1:length(squares_ind)){
		map_layout[squares_ind[[i]]] <- i
	}
	map_layout[is.na(map_layout)] <- 0
	
	return(map_layout)
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

# ---- Proportion of Strata ----
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


# ---- Community Master Data Set ----
comm_master <- merge(beta_div_dt, processed_dt, all=TRUE) # community-level data

# ---- Species Master Data Set ----
spp_master <- merge(detect_ce_dt, propStrat, all=TRUE) # species-specific data
spp_master[,has_stretches:=all(c(0,1)%in%una(present)), by=c("reg","spp")]
stretches <- spp_master[(has_stretches),data.table(year, as.data.table(event_stretches(.SD))), keyby=c("reg","spp")]
spp_master <- merge(spp_master, stretches, all=TRUE, by=c("reg","spp","year"))

spp_master[stretch_id==-1 | hybrid_part==2, stretch_type:="pre_ext"]
spp_master[stretch_id==-2 | hybrid_part==1, stretch_type:="post_col"]
spp_master[!is.na(stretch_type),event_year:=c(post_col=min(year), pre_ext=max(year))[stretch_type[1]],by=c("reg","spp","stretch_type", "stretch_id", "hybrid_part")]
spp_master[!is.na(stretch_type), stretch_length:=(lu(year)-(stretch_type=="pre_ext")), by=c("reg","spp","stretch_type", "event_year")]


# ---- Map Data ----
mapDat <- make_mapDat(p)


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
# setup figure
dev.new(width=6.8, height=7)
# pdf("~/Desktop/richness_occupancySpark.pdf", width=6.8, height=8)
par(mar=c(1.5,1.5,0.5,0.1), oma=c(0.1,0.1,0.1,0.1), cex=1, ps=8, mgp=c(0.65, 0.1, 0), tcl=-0.1)

u_regs <- spp_master[,una(reg)]
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


# loop through each region
u_regs <- spp_master[,una(reg)]
for(r in 1:length(u_regs)){
	t_spp_master <- spp_master[reg==u_regs[r]]
	setorder(t_spp_master, year)
		# plot richness time series
		plt_stretch <- t_spp_master[,list(richness=sum(present)),keyby="year"]
		c1 <- mean(par("cin"))/par("pin")
		ylim <- plt_stretch[,range(richness)+c(-c1[2],c1[2])*diff(par()$usr[3:4])*c(0.5,1)]
		xlim <- plt_stretch[,range(year)+c(-c1[1],c1[1])*diff(par()$usr[1:2])*0.15]
		plot(plt_stretch, type='l', main=u_regs[r], ylim=ylim, xlim=xlim, lwd=2) # use msom?
	
		# loop through each richness year
		u_ev_yrs <- sort(t_spp_master[,una(event_year, na.rm=TRUE)])
		for(y in 1:length(u_ev_yrs)){
			t_ev_yr <- u_ev_yrs[y]
			ty_spp_master <- t_spp_master[event_year==t_ev_yr]
			t_r <- t_spp_master[year==t_ev_yr,sum(present)] # or could use msom
	
			# rescale proportion by species
			ty_spp_master[,propStrata:=(propStrata-min(propStrata)), by="spp"]
			ty_spp_master[,propStrata:=propStrata/max(propStrata), by="spp"]
	
			t_plot <- ty_spp_master[,j={
				t_plot <- .SD[,list(year,propStrata)]
				setorder(t_plot, year)
				coloD <- c("blue","red")[(stretch_id==-1 | hybrid_part==2)+1L]
				colo <- adjustcolor(coloD, 0.25)
				t_plot[,list(year=year, propStrata=propStrata, stretch_type=stretch_type, colo=colo, coloD=(coloD))]
			},by=c("spp","stretch_id","hybrid_part")]
			
			mu_plot <- t_plot[,list(propStrata=mean(propStrata[is.finite(propStrata)]), .N),by=c("year","colo","coloD", "stretch_type")]
			setorder(mu_plot, year, coloD)
			ust <- mu_plot[,unique(stretch_type)]
			for(st in 1:length(ust)){
				x <- mu_plot[stretch_type==ust[st],year]
				y <- mu_plot[stretch_type==ust[st],propStrata]
				if(length(x)>=3 & sum(is.finite(y))>=3){
					if(sum(is.finite(y))>=3){
						y <- fitted(lm(y~x, data.frame(x,y))) #fitted(loess(y~x, data.frame(x,y))) # spline(x, y, n=length(x))$y 
					}
					y_align <- c(pre_ext="right", post_col="left")[ust[st]]
					x_align <- c(pre_ext="right", post_col="left")[ust[st]]
					col <- mu_plot[stretch_type==ust[st],unique(coloD)]
					ax_sides <- list(pre_ext=c(1,2),post_col=c(1,4))[[ust[st]]]
					acol <- mu_plot[stretch_type==ust[st],unique(colo)]
					x_cex <- 0.75 #t_spp_master[,1/lu(year)*10]
					y_cex <- x_cex
				
					rbLib::sparklines(x, y, x_pt=t_ev_yr, y_pt=t_r, x_align=x_align, y_align=y_align, col=col, ax_sides=ax_sides, lwd=0.75, awd=0.5, acol=acol, x_cex=x_cex, y_cex=y_cex)	
				}
			}
		}
}
# dev.off()



# ====================================================================
# = Regressions of propStrata x Detection x Richness x BetaDiversity =
# ====================================================================
# ---- propStrata vs Detectability ----
dev.new()
par(mfrow=c(3,3))
spp_master[,j={plot(propStrata,plogis(detect_mu), main=reg[1]);abline(a=0,b=1,col='blue')},by="reg"]


# ---- propStrata vs Year from Extinction/ Colonization ----
# dev.new()
# par(mfrow=c(3,3))
# spp_master[!is.na(stretch_type), plot(ext_dist_sign, propStrata, main=reg[1]), by="reg"]

dev.new()
par(mfrow=c(3,3))
spp_master[!is.na(stretch_type), j={
	tp <- .SD[,list(ext_dist_sign=ext_dist_sign,propStrata=propStrata-min(propStrata)),by="spp"]
	tp <- tp[,list(ext_dist_sign=ext_dist_sign,propStrata=propStrata/max(propStrata)),by="spp"]
	tp[stretch_length>3,plot(ext_dist_sign, propStrata, main=reg[1])]
}, by="reg"]


# ---- Beta Diversity vs Richness ----
dev.new()
par(mfrow=c(3,3))
comm_master[,plot(beta_div_obs,reg_rich,main=reg[1]),by="reg"]

dev.new()
par(mfrow=c(3,3))
comm_master[,plot(diff(beta_div_obs),diff(reg_rich),main=reg[1]),by="reg"]

dev.new()
par(mfrow=c(3,3))
comm_master[,j={ccf(beta_div_obs,reg_rich,main=reg[1]);NULL},by="reg"]

dev.new()
par(mfrow=c(3,3))
comm_master[,j={ccf(diff(beta_div_obs),diff(reg_rich),main=reg[1]);NULL},by="reg"]


# =======================================================
# = Colonization, Extinction, Richness Maps and Scatter =
# =======================================================
# ---- test map plotting ----
# dev.new(height=3, width=7)
# par(mar=c(1.5,1.5,0.5,0.5), mgp=c(0.75,0.1,0), tcl=-0.1,ps=8, cex=1)
# layout(map_layout)
# u_regs <- mapDat[,unique(reg)]
# for(r in 1:lu(u_regs)){
# 	# mapDat[reg==u_regs[r], plot(lon,lat, col=as.factor(u_regs)[r], bty='l', xlim=map_xlims[,r], ylim=map_ylims[,r])]
# 	mapDat[reg==u_regs[r], plot(lon,lat, col=as.factor(u_regs)[r], bty='l')]
# 	map(add=TRUE)
# }

# ---- plots ----
map_names <- c("Average Richness", "Richness Variability", "Colonization Rate", "Colonizations per Species", "Extinction Rate", "Extinctions per Species")
map_expr <- list(
	bquote(avgRich),
	bquote(sdRich),
	bquote(n_spp_col_weighted),
	bquote(n_spp_col_weighted/avgRich),
	bquote(n_spp_ext_weighted),
	bquote(n_spp_ext_weighted/avgRich)
)

map_layout <- trawl_layout()

for(mp in 1:length(map_names)){
	dev.new(height=3, width=7)
	par(mar=c(1.5,1.5,0.5,0.5), mgp=c(0.75,0.1,0), tcl=-0.1,ps=8, cex=1, oma=c(0.5,0.5,1,0.1))
	layout(map_layout)
	u_regs <- mapDat[,unique(reg)]
	for(r in 1:lu(u_regs)){
		# mapDat[reg==u_regs[r], plot(lon,lat, col=as.factor(u_regs)[r], bty='l', xlim=map_xlims[,r], ylim=map_ylims[,r])]
		mapDat[reg==u_regs[r], plot_space(lon,lat, eval(map_expr[[mp]]), bty='l')]
		map(add=TRUE, fill=TRUE, col="white")
	}
	mtext(map_names[mp], side=3, outer=TRUE, font=2, line=0)
}


# ---- scatter ----
dev.new()
par(mfrow=c(3,3))
mapDat[,plot(n_spp_col_weighted, avgRich, main=reg[1]), by="reg"]

dev.new()
par(mfrow=c(3,3))
mapDat[,plot(n_spp_ext_weighted, avgRich, main=reg[1]), by="reg"]



# ===================================================
# = MSOM Parameter Response Metrics (opt, tol, max) =
# ===================================================
psi.opt <- function(b1,b2){-b1/(2*b2)}
psi.tol <- function(b2){1/sqrt(-2*b2)}
psi.max <- function(b0,b1,b2){1/(1+exp((b1^2)/(4*b2)-b0))}

r = 1
t_ua <- p[[r]]$alpha_unscale

t_ua[,list(psi.opt),by=c("spp","year")]








