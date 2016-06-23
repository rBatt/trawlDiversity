library("rbLib")
library("vegan")
library("maps")
library("spatstat")
library("fields")
library("raster")
library("spdep")
library("trawlDiversity")

setwd("~/Documents/School&Work/pinskyPost/trawl/")

load("trawlDiversity/pkgBuild/results/processedMsom/p.RData")


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
	prop_strat_dt <- merge(prop_strat_dt, p[[i]]$processed[,list(reg,year,bt_ann)], by=c("year"), all=TRUE)
	setcolorder(prop_strat_dt, c("reg","year","spp","propStrata","bt_ann"))
	
	propStrat[[i]] <- prop_strat_dt #data.table(reg=t_reg, prop_strat_dt)
}
propStrat <- rbindlist(propStrat)[reg!="wcann"]
setkey(propStrat, reg, year, spp)


# ---- Proportion of Tows ----
propTows <- list()
for(i in 1:length(p)){
	t_reg <- p[[i]]$processed[,una(reg)]
	t_rd <- p[[i]]$rd
	
	ssky_tbl <- apply(t_rd[,table(spp, stratum, K, year)>0], c(1,2,4), sum)
	n_tows_sy <- t_rd[,apply(table(stratum,year,K)>0,1:2,sum)]
	n_tows_y <- colSums(n_tows_sy)
	spp_tows <- ssky_tbl
	spp_tows[] <- rep(n_tows_sy, each=nrow(spp_tows))
	prop_tow_ssy <- ssky_tbl/spp_tows
	
	prop_tow_tot <- apply(prop_tow_ssy, c(1,3), mean, na.rm=TRUE)
	prop_tow_occ <- apply(prop_tow_ssy, c(1,3), function(x)mean(x[x>0], na.rm=TRUE)) # function(x)sum(x[x>0])/max(1, sum(x>0))
	prop_tow_occ_multi <- apply(prop_tow_ssy*(spp_tows>1), c(1,3), function(x)mean(x[x>0], na.rm=TRUE))
	
	tbl2dt <- function(tbl, name){
		tbl_dt_wide <- data.table(spp=rownames(tbl), as.data.table(tbl))
		tbl_dt <- data.table:::melt.data.table(tbl_dt_wide, id.vars="spp", variable.name="year", value.name=name)
		tbl_dt[,year:=as.integer(as.character(year))]
		tbl_dt[,reg:=t_reg]
		return(tbl_dt)
	}
	
	propTow_tot <- tbl2dt(prop_tow_tot, "propTow_tot")
	propTow_occ <- tbl2dt(prop_tow_occ, "propTow_occ")
	propTow_occMulti <- tbl2dt(prop_tow_occ_multi, "propTow_occMulti")
	
	propTow <- merge(propTow_tot,propTow_occ, by=c("reg","year","spp"), all=TRUE)
	propTow <- merge(propTow,propTow_occMulti, by=c("reg","year","spp"), all=TRUE)
	setcolorder(propTow, c("reg","year","spp","propTow_tot","propTow_occ","propTow_occMulti"))
	propTow[,c("propTow_tot_avg","propTow_occ_avg","propTow_occMulti_avg"):=lapply(list(propTow_tot,propTow_occ,propTow_occMulti), meanna, na.rm=TRUE), by=c("reg","spp")]
	
	propTows[[i]] <- propTow #data.table(reg=t_reg, prop_strat_dt)
}
propTows <- rbindlist(propTows)[reg!="wcann"]
setkey(propTows, reg, year, spp)


# ---- Response Metrics ----
if(file.exists("~/Documents/School&Work/pinskyPost/trawl/trawlDiversity/pkgBuild/results/resp_metrics.RData")){
	load("~/Documents/School&Work/pinskyPost/trawl/trawlDiversity/pkgBuild/results/resp_metrics.RData")
}else{
	u_regs <- sapply(p, function(x)x$processed[,una(reg)])
	resp_metrics <- list()
	for(r in 1:length(p)){
		t_reg <- u_regs[r]
		if(t_reg == "wcann"){
			next
		}
		t_ua <- p[[r]]$alpha_unscale
		X_bt <- p[[r]]$bt
	
		t_resp_metrics <- t_ua[,j={
			ty <- year[1]
			t_resp <- get_response(.SD, X=X_bt) # msom predictions for probability present based on depth and temp
			t_opt_resps <- t_resp[resp>=quantile(resp, 0.95, na.rm=TRUE), ] # given environment, top 5% of samples to contain the spp
			t_opt_depth <- t_opt_resps[,median(depth)] # median of depths in 5% most likely to contain the spp
			b0_at_od <- .SD[,a1+a4*t_opt_depth+a5*t_opt_depth^2] # intercept given depth is ~optimal depth
			t_opt <- median(.SD[,psi.opt(a2,a3)], na.rm=TRUE) # optimal bottom temp
			t_tol <- median(.SD[,psi.tol(a3)], na.rm=TRUE) # bottom temp tolerance
			t_max <- median(.SD[,psi.max(a1,a2,a3)], na.rm=TRUE) # max probability present
			data.table(bt_opt=t_opt, bt_tol=t_tol, prob_max=t_max, depth_base_opt=t_opt_depth)
		},by=c("spp","year")]

		# The idea behind btemp_atOptDepth is to find the annual bottom temperature for the strata that are at the depths typically occupied by each species
		# First, compare all depths ever seen in each strata to the "optimal depth" values (optimal depth = median depth in top 5% of predicted occurrences for each year)
		# Second, for each stratum, find its depth value (depths can vary among years depending on sampling location w/in stratum) closest to any optimum (there is an optimum for each year)
		# Third, rank the strata by their smallest deviation from the optimal depth value (rank is according to its best depth in any year)
		# Fourth, subset/extract the top 10% best strata according to their rank in step 3
		# Fifth, in each year, calculate the mean bottom temperature in the strata identified in step 4
		# For each species in each year, we now have the average temperature in the strata that, according to their depths, are most likely to be occupied by the species
		# A similar idea would be to take the annual average temperature across strata every occupied by a species; however, this might not be restrictive enough to really focus in on the "best" locations. But if the above approach is too complex for a paper, maybe something along those lines would work.
		btemp_atOptDepth <- t_resp_metrics[,j={
			strata_devs <- X_bt[, list(bt, optDepth_dev=abs(depth-depth_base_opt)),by=c("stratum","year")] # deviations from opt depth
			n_strat_samp <- X_bt[,ceiling(lu(stratum)*0.1)] # 10% = n strata
			ODS <- strata_devs[,list(optDepth_dev_best=(optDepth_dev[which.min((optDepth_dev))])),by="stratum"] # smallest deviation; ODS = optimal depth strata
			setorder(ODS, optDepth_dev_best) # order strata by their smallest deviation
			best_strata <- ODS[1:pmin(n_strat_samp, nrow(ODS)), stratum] # select top 10% (or all, if somehow less available [shouldn't happen])
			btemp_ODS <- strata_devs[stratum%in%best_strata,list(btemp_ODS=mean(bt, na.rm=TRUE), nODS=lu(stratum)),by="year"] # get mean bottom temperature in best strata for each year
			btemp_ODS
		},by="spp"]
	
		# merge btODS w/ response metrics, all year-spp combos
		setkey(btemp_atOptDepth, spp,year)
		btemp_atOptDepth <- merge(btemp_atOptDepth, btemp_atOptDepth[,CJ(spp=unique(spp),year=unique(year))], by=c("spp","year"), all=TRUE)
		btemp_atOptDepth[is.na(nODS), nODS:=0]
		resp_metrics[[r]] <- merge(t_resp_metrics, btemp_atOptDepth, by=c("spp","year"), all=TRUE)
		resp_metrics[[r]][,reg:=t_reg]
		print(r); flush.console()
	}
	resp_metrics <- rbindlist(resp_metrics)
	save(resp_metrics, file="~/Documents/School&Work/pinskyPost/trawl/trawlDiversity/pkgBuild/results/resp_metrics.RData")
}


# ================================
# = Merge into Data Sets to Save =
# ================================
# ---- Species Master Data Set ----
spp_master <- merge(detect_ce_dt, propStrat, all=TRUE) # species-specific data
spp_master[,has_stretches:=all(c(0,1)%in%una(present)), by=c("reg","spp")]
stretches <- spp_master[(has_stretches),data.table(year, as.data.table(event_stretches(.SD))), keyby=c("reg","spp")]
spp_master <- merge(spp_master, stretches, all=TRUE, by=c("reg","spp","year"))

spp_master[stretch_id==-1 | hybrid_part==2, stretch_type:="pre_ext"]
spp_master[stretch_id==-2 | hybrid_part==1, stretch_type:="post_col"]
spp_master[!is.na(stretch_type),event_year:=c(post_col=min(year), pre_ext=max(year))[stretch_type[1]],by=c("reg","spp","stretch_type", "stretch_id", "hybrid_part")]
spp_master[!is.na(stretch_type), stretch_length:=(lu(year)-(stretch_type=="pre_ext")), by=c("reg","spp","stretch_type", "event_year")]

spp_master <- merge(spp_master, resp_metrics, by=c("reg","spp","year"), all=TRUE)

spp_master[,c("bt_opt_avg","bt_tol_avg","detect_mu_avg"):=list(bt_opt_avg=mean(bt_opt, na.rm=TRUE), bt_tol_avg=mean(bt_tol, na.rm=TRUE), detect_mu_avg=mean(detect_mu, na.rm=TRUE)),by=c('reg','spp')]

define_ce_categ <- function(X){
	ext <- X[,ext]
	col <- X[,col]
	if(all(col==0) & all(ext==0)){
		return("neither")
	}
	if(all(col==0) & any(ext==1)){
		return("leaver")
	}
	if(any(col==1) & all(ext==0)){
		return("colonizer")
	}
	if(any(col==1) & any(ext==1)){
		return("both")
	}
}
spp_master[,ce_categ:=define_ce_categ(.SD),by=c("reg","spp")]

spp_master <- merge(spp_master, propTows, by=c("reg","spp","year"), all=TRUE)

# ---- Community Master Data Set ----
comm_metrics <- spp_master[present==1,j={
	lapply(.SD, mean, na.rm=TRUE)
},by=c("reg","year"), .SDcols=c("bt_opt_avg","bt_tol_avg","detect_mu_avg","detect_mu")]

comm_tows <- spp_master[present==1,j={
	lapply(.SD, mean, na.rm=TRUE)
},by=c("reg","year"), .SDcols=c("propTow_tot","propTow_occ","propTow_occMulti","propTow_tot_avg", "propTow_occ_avg", "propTow_occMulti_avg")]

comm_master <- merge(beta_div_dt, processed_dt, all=TRUE) # community-level data
comm_master <- merge(comm_master, comm_metrics, by=c("reg","year"), all=TRUE)
comm_master <- merge(comm_master, spp_master[present==1,list(propStrata_avg=mean(propStrata)),by=c("reg","year")], by=c("reg","year"), all=TRUE)
comm_master <- merge(comm_master, comm_tows, by=c("reg","year"), all=TRUE)

# ---- Map Data ----
mapDat <- make_mapDat(p)

# ---- Function to help outline a region ----
regOutline <- function(X){
	dev.new()
	outlines <- list()

	rs <- X[,una(reg)]
	nr <- length(rs)
	for(r in 1:nr){
		td <- X[reg==rs[r]]
		rast <- raster(xmn=td[,min(lon)-1], xmx=td[,max(lon)+1], ymn=td[,min(lat)-1], ymx=td[,max(lat)+1], res=c(0.5,0.5))
		td[,plot(rasterize(cbind(x=lon,y=lat), rast))]
		td[,points(lon, lat, cex=5)]
		map(add=TRUE, fill=F)
		to <- locator(type='o', pch=20, col='blue')
		outlines[[r]] <- data.table(reg=rs[r], lonP=to$x, latP=to$y)
	}
	
	return(rbindlist(outlines))
}
# outlines <- regOutline(mapDat)
# save(outlines, file="~/Documents/School&Work/pinskyPost/trawl/trawlDiversity/data/outlines.RData")


# ---- make window polygon for ppp ----
make_owin <- function(X, outlines){
	rs <- X[,una(reg)]
	nr <- length(rs)
	X_owin <- list()
	for(r in 1:nr){
		td <- X[reg==rs[r]]
		o <- outlines[reg==rs[r], list(x=lonP, y=latP)]
		o_r <- outlines[reg==rs[r], list(x=rev(lonP), y=rev(latP))]
	
		xr <- td[,range(c(lon,o[,x]))]
		yr <- td[,range(c(lat,o[,y]))]
		ow_exp1 <- bquote(owin(xr, yr, poly=o)) # use when I correctly traced outline counterclockwise
		ow_exp2 <- bquote(owin(xr, yr, poly=o_r)) # reversed (use if traced outline in wrong direction)
	
		X_owin[[rs[r]]] <- tryCatch(td[, eval(ow_exp1)], error=function(cond)td[, eval(ow_exp2)])
	}
	return(X_owin)
}
mapOwin <- make_owin(mapDat, outlines)
lapply(mapOwin, plot)


# ---- do autocorrelation ----
# dev.new(width=7, height=3)
png("~/Desktop/FigureS4_nb_moranI.png", width=7, height=3, units='in', res=200)
map_layout <- trawl_layout()
par(mar=c(0.25,0.25,0.25,0.25), mgp=c(0.25,0.075,0), tcl=-0.1, ps=8, cex=1, oma=c(0.1,0.1,0.1,0.1))
layout(map_layout)
pretty_reg <- c("ebs"="E. Bering Sea", "ai"="Aleutian Islands", "goa"="Gulf of Alaska", "wctri"="West\nCoast\nUS", "gmex"="Gulf of Mexico", "sa"="Southeast US", "neus"="Northeast US", "shelf"="Scotian Shelf", "newf"="Newfoundland")

X <- copy(mapDat)
rs <- X[,una(reg)]
nr <- length(rs)
localAC <- list()
for(r in 1:nr){
	td <- X[reg==rs[r]]
	localAC[[rs[r]]]$td <- td
	
	locs_ll <- as.matrix(td[,list(lon,lat)])
	locs <- simplify2array(ll2km(td[,lon], td[,lat]))[,2:1]
	
	nn1 <- knn2nb(knearneigh(locs, k=1))
	max2NDist <- max(unlist(nbdists(nn1, locs)))
	localAC[[rs[r]]]$max2NDist <- max2NDist
	localAC[[rs[r]]]$nb <- dnearneigh(locs, d1=0, d2=max2NDist) # graph2nb(gabrielneigh(locs))
	
	ml_out <- localmoran(td[,n_spp_col_weighted], nb2listw(localAC[[rs[r]]]$nb), p.adjust.method="BH")
	localAC[[rs[r]]]$I <- ml_out
	
	# mod <- lm(n_spp_col_weighted~1, data=td)
	# ml_ex <- localmoran.exact(mod, nb=dnearneigh(locs, d1=0, d2=125))
	# ml_sad <- localmoran.sad(mod, nb=dnearneigh(locs, d1=0, d2=125))
	# ex_p <- sapply(ml_ex, function(x)x$p.value)
	# sad_p <- sapply(ml_sad, function(x)x$p.value)
	# ml_p <- ml_out[,5] #sapply(ml_out, function(x)x$p.value)
	# plot(data.frame(ex_p, sad_p, ml_p))
	# abline(a=0, b=1)
	
	plot(mapOwin[[r]], coords=locs_ll, add=F, main="")
	box()
	if(rs[r]=='wctri'){
		mtext(pretty_reg[rs[r]], side=3, line=-3)
	}else{
		mtext(pretty_reg[rs[r]], side=3, line=-1)
	}
	plot(localAC[[rs[r]]]$nb, locs_ll, add=T, col=adjustcolor('black', 0.5), cex=0.8, lwd=0.5)
	sig_lac <- ml_out[,5]<0.05
	if(any(sig_lac)){
		zl <- range(ml_out[sig_lac,1], na.rm=TRUE)
		t_col <- zCol(256,ml_out[sig_lac,1])
		map_col <- zCol(6, 1:6)
		print(max(ml_out[,1]))
		points(x=locs_ll[sig_lac,1], y=locs_ll[sig_lac,2], bg=t_col, pch=21, cex=1.1)
		switch(rs[r],
			ebs = mapLegend(x=0.05, y=0.25, h=0.375, w=0.025, zlim=zl, cols=map_col, lab.cex=1),
			ai = mapLegend(x=0.985, y=0.3, w=0.02, h=0.5, zlim=zl, cols=map_col, lab.cex=1),
			goa = mapLegend(x=0.985, y=0.15, w=0.02,  zlim=zl, cols=map_col, lab.cex=1),
			wctri = mapLegend(x=0.1, y=0.125, w=0.07, zlim=zl, cols=map_col, lab.cex=1),
			gmex = mapLegend(x=0.95, y=0.2, h=0.375, zlim=zl, cols=map_col, lab.cex=1),
			sa = mapLegend(x=0.95, y=0.15, zlim=zl, cols=map_col, lab.cex=1),
			neus = mapLegend(x=0.95, y=0.15, zlim=zl, cols=map_col, lab.cex=1),
			shelf = mapLegend(x=0.95, y=0.15, zlim=zl, cols=map_col, lab.cex=1),
			newf = mapLegend(x=0.05, y=0.15, h=0.25, zlim=zl, cols=map_col, lab.cex=1)
		)
	}else{
		# plot(x=locs[,1], y=locs[,2], col='blue')
	}
	# switch(rs[r],
	# 	ebs = legend("topright", legend="A", bty='n', text.font=2, inset=c(-0.02,-0.15), cex=1.25),
	# 	ai = legend("topleft", legend="C", bty='n', text.font=2, inset=c(-0.065,-0.45), cex=1.25, xpd=T),
	# 	goa = legend("topleft", legend="B", bty='n', text.font=2, inset=c(-0.065,-0.06), cex=1.25),
	# 	wctri = legend("top", legend="E", bty='n', text.font=2, inset=c(0,0.15), cex=1.25),
	# 	gmex = legend("topleft", legend="G", bty='n', text.font=2, inset=c(-0.175,-0.12), cex=1.25),
	# 	sa = legend("topleft", legend="I", bty='n', text.font=2, inset=c(-0.15,-0.075), cex=1.25),
	# 	neus = legend("topleft", legend="H", bty='n', text.font=2, inset=c(-0.125,-0.05), cex=1.25),
	# 	shelf = legend("topleft", legend="D", bty='n', text.font=2, inset=c(-0.1,-0.125), cex=1.25),
	# 	newf = legend("topright", legend="F", bty='n', text.font=2, inset=c(-0.01,-0.05), cex=1.25)
	# )
	
	# # plots
	# dev.new(width=5, height=5)
	# par(mfrow=c(2,2), mar=c(2,2,1,0.5), ps=8, mgp=c(1,0.1,0), tcl=-0.1, cex=1)
	#
	# # moran's i vs neighborhood size
	# plot(dists, t_mt, type='o', main="Sensitivity to Neighbor Distance", xlab="Max Distance to Neighbor (km)", ylab="Moran's I")
	#
	# # correlogram
	# if(min(card(nb))>2){
	#   correlo <- sp.correlogram(nb, td[,n_spp_col_weighted], order = min(card(nb)), method = "I", style = "B")
	# 	plot(correlo)
	# }else{
	# 	plot(1, type='n', xaxt='n', yaxt='n', ylab='', xlab='')
	# }
	#
	# # value at focal vs at lag 1 (moran plot)
	# moran.plot(td[,n_spp_col_weighted], localAC[[rs[r]]]$nb)
	#
	# # spatial sca
	# if(any(ml_out[,5]<0.05)){
	# 	plot_space(x=locs[,1], y=locs[,2], z=ml_out[,1]*(ml_out[,5]<0.05), scatter=TRUE)
	# }else{
	# 	plot(x=locs[,1], y=locs[,2], col='blue')
	# }
}
dev.off()




# =================================
# = Save Data Ojbects for Package =
# =================================
save(mapDat, file="trawlDiversity/data/mapDat.RData")
save(spp_master, file="trawlDiversity/data/spp_master.RData")
save(comm_master, file="trawlDiversity/data/comm_master.RData")

