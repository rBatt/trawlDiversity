
library("trawlDiversity")
library("rbLib")
library("vegan")
library("maps")
library("spatstat")
library("fields")


# setwd("~/Documents/School&Work/pinskyPost/trawl/")
# load("trawlDiversity/pkgBuild/results/processedMsom/p.RData")

# ---- Generic ----
regs <- comm_master[,una(reg)] #sapply(p, function(x)x$processed[,una(reg)])
pretty_reg <- c("ebs"="E. Bering Sea", "ai"="Aleutian Islands", "goa"="Gulf of Alaska", "wctri"="West Coast US", "gmex"="Gulf of Mexico", "sa"="Southeast US", "neus"="Northeast US", "shelf"="Scotian Shelf", "newf"="Newfoundland")
pretty_col <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','navy','#a65628','salmon','#999999')
names(pretty_col) <- names(pretty_reg)

# ====================
# = Richness Figures =
# ====================

# # ---- Time Series Naive Richness ----
# dev.new()
# par(mfrow=c(3,3))
# comm_master[,j={plot(year, naive_rich, type="o", main=una(reg))},by=c("reg")]

# ---- Time Series MSOM Richness ----
dev.new()
par(mfrow=c(3,3))
comm_master[,j={plot(year, reg_rich, type="o", main=una(reg))},by=c("reg")]

# ---- Scatter Richness Naive vs MSOM ----
dev.new()
par(mfrow=c(3,3))
comm_master[,j={plot(naive_rich, reg_rich, main=una(reg)); abline(a=0, b=1)},by=c("reg")]


# ==========================
# = Beta Diversity Figures =
# ==========================
#
# # ---- Spatial Variance Time Series Naive ----
# dev.new()
# par(mfrow=c(3,3))
# comm_master[,j={plot(year, beta_div_obs, type="o", main=una(reg))}, by=c("reg")]
#
# # ---- Spatial Variance Time Series MSOM ----
# dev.new()
# par(mfrow=c(3,3))
# comm_master[,j={plot(year, beta_div_mu, type="o", main=una(reg))}, by=c("reg")]
#
# # ---- Scatter Beta Diversity: Obs vs MSOM ----
# dev.new()
# par(mfrow=c(3,3))
# comm_master[,j={plot(beta_div_obs, beta_div_mu, main=una(reg)); abline(a=0, b=1)},by=c("reg")]


# =========================
# = Detectability Figures =
# =========================

# # ---- Time Series of Community Average Detectability ----
dev.new()
par(mfrow=c(3,3))
comm_master[,j={
	plot(year, detect_mu, type="o", main=una(reg))
}, by=c("reg")]

dev.new()
par(mfrow=c(3,3))
comm_master[,j={
	plot(year, detect_mu_avg, type="o", main=una(reg))
}, by=c("reg")]

# # ---- Time Series of Species Detectability ----
# dev.new()
# par(mfrow=c(3,3))
# for(r in 1:length(regs[regs!="wcann"])){
# 	t_reg <- regs[regs!="wcann"][r]
# 	t_dt <- spp_master[reg==t_reg]
#
# 	xlim <- range(t_dt[,year])
# 	ylim <- range(plogis(t_dt[,detect_mu]), na.rm=TRUE)
#
# 	us <- una(t_dt[,spp])
# 	t_col <- adjustcolor("black", alpha=0.25)
# 	for(s in 1:lu(us)){
# 		xy <- t_dt[spp==us[s], list(year, detectability=plogis(detect_mu))]
# 		if(s==1){
# 			plot(xy, xlim=xlim, ylim=ylim, type="l", lwd=0.5, main=t_reg, col=t_col)
# 		}else{
# 			lines(xy, lwd=0.5, col=t_col)
# 		}
# 	}
#
# }

# ---- Density Plots of Species-Specific Trends in Detectability ----
dev.new()
par(mfrow=c(3,3))
spp_master[,j={
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
# # ---- Time series of Community Average of Proportion Strata ----
# dev.new()
# par(mfrow=c(3,3))
# comm_master[,plot(year, propStrata_avg, type='o', main=reg[1]),by='reg']

# # ---- Time Series of Community Average Deviation from 50% Strata Occupied ----
# dev.new()
# par(mfrow=c(3,3))
# spp_master[,j={
# 	comm_propStrat <- .SD[,list(propStrat_mu_dev50=mean(abs(propStrata - 0.5))),keyby=c("year")]
# 	plot(comm_propStrat, type="o", main=una(reg))
# }, by=c("reg")]
#
# # ---- Time Series of Proportion of Strata Occupied by Each Species ----
# dev.new()
# par(mfrow=c(3,3))
# for(r in 1:length(regs[regs!="wcann"])){
# 	t_reg <- regs[regs!="wcann"][r]
# 	t_dt <- spp_master[reg==t_reg]
#
# 	xlim <- range(t_dt[,year])
# 	ylim <- range(t_dt[,propStrata], na.rm=TRUE)
#
# 	us <- una(t_dt[,spp])
# 	t_col <- adjustcolor("black", alpha=0.25)
# 	for(s in 1:lu(us)){
# 		xy <- t_dt[spp==us[s], list(year, propStrata)]
# 		if(s==1){
# 			plot(xy, xlim=xlim, ylim=ylim, type="l", lwd=0.5, main=t_reg, col=t_col)
# 		}else{
# 			lines(xy, lwd=0.5, col=t_col)
# 		}
# 	}
# }
#
# # ---- Density Plots of Species-Specific Trends in Proportion Strata Occupied ----
# dev.new()
# par(mfrow=c(3,3))
# spp_master[,j={
# 	t_slopes <- .SD[,j={as.numeric(timeSlope(year, propStrata))},by=c("spp")]
#
# 	ps_cat <- function(x){
# 		all_below <- all(x < 0.5 | is.na(x))
# 		if(all_below){return("below")}
#
# 		all_above <- all(x > 0.5 | is.na(x))
# 		if(all_above){return("above")}
#
# 		return("both")
# 	}
# 	spp_cats <- .SD[,list(spp_cat=ps_cat(propStrata)), by="spp"]
# 	above_spp <- spp_cats[spp_cat=="above",una(spp)]
# 	below_spp <- spp_cats[spp_cat=="below",una(spp)]
# 	both_spp <- spp_cats[spp_cat=="both",una(spp)]
#
# 	d2 <- function(x){
# 		if(all(is.na(x))){
# 			return(list(x=NA, y=NA))
# 		}else{
# 			return(density(x, na.rm=TRUE))
# 		}
# 	}
#
# 	dens_above <- t_slopes[spp%in%above_spp, d2(V1)]
# 	dens_below <- t_slopes[spp%in%below_spp, d2(V1)]
# 	dens_both <- t_slopes[spp%in%both_spp, d2(V1)]
# 	d_list <- list("above"=dens_above, "below"=dens_below, "both"=dens_both)
#
# 	ylim <- range(sapply(d_list, function(x)range(x$y)))
# 	xlim <- range(sapply(d_list, function(x)range(x$x)))
#
# 	plot(d_list[["above"]][c("x","y")], col="red", xlim=xlim, ylim=ylim, type="l", xlab="Trend in Proportion Strata Occupied", ylab="Density", main=una(reg))
# 	lines(d_list[["below"]][c("x","y")], col="blue")
# 	lines(d_list[["both"]][c("x","y")], col="black")
# 	abline(v=0, lty="dashed")
#
# 	if(reg[1]=="ebs"){
# 		legend("topright",legend=c("spp above 50%", "spp that cross", "spp below 50%"), text.col=c("red", "black", "blue"), bty="n")
# 	}
#
# 	NULL
#
# }, by=c("reg")]
#
# # ---- Scatter Plot of Trend vs Mean Strata Occupied ----
# dev.new()
# par(mfrow=c(3,3))
# spp_master[,j={
# 	t_slopes_means <- .SD[,j={list(slope=as.numeric(timeSlope(year, propStrata)), mean=mean(propStrata, na.rm=TRUE))},by=c("spp")]
# 	t_slopes_means[,plot(slope, mean, main=una(reg), ylab="Species Mean % Strata", xlab="Species Slope in % Strata")]
# 	abline(v=0, h=0.5, col="blue")
# 	NULL
# }, by=c("reg")]


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


# # ---- Beta Diversity vs Richness ----
# dev.new()
# par(mfrow=c(3,3))
# comm_master[,plot(beta_div_obs,reg_rich,main=reg[1]),by="reg"]
#
# dev.new()
# par(mfrow=c(3,3))
# comm_master[,plot(diff(beta_div_obs),diff(reg_rich),main=reg[1]),by="reg"]
#
# dev.new()
# par(mfrow=c(3,3))
# comm_master[,j={ccf(beta_div_obs,reg_rich,main=reg[1]);NULL},by="reg"]
#
# dev.new()
# par(mfrow=c(3,3))
# comm_master[,j={ccf(diff(beta_div_obs),diff(reg_rich),main=reg[1]);NULL},by="reg"]


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
# 	maps::map(add=TRUE)
# }

# ---- plots ----
# map_names <- c("Average Richness", "Richness Variability", "Colonization Rate", "Colonizations per Species", "Extinction Rate", "Extinctions per Species")
map_names <- c("Average Richness", "Colonizations per Species")
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


# # ---- scatter ----
# dev.new()
# par(mfrow=c(3,3))
# mapDat[,plot(n_spp_col_weighted, avgRich, main=reg[1]), by="reg"]
#
# dev.new()
# par(mfrow=c(3,3))
# mapDat[,plot(n_spp_ext_weighted, avgRich, main=reg[1]), by="reg"]



# ===================================================
# = MSOM Parameter Response Metrics (opt, tol, max) =
# ===================================================
# # these graphs were originally only meant for 1 region at a time, so not too useful
# dev.new(width=10, height=10)
# par(mfrow=auto.mfrow(spp_master[,lu(spp)]), mar=c(0.5,0.5,0.5,0.5), cex=1, ps=6, mgp=c(0.5,0.1,0), tcl=-0.1, oma=c(0.1,0.1,1,0.1))
# opt_lims <- spp_master[,range(bt_opt, na.rm=TRUE)]
# spp_master[,plot(density(bt_opt, from=opt_lims[1], to=opt_lims[2], na.rm=TRUE), main=spp[1]),by="spp"]
# mtext(paste(t_reg, "Optimal Bottom Temperature"), side=3, outer=TRUE, line=0)
#
# dev.new(width=10, height=10)
# par(mfrow=auto.mfrow(spp_master[,lu(spp)]), mar=c(0.5,0.5,0.5,0.5), cex=1, ps=6, mgp=c(0.5,0.1,0), tcl=-0.1, oma=c(0.1,0.1,1,0.1))
# tol_lims <- spp_master[,range(bt_tol, na.rm=TRUE)]
# spp_master[,plot(density(bt_tol, from=tol_lims[1], to=tol_lims[2], na.rm=TRUE), main=spp[1]),by="spp"]
# mtext(paste(t_reg, "Tolerance of Bottom Temperature"), side=3, outer=TRUE, line=0)

# # ---- time series of thermal optimum (using that year's model only) ----
# dev.new()
# par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
# spp_master[, plot(.SD[,list(avg_thermal_opt=mean(bt_opt, na.rm=TRUE)),by="year"],main=reg[1], type='o'), by=c("reg")]
#
# # ---- time series of thermal optimum (long-term average per species, only species present in that year) ----
# dev.new()
# par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
# comm_master[,plot(year, bt_opt_avg, type="o", main=reg[1]),by="reg"]
#
# # ---- difference between observed temperature and optimum, for each year, compare present & absent species ----
# dev.new()
# par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
# spp_master[,j={
# 	# temp_dev <- .SD[, list(opt_temp_deviation=mean(bt_opt_avg - btemp_ODS, na.rm=TRUE)),keyby=c('year','present')]
# 	# temp_dev[,plot(year, abs(opt_temp_deviation), col=c("red","blue")[(present+1)], main=reg[1])]
# 	# temp_dev[present==1,lines(year, abs(opt_temp_deviation), col=c("red","blue")[(present+1)], main=reg[1])]
# 	# temp_dev[present==0,lines(year, abs(opt_temp_deviation), col=c("red","blue")[(present+1)], main=reg[1])]
# 	# .SD[,plot(year, abs(bt_opt_avg - btemp_ODS), col=factor(present))]
#
# 	# oops, but kinda interesting (below)
# 	# instead of species optimum, i plotted region temperature vs temperature in the best depths (best for each spp)
# 	# temp_dev <- .SD[, list(opt_temp_deviation=mean(bt_ann - btemp_ODS, na.rm=TRUE)),keyby=c('year','present')]
# 	# temp_dev[,plot(year, abs(opt_temp_deviation), col=c("red","blue")[(present+1)], main=reg[1])]
# 	# temp_dev[present==1,lines(year, abs(opt_temp_deviation), col=c("red","blue")[(present+1)], main=reg[1])]
# 	# temp_dev[present==0,lines(year, abs(opt_temp_deviation), col=c("red","blue")[(present+1)], main=reg[1])]
# 	# if(reg[1]=="ai"){legend("topleft",legend=c("present","absent"), text.col=c("blue","red"), inset=c(-0.15,-0.05))}
#
# 	temp_dev <- .SD[, list(opt_temp_deviation=mean(bt_opt_avg - bt_ann, na.rm=TRUE)),keyby=c('year','present')]
# 	temp_dev[,plot(year, abs(opt_temp_deviation), col=c("red","blue")[(present+1)], main=reg[1])]
# 	temp_dev[present==1,lines(year, abs(opt_temp_deviation), col=c("red","blue")[(present+1)], main=reg[1])]
# 	temp_dev[present==0,lines(year, abs(opt_temp_deviation), col=c("red","blue")[(present+1)], main=reg[1])]
# 	if(reg[1]=="ai"){legend("topleft",legend=c("present","absent"), text.col=c("blue","red"), inset=c(-0.15,-0.05))}
# },by=c('reg')]
#
# # ---- scatter plot of annual btemp and opt temp ----
# dev.new()
# par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
# spp_master[present==1, j={plot(.SD[,list(bt_ann=mean(bt_ann,na.rm=T), mean_bt_opt_avg=mean(bt_opt_avg,na.rm=T)),by='year'][,list(bt_ann,mean_bt_opt_avg)],main=reg[1]);abline(a=0,b=1)},by=c('reg')]
#
# # ---- boxplots of optimal temperature for C/E categories ----
# dev.new()
# par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
# spp_master[,j={
# 	dt <- .SD[,list(bt_opt_avg=mean(bt_opt_avg),ce_categ=una(ce_categ)),by=c("spp")]
# 	dtbp <- boxplot(bt_opt_avg~ce_categ, main=reg[1], data=dt, ylab="bt_opt_avg")
# 	NULL
# },by=c("reg")]
#
# # ---- boxplots of temperature tolerance for C/E categories ----
# dev.new()
# par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
# spp_master[,j={
# 	dt <- .SD[,list(bt_tol_avg=mean(bt_tol_avg),ce_categ=una(ce_categ)),by=c("spp")]
# 	dtbp <- boxplot(bt_tol_avg~ce_categ, main=reg[1], data=dt, ylab="bt_tol_avg")
# 	NULL
# },by=c("reg")]
#
# # ---- scatter plot of prevalence (% strata) vs temp tolerance ----
# # actually shows a relationship
# dev.new()
# par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
# spp_master[,j={
# 	dt <- .SD[,list(bt_tol_avg=mean(bt_tol_avg),mean_prevalence=mean(propStrata)),by=c("spp")]
# 	dt[,plot(bt_tol_avg, mean_prevalence, main=reg[1])]
# 	NULL
# },by=c("reg")]
#
# # ---- scatter plot of prevalence vs deviation between thermal optimum and ODS temp ----
# dev.new()
# par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
# spp_master[,j={
# 	dt <- .SD[,list(bt_opt_ods_dev=abs(bt_opt_avg-btemp_ODS),prevalence=(propStrata)),by=c("spp")]
# 	dt[,plot(bt_opt_ods_dev, prevalence, main=reg[1])]
# 	NULL
# },by=c("reg")]


# ================================
# = Lucky Detection and Richness =
# ================================
dev.new()
par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
comm_master[,j={
	plot(plogis(detect_mu_avg), reg_rich, main=reg[1], col=zCol(256,year))
},by="reg"]


dev.new()
par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
comm_master[,j={
	plot(plogis(detect_mu), reg_rich, main=reg[1], col=zCol(256,year))
},by="reg"]


# dev.new()
# par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
# comm_master[,j={
# 	plot(plogis(detect_mu_avg), naive_rich, main=reg[1], col=zCol(256,year))
# },by="reg"]

dev.new()
par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
spp_master[,j={
	dt <- .SD[present==1,list(detect_mu_avg=mean(detect_mu_avg, na.rm=TRUE), prevalence=(propStrata)),by=c("spp")]
	dt[,plot(prevalence, plogis(detect_mu_avg), main=reg[1])]
	abline(a=0,b=1,col='blue')
	NULL
},by=c("reg")]

dev.new()
par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
comm_master[,j={
	plot(propStrata_avg, reg_rich, main=reg[1], col=zCol(256,year))
},by=c("reg")]



# ---- annual Kmax ----
# load("trawlDiversity/pkgBuild/results/processedMsom/p.RData")
# u_regs <- regs
# annual_Kmax <- list()
# for(r in 1:length(u_regs)){
# 	annual_Kmax[[r]] <- p[[r]]$rd[,max(Kmax),by=c('stratum','year','reg')][,list(Kmax_mu=mean(V1)),by=c("reg","year")]
# }
# annual_Kmax <- rbindlist(annual_Kmax)
# annual_Kmax <- annual_Kmax[reg!='wcann']
# setkey(annual_Kmax, reg, year)
# dev.new()
# par(mfrow=c(3,3))
# annual_Kmax[,plot(year, Kmax_mu, type='o', main=reg[1]),by='reg']


# ---- detectability by c/e category ----
dev.new()
par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), cex=1, ps=8, mgp=c(0.75,0.1,0), tcl=-0.1)
spp_master[,j={.SD[!duplicated(paste(spp)),j={boxplot(plogis(detect_mu_avg)~ce_categ, ylab="plogis(detect_mu_avg)");NULL}];NULL},by=c("reg")]

# ---- detectability vs years since first observed ----
dev.new()
par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), cex=1, ps=8, mgp=c(0.75,0.1,0), tcl=-0.1)
spp_master[,,by='reg']
spp_master[,j={
	ylim <- .SD[,range(detect_mu)]
	xlim <- .SD[,range(yrs_since_1st)]
	plot(yrs_since_1st,detect_mu, main=reg[1], col=adjustcolor('black',0.05), pch=20)
	# .SD[,j={
	# 	if(.N>3){
	# 		lines(yrs_since_1st[is.finite(detect_mu)],fitted(lm(detect_mu~yrs_since_1st)), col=adjustcolor('blue',0.1))
	# 		NULL
	# 	}
	# },by='spp']
	detect_mods <- .SD[,j={
		if(.N>3){
			mod <- lm(detect_mu~yrs_since_1st)
			data.table(slope=coef(mod)[2], detect_mu_hat=fitted(mod), yrs_since_1st=yrs_since_1st[is.finite(detect_mu)])
		}
	},by='spp']
	detect_mods[,slope_col:=adjustcolor(zCol(256, slope),0.5)]
	# detect_mods[,slope_col:=adjustcolor(zCol(256, c(-max(abs(slope)),max(abs(slope)),slope))[-c(1:2)],0.5)]
	detect_mods[,lines(yrs_since_1st, detect_mu_hat, col=slope_col),by='spp']
},by='reg']



