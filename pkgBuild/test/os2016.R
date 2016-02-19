

library("trawlDiversity")
library("rbLib")
library("maps")

regs <- c("ebs", "ai", "goa", "wctri", "wcann", "gmex", "sa", "neus", "shelf", "newf")
pretty_reg <- c("ebs"="E. Bering Sea", "ai"="Aleutian Islands", "goa"="Gulf of Alaska", "wc"="West Coast US", "gmex"="Gulf of Mexico", "sa"="Southeast US", "neus"="Northeast US", "shelf"="Scotian Shelf", "newf"="Newfoundland")
pr <- names(pretty_reg)

pretty_col <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','navy','#a65628','salmon','#999999')
names(pretty_col) <- pr


# ====================================
# = Get Trimmed Data for Each Region =
# ====================================
data_in_regs <- list()
for(r in 1:length(regs)){
	t_reg <- regs[r]
	t_trimmed <- trim_msom(t_reg, gridSize=0.5, depthStratum=100, tolFraction=0.15, grid_stratum=TRUE, plot=FALSE)
	setkey(t_trimmed, year, stratum, haulid, spp)
	data_in_regs[[regs[r]]] <- t_trimmed
}


data_all <- rbindlist(data_in_regs, fill=TRUE)
data_all <- data_all[reg!="wcann"]
data_all[reg=="wctri", reg:="wc"]

data_all <- data_all[abund>0 & !is.na(abund)]


# =======================================
# = Get Bottom Temperature and Richness =
# =======================================
bt_all <- bt_metrics(pr, data_regs=data_all)
rich <- run_sac(pr, data_regs=data_all)
rich_bt <- merge(bt_all, rich, by=c("reg","year"), all=TRUE)

rich_bt <- rich_bt[reg!="wcann"]
rich_bt[reg=="wctri", reg:="wc"]

rich_bt[,reg_order:=as.numeric(factor(reg, levels=pr, ordered=TRUE))]
setorder(rich_bt, reg_order)
rich_bt[,reg_order:=NULL]


# =================================================
# = Temporal Patterns in Richness and Temperature =
# =================================================

# ---- Time Series of Richness: Separate Panels ----
dev.new()
par(mfrow=auto.mfrow(length(pr)), oma=c(0,0,1,0))
rich_bt[,plot(year, richness, type="o", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("Richness", side=3, line=-0.5, outer=TRUE)

# ---- Time Series of Mean Temperature: Separate Panels ----
dev.new(width=16, height=4)
# par(mfrow=auto.mfrow(length(pr)), oma=c(0,0,1,0))
par(mfrow=c(1,9), mar=c(0.75,0.75,0.25,0.25), oma=c(1,1,1,0), ps=8, cex=1, mgp=c(1,0.1,0), tcl=-0.1)
rich_bt[,
	j={
		plot(year, strat_mean, type="o", col=pretty_col[una(reg)], lwd=0.5, main=pretty_reg[una(reg)], xlab="", ylab="")
		lines(year, strat_mean_6yr_mean, col=pretty_col[una(reg)], lwd=1)
		lines(year, strat_mean_9yr_mean, col=pretty_col[una(reg)], lwd=1.5)
	}, by=c("reg")
]
mtext("Mean Bottom Temperature", side=3, line=-0.5, outer=TRUE)

# ---- Time Series of Temperature Spatial Standard Deviation: Separate Panels ----
dev.new()
par(mfrow=auto.mfrow(length(pr)), oma=c(0,0,1,0))
rich_bt[,plot(year, strat_sd, type="o", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("SD of Bottom Temperature", side=3, line=-0.5, outer=TRUE)

# ---- Time Series of Temperature Temporal Standard Deviation: Separate Panels ----
dev.new()
par(mfrow=auto.mfrow(length(pr)), oma=c(0,0,1,0))
rich_bt[,plot(year, strat_mean_9yr_sd, type="o", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("Temporal SD (9 yrs) of Bottom Temperature", side=3, line=-0.5, outer=TRUE)


# ================================================
# = Compare Richness and Mean Bottom Temperature =
# ================================================

# ---- Richness vs Mean Bottom Temp: Separate Panels ----
dev.new()
par(mfrow=auto.mfrow(length(pr)), oma=c(0,0,1,0))
rich_bt[,plot(strat_mean, richness, type="p", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("Richness vs Mean Temperature", side=3, line=-0.5, outer=TRUE)

# ---- Detrended Richness vs Detrended Mean Bottom Temp: Separate Panels ----
dev.new()
par(mfrow=auto.mfrow(length(pr)), oma=c(0,0,1,0))
rich_bt[,
	j={
		ttemp <- residuals(lm(strat_mean~year, na.action=na.exclude))
		trich <- residuals(lm(richness~year, na.action=na.exclude))
		plot(ttemp, trich, type="p", col=pretty_col[una(reg)], lwd=1, main=pretty_reg[una(reg)])

	}, by=c("reg")
]
mtext("Richness vs Mean Temperature (Detrended)", side=3, line=-0.5, outer=TRUE)

# ---- Regional Slopes from Richness vs Mean Bottom Temp ----
# dev.new(width=7, height=7)
png("~/Documents/School&Work/Presentations/OS2016/Figures/richness_temp_slopes.png", width=9, height=9, units="in", res=800)
par(mar=c(3.5,3.75,0.5,0.5), oma=c(0.2,0.2,0.2,0.2), ps=18, cex=1, mgp=c(2.5,0.65,0), tcl=-0.25)
reg_slope_bt_rich <- rich_bt[,list(
			bt_slope=coef(lm(strat_mean~year, na.action=na.exclude))[2],
			bt_slope_se=summary(lm(strat_mean~year, na.action=na.exclude))$coef[2,2],
			rich_slope=coef(lm(richness~year, na.action=na.exclude))[2],
			rich_slope_se=summary(lm(richness~year, na.action=na.exclude))$coef[2,2]
		), by=c("reg")
]
xlim <- reg_slope_bt_rich[,range(c(bt_slope-bt_slope_se, bt_slope+bt_slope_se))]
ylim <- reg_slope_bt_rich[,range(c(rich_slope-rich_slope_se, rich_slope+rich_slope_se))]
reg_slope_bt_rich[,plot(bt_slope, rich_slope, pch=20, cex=4, col=pretty_col[reg], xlim=xlim, ylim=ylim, ylab="Richness Trend", xlab="Temperature Trend")]
abline(v=0, lty="dashed")
abline(h=0, lty="dashed")
reg_slope_bt_rich[,arrows(x0=bt_slope-bt_slope_se, x1=bt_slope+bt_slope_se, y0=rich_slope, code=3, angle=90, length=0.1, col=pretty_col[reg])]
reg_slope_bt_rich[,arrows(x0=bt_slope, y0=rich_slope-rich_slope_se, y1=rich_slope+rich_slope_se, code=3, angle=90, length=0.1, col=pretty_col[reg])] 
dev.off()



# =============================================
# = Compare Richness and Temperature Variance =
# =============================================
# ---- Richness vs Spatial SD Bottom Temp: Separate Panels ----
dev.new()
par(mfrow=auto.mfrow(length(pr)), oma=c(0,0,1,0))
rich_bt[!is.na(strat_sd),plot(strat_sd, richness, type="p", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("Richness vs Spatial SD Temperature", side=3, line=-0.5, outer=TRUE)

dev.new()
par(mfrow=auto.mfrow(length(pr)), oma=c(0,0,1,0))
rich_bt[!is.na(strat_mean_6yr_sd),plot(strat_mean_6yr_sd, richness, type="p", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("Richness vs Temporal SD (6 yrs) Temperature", side=3, line=-0.5, outer=TRUE)

dev.new()
par(mfrow=auto.mfrow(length(pr)), oma=c(0,0,1,0))
rich_bt[!is.na(strat_mean_9yr_sd),plot(strat_mean_9yr_sd, richness, type="p", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("Richness vs Temporal SD (9 yrs) Temperature", side=3, line=-0.5, outer=TRUE)


# ---- Richness Anomaly vs Temporal SD (9 yrs) Bottom Temp: All Together ----
dev.new()
par(oma=c(0,0,0,0))
rich_bt[,richness_anomaly:=richness-mean(richness, na.rm=TRUE), by=c("reg")]
rich_bt[!is.na(strat_mean_9yr_sd),plot(strat_mean_9yr_sd, richness_anomaly, type="p", col=pretty_col[reg], pch=20, cex=1.25)]
rich_bt[!is.na(strat_mean_9yr_sd),j={
	v <- strat_mean_9yr_sd
	v2 <- strat_mean_9yr_sd^2
	mod <- lm(richness_anomaly~v+v2)
	vals <- seq(min(strat_mean_9yr_sd), max(strat_mean_9yr_sd), by=0.01)
	pred <- predict.lm(mod, newdata=data.frame(v=vals, v2=vals^2))
	# lwd <- ifelse(nrow(.SD)>=10, 2, 0.5)
	f_stat <- summary(mod)$fst
	lwd <- ifelse(pf(f_stat[1], f_stat[2], f_stat[3], lower.tail=FALSE) < 0.05, 5, 2)
	# lwd=4
	lines(vals, pred, col=pretty_col[reg], lwd=lwd)
}, by=c("reg")]
mtext("Richness Anomaly vs Temporal SD (9 yrs) Temperature", side=3, line=1, outer=F)




# ==================
# = Poster Figures =
# ==================

# ---- Map with Richness and Temperature Time Series ----
# Map
ll_all <- data_all[!duplicated(stratum),list(lon=mean(lon), lat=mean(lat), reg=reg), by=c("reg","stratum")]

nr = 60
nc = 150
l_mat <- matrix(NA, nrow=nr, nc=nc)
l_mat[c(1,nr,1,nr),c(1,1,nc,nc)] <- 1 # where map goes

l_mat[c(4, 19, 4, 19), c(5,5,20,20)] <- 2 # where ebs rich goes
l_mat[c(32,47,32,47), c(11,11,26,26)] <- 3 # where ai rich goes
l_mat[c(25,40,25,40), c(32,32,47,47)] <- 4 # where goa rich goes
l_mat[c(38,53,38,53), c(53,53,68,68)] <- 5 # where wc rich goes
l_mat[c(34, 49, 34, 49), c(80,80,95,95)] <- 6 # where gmex rich goes
l_mat[c(30, 45, 30, 45), c(98, 98, 113, 113)] <- 7 # where sa rich goes
l_mat[c(40, 55, 40, 55), c(123, 123, 138, 138)] <- 8 # where neus rich goes
l_mat[c(12, 27, 12, 27), c(110, 110, 125, 125)] <- 9 # where shelf rich goes
l_mat[c(3, 18, 3, 18), c(132,132,147,147)] <- 10 # where newf temp goes
l_mat[is.na(l_mat)] <- 0


# dev.new(width=29, height=10)
png("~/Documents/School&Work/Presentations/OS2016/Figures/map_insetRichnessTrend.png", width=29, height=10, units="in", res=800)
layout(l_mat)
par(mar=c(1.75,2,0.5,0.5), cex=1, ps=18, mgp=c(1,0.65,0), tcl=-0.25, bg="white", bty="o", lwd=2)
ylim <- range(ll_all[,lat]) + c(0, 8)
xlim <- range(ll_all[,lon]) + c(-1, 0)
ll_all[,plot(lon, lat, pch=19, col=pretty_col[reg], ylim=ylim, xlim=xlim, xaxt="s",yaxt="s", xlab="", ylab="", cex=4)]
map(add=TRUE, lwd=0.5)
par(mar=c(0.75,0.75,0.25,0.25), cex=1, ps=18, mgp=c(1,0.5,0), tcl=-0.25, bg="white", bty="o")

for(r in 1:9){
	rich_bt[reg==pr[r],plot(year, richness, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")]
	pu <- par("usr")	
	prop <- par("plt")
	
	dx <- diff(pu[1:2])/diff(prop[1:2])
	x0 <- pu[1] - prop[1]*dx
	x1 <- pu[2] + (1 - prop[2])*dx
	
	dy <- diff(pu[3:4])/diff(prop[3:4])
	y0 <- pu[3] - prop[3]*dy
	y1 <- pu[4] + (1 - prop[4])*dy
	
	if(pr[r]!="ebs"){rect(x0,y0,x1,y1,col = "white", border=NA, xpd=TRUE)}
	box(bty="l")
	grid(lty="solid",lwd=0.35)
	rich_bt[reg==pr[r],lines(year, richness, type="l", col=pretty_col[una(reg)], lwd=5)]
	axis(side=1, labels=TRUE)
	axis(side=2, labels=TRUE)
	# mtext("Richness", side=2, line=1)
}
dev.off()





# ---- Time series of temperature: spatial mean, spatial sd, temporal sd ----
# dev.new(width=23, height=5)
png("~/Documents/School&Work/Presentations/OS2016/Figures/bottom_temperature_time_series.png", width=23, height=5, units="in", res=800)
par(mfrow=c(3,9), mar=c(1.5,2.5,0.25,0.5), oma=c(1.0,1.5,1.5,0.1), ps=18, cex=1, mgp=c(1,0.65,0), tcl=-0.25)

rich_bt[, # temperature spatial mean
	j={
		plot(year, strat_mean, type="n", col=pretty_col[una(reg)], lwd=2, xlab="", ylab="", xaxt="n")
		axis(side=1, labels=FALSE)
		mtext(pretty_reg[una(reg)], side=3, line=0.1, font=2, cex=1.15)
		grid(lty="solid", lwd=0.35)
		lines(year, strat_mean, col=pretty_col[una(reg)], lwd=2)
		# lines(year, strat_mean_6yr_mean, col=pretty_col[una(reg)], lwd=1.25)
		# lines(year, strat_mean_9yr_mean, col=pretty_col[una(reg)], lwd=2)
		if(reg=="ebs"){mtext("Spatial Mean", side=2, line=2.5, outer=FALSE)}
	}, by=c("reg")
]

rich_bt[, # temperature spatial sd
	j={
		plot(year, strat_sd, type="n", col=pretty_col[una(reg)], lwd=2, xlab="", ylab="", xaxt="n")
		axis(side=1, labels=FALSE)
		grid(lty="solid",lwd=0.35)
		lines(year, strat_sd, col=pretty_col[una(reg)], lwd=2)
		if(reg=="ebs"){mtext("Spatial SD", side=2, line=2.5, outer=FALSE)}
	}, by=c("reg")
]

rich_bt[, # temperature spatial sd
	j={
		plot(year, strat_mean_9yr_sd, type="n", col=pretty_col[una(reg)], lwd=2, xlab="", ylab="")
		grid(lty="solid",lwd=0.35)
		lines(year, strat_mean_9yr_sd, col=pretty_col[una(reg)], lwd=2)
		if(reg=="ebs"){mtext("Temporal SD", side=2, line=2.5, outer=FALSE)}
	}, by=c("reg")
]
dev.off()


# ---- Time series of just mean BT ----
png("~/Documents/School&Work/Presentations/OS2016/Figures/single_bottom_temperature_time_series.png", width=29, height=3, units="in", res=800)
par(mfrow=c(1,9), mar=c(1.0,2.5,0.25,0.5), oma=c(1.0,1.5,1.5,0.1), ps=18, cex=1, mgp=c(1,0.65,0), tcl=-0.25)

rich_bt[, # temperature spatial mean
	j={
		plot(year, strat_mean, type="n", col=pretty_col[una(reg)], lwd=2, xlab="", ylab="", xaxt="s")
		axis(side=1, labels=FALSE)
		mtext(pretty_reg[una(reg)], side=3, line=0.25, font=2, cex=1.25)
		grid(lty="solid", lwd=0.35)
		lines(year, strat_mean, col=pretty_col[una(reg)], lwd=1, type="o", cex=0.75)
		lines(year, strat_mean_6yr_mean, col=pretty_col[una(reg)], lwd=1.5)
		lines(year, strat_mean_9yr_mean, col=pretty_col[una(reg)], lwd=2)
		if(reg=="ebs"){mtext("Temperature", side=2, line=2.5, outer=FALSE)}
	}, by=c("reg")
]
dev.off()

# ---- Standardized Richness vs Spatial Stdev. of Temperature ----
# dev.new(width=7, height=7)
png("~/Documents/School&Work/Presentations/OS2016/Figures/stdrz_rich_vs_spatial_stdev_temp.png", width=9, height=9, units="in", res=800)
par(mar=c(3.5,3.75,0.5,0.5), oma=c(0.2,0.2,0.2,0.2), ps=18, cex=1, mgp=c(2.5,0.65,0), tcl=-0.25)
rich_bt[,richness_standard:=richness_anomaly/sd(richness, na.rm=TRUE), by=c("reg")]
rich_bt[!is.na(strat_sd),plot(strat_sd, richness_standard, type="p", col=pretty_col[reg], pch=20, cex=1.25, xlab="Spatial Stdev. of Temperature", ylab="Standardized Richness")]
rich_bt[!is.na(strat_sd),j={
	v <- strat_sd
	# v2 <- strat_sd^2
	mod <- lm(richness_standard~v)
	# mod <- lm(richness_standard~v+v2)
	vals <- seq(min(strat_sd), max(strat_sd), by=0.01)
	pred <- predict.lm(mod, newdata=data.frame(v=vals, v2=vals^2))
	# lwd <- ifelse(nrow(.SD)>=10, 2, 0.5)
	f_stat <- summary(mod)$fst
	lwd <- ifelse(pf(f_stat[1], f_stat[2], f_stat[3], lower.tail=FALSE) < 0.05, 5, 2)
	# lwd=4
	lines(vals, pred, col=pretty_col[reg], lwd=lwd)
}, by=c("reg")]
dev.off()


# ---- Standardized Richness vs Temporal Stdev. of Temperature ----
# dev.new(width=7, height=7)
png("~/Documents/School&Work/Presentations/OS2016/Figures/stdrz_rich_vs_temporal_stdev_temp.png", width=9, height=9, units="in", res=800)
par(mar=c(3.5,3.75,0.5,0.5), oma=c(0.2,0.2,0.2,0.2), ps=18, cex=1, mgp=c(2.5,0.65,0), tcl=-0.25)
rich_bt[,richness_standard:=richness_anomaly/sd(richness, na.rm=TRUE), by=c("reg")]
rich_bt[!is.na(strat_mean_6yr_sd),plot(strat_mean_6yr_sd, richness_standard, type="p", col=pretty_col[reg], pch=20, cex=1.25, xlab="Temporal Stdev. of Temperature", ylab="Standardized Richness")]
rich_bt[!is.na(strat_mean_6yr_sd),j={
	v <- strat_mean_6yr_sd
	v2 <- strat_mean_6yr_sd^2
	mod <- lm(richness_standard~v+v2)
	vals <- seq(min(strat_mean_6yr_sd), max(strat_mean_6yr_sd), by=0.01)
	pred <- predict.lm(mod, newdata=data.frame(v=vals, v2=vals^2))
	# lwd <- ifelse(nrow(.SD)>=10, 2, 0.5)
	f_stat <- summary(mod)$fst
	lwd <- ifelse(pf(f_stat[1], f_stat[2], f_stat[3], lower.tail=FALSE) < 0.05, 5, 2)
	# lwd=4
	lines(vals, pred, col=pretty_col[reg], lwd=lwd)
}, by=c("reg")]
dev.off()



# ---- Hypothesis Figures ----
png("~/Documents/School&Work/Presentations/OS2016/Figures/hypothesis_niche_inheritence.png", width=5, height=5, units="in", res=600)
par(mar=c(3.5,3.75,0.5,0.5), oma=c(0.2,0.2,0.2,0.2), ps=18, cex=1, mgp=c(1.5,0.65,0), tcl=-0.25)
set.seed(1)
fx <- rnorm(100)
fy <- fx + rnorm(100, sd=0.15)
plot(fx, fy, type="n", pch=20, cex=1, xlab="Temperature or its Spatial Stdev.", ylab="Richness", xaxt="n", yaxt="n")
segments(x0=-2, y0=-2, x1=2, y1=2, lwd=3, lty="dashed")
dev.off()

png("~/Documents/School&Work/Presentations/OS2016/Figures/hypothesis_storage_effect.png", width=5, height=5, units="in", res=600)
par(mar=c(3.5,3.75,0.5,0.5), oma=c(0.2,0.2,0.2,0.2), ps=18, cex=1, mgp=c(1.5,0.65,0), tcl=-0.25)
set.seed(1)
fx <- sort(rnorm(100))
shift_fx <- min(fx)
fx <- fx - shift_fx
fx2 <- -(fx + shift_fx)^2
fy <- fx + fx2 + rnorm(100, sd=0.15)
plot(fx, fy, type="n", pch=20, cex=1, xlab="Temporal Stdev. of Temperature", ylab="Richness", xaxt="n", yaxt="n")
lines(fx, predict(lm(fy~fx+fx2)), lwd=3, lty="dashed")
dev.off()


# ==============
# = BS Checker =
# ==============
# ---- Time of Year Changing? ----
dev.new()
par(mfrow=auto.mfrow(data_all[,lu(reg)]))
data_all[, j={
	mu <- .SD[,list("Day of year"=mean(yday(datetime))),by="year"]
	ma <- .SD[,list("Day of year"=quantile(yday(datetime),0.95)),by="year"]
	mi <- .SD[,list("Day of year"=quantile(yday(datetime),0.05)),by="year"]
	ylim <- range(c(mu[[2]],ma[[2]],mi[[2]]))
	plot(mu, type="l", main=una(reg), ylim=ylim, lwd=2)
	lines(ma, type="l", lwd=1)
	lines(mi, type="l", lwd=1)
}, by="reg"]

dev.new(width=5, height=5)
# pdf("~/Desktop/trawl_sampling_timing.pdf",width=4, height=4)
par(mfrow=auto.mfrow(data_all[,lu(reg)]), mar=c(1.75,1.75,0.5,0.5), oma=c(0.1,0.1,1.5,0), mgp=c(0.75,0.1,0),tcl=-0.1, cex=1, ps=8)
zCol <- function(nCols, Z){
	cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(nCols)
	colVec_ind <- cut(Z, breaks=nCols)
	colVec <- cols[colVec_ind]
}
data_all[!duplicated(datetime), j={
	yd <- yday(datetime)
	mi <- min(yd)
	ma <- max(yd)
	d <- .SD[,density(yday(datetime), from=mi, to=ma)[c("x","y")],by="year"]
	ylim <- range(d[,y])

	uyrs <- d[,sort(una(year))]
	yr_cols <- adjustcolor(zCol(length(uyrs), uyrs), 0.75)
	d[year==uyrs[1], plot(x,y, ylim=ylim, xlab="Day of year", ylab="density", type="l", main=una(reg), col=yr_cols[1])]
	for(y in 2:length(uyrs)){
		tcol <- yr_cols[y]
		d[year==uyrs[y], lines(x, y, col=tcol)]
	}

}, by="reg"]
mtext("Distriubtion of Sampling Timing among Years (early years blue, late red)", side=3, outer=TRUE, line=0.5)
# dev.off()


# ---- Number of strata per year ----
dev.new()
par(mfrow=auto.mfrow(data_all[,lu(reg)]))
data_all[, j={
	plot(.SD[,list("N Strata"=lu(stratum)), by="year"], type="l", main=una(reg), lwd=2)
}, by="reg"]



