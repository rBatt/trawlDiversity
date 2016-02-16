library("trawlDiversity")
library("rbLib")

regs <- c("ebs", "ai", "goa", "wctri", "wcann", "gmex", "sa", "neus", "shelf", "newf")


# ====================================
# = Get Trimmed Data for Each Region =
# ====================================
data_in_regs <- list()
for(r in 1:length(regs)){

	t_reg <- regs[r]

	data_in_all0 <- trim_msom(t_reg, gridSize=0.5, depthStratum=100, tolFraction=0.15, grid_stratum=TRUE, plot=FALSE)

	nSite_min <- data_in_all0[,lu(stratum), by="year"][,min(V1)]

	data_in_all <- data_in_all0
	setkey(data_in_all, year, stratum, haulid, spp)
	
	data_in_regs[[regs[r]]] <- data_in_all
}


# =======================================
# = Get Bottom Temperature and Richness =
# =======================================
bt_all <- bt_metrics(regs, data_regs=data_in_regs)
rich <- run_sac(regs, data_regs=data_in_regs)
rich_bt <- merge(bt_all, rich, by=c("reg","year"), all=TRUE)

rich_bt <- rich_bt[reg!="wcann"]
rich_bt[reg=="wctri", reg:="wc"]


# =================================================
# = Temporal Patterns in Richness and Temperature =
# =================================================
pretty_reg <- c("ebs"="E. Bering Sea", "ai"="Aleutian Islands", "goa"="Gulf of Alaska", "wc"="West Coast US", "gmex"="Gulf of Mexico", "sa"="Southeast US", "neus"="Northeast US", "shelf"="Scotian Shelf", "newf"="Newfoundland")
rn <- names(pretty_reg)
pretty_col <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','navy','#a65628','salmon','#999999')
names(pretty_col) <- rn


# ---- Time Series of Richness: Separate Panels ----
dev.new()
par(mfrow=auto.mfrow(length(rn)), oma=c(0,0,1,0))
rich_bt[,plot(year, richness, type="o", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("Richness", side=3, line=-0.5, outer=TRUE)


# ---- Time Series of Mean Temperature: Separate Panels ----
dev.new()
par(mfrow=auto.mfrow(length(rn)), oma=c(0,0,1,0))
rich_bt[,
	j={
		plot(year, strat_mean, type="o", col=pretty_col[una(reg)], lwd=0.5, main=pretty_reg[una(reg)])
		lines(year, strat_mean_6yr_mean, col=pretty_col[una(reg)], lwd=1)
		lines(year, strat_mean_9yr_mean, col=pretty_col[una(reg)], lwd=1.5)
	}, by=c("reg")
]
mtext("Mean Bottom Temperature", side=3, line=-0.5, outer=TRUE)


# ---- Time Series of Temperature Standard Deviation: Separate Panels ----
dev.new()
par(mfrow=auto.mfrow(length(rn)), oma=c(0,0,1,0))
rich_bt[,plot(year, strat_sd, type="o", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("SD of Bottom Temperature", side=3, line=-0.5, outer=TRUE)



# ---- Time Series of Temperature Rolling Window Standard Deviation: Separate Panels ----
dev.new()
par(mfrow=auto.mfrow(length(rn)), oma=c(0,0,1,0))
rich_bt[,plot(year, strat_mean_6yr_sd, type="o", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("Temporal SD (6 yrs) of Bottom Temperature", side=3, line=-0.5, outer=TRUE)

dev.new()
par(mfrow=auto.mfrow(length(rn)), oma=c(0,0,1,0))
rich_bt[,plot(year, strat_mean_9yr_sd, type="o", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("Temporal SD (9 yrs) of Bottom Temperature", side=3, line=-0.5, outer=TRUE)




# ================================================
# = Compare Richness and Mean Bottom Temperature =
# ================================================

# ---- Richness vs Mean Bottom Temp: Separate Panels ----
dev.new()
par(mfrow=auto.mfrow(length(rn)), oma=c(0,0,1,0))
rich_bt[,plot(strat_mean, richness, type="p", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("Richness vs Mean Temperature", side=3, line=-0.5, outer=TRUE)


# ---- Detrended Richness vs Detrended Mean Bottom Temp: Separate Panels ----
dev.new()
par(mfrow=auto.mfrow(length(rn)), oma=c(0,0,1,0))
rich_bt[,
	j={
		ttemp <- residuals(lm(strat_mean~year, na.action=na.exclude))
		trich <- residuals(lm(richness~year, na.action=na.exclude))
		plot(ttemp, trich, type="o", col=pretty_col[una(reg)], lwd=1, main=pretty_reg[una(reg)])
		
	}, by=c("reg")
]
mtext("Richness vs Mean Temperature (Detrended)", side=3, line=-0.5, outer=TRUE)


# ---- Richness vs Mean Bottom Temp: All Together ----
dev.new()
par(oma=c(0,0,0,0))
rich_bt[,plot(strat_mean, richness, type="p", col=pretty_col[reg])]
rich_bt[,lines(strat_mean, predict(lm(richness~strat_mean)), col=pretty_col[una(reg)]), by=c("reg")]
mtext("Richness vs Mean Temperature", side=3, line=1, outer=F)


# ---- Detrended Richness vs Detrended Mean Bottom Temp: All Together ----
dev.new()
par(oma=c(0,0,0,0))
rich_bt[,c("det_strat_mean","det_rich") := list(
			residuals(lm(strat_mean~year, na.action=na.exclude)),
			residuals(lm(richness~year, na.action=na.exclude))
		), by=c("reg")
]
rich_bt[,plot(det_strat_mean, det_rich, type="p", col=pretty_col[reg])]
rich_bt[,lines(det_strat_mean, predict(lm(det_rich~det_strat_mean)), col=pretty_col[una(reg)]), by=c("reg")]
mtext("Richness vs Mean Temperature (Detrended)", side=3, line=1, outer=F)


# ---- Regional Slopes from Richness vs Mean Bottom Temp ----
dev.new()
par(oma=c(0,0,0,0))
reg_slope_bt_rich <- rich_bt[,list(
			bt_slope=coef(lm(strat_mean~year, na.action=na.exclude))[2],
			bt_slope_se=summary(lm(strat_mean~year, na.action=na.exclude))$coef[2,2],
			rich_slope=coef(lm(richness~year, na.action=na.exclude))[2],
			rich_slope_se=summary(lm(richness~year, na.action=na.exclude))$coef[2,2]
		), by=c("reg")
]
xlim <- reg_slope_bt_rich[,range(c(bt_slope-bt_slope_se, bt_slope+bt_slope_se))]
ylim <- reg_slope_bt_rich[,range(c(rich_slope-rich_slope_se, rich_slope+rich_slope_se))]
reg_slope_bt_rich[,plot(bt_slope, rich_slope, pch=20, cex=3, col=pretty_col[reg], xlim=xlim, ylim=ylim)]
abline(v=0, lty="dashed")
abline(h=0, lty="dashed")
reg_slope_bt_rich[,arrows(x0=bt_slope-bt_slope_se, x1=bt_slope+bt_slope_se, y0=rich_slope, code=3, angle=90, length=0.1, col=pretty_col[reg])]
reg_slope_bt_rich[,arrows(x0=bt_slope, y0=rich_slope-rich_slope_se, y1=rich_slope+rich_slope_se, code=3, angle=90, length=0.1, col=pretty_col[reg])] 

mtext("Richness Slope vs Temperature Slope", side=3, line=1, outer=F)




# =============================================
# = Compare Richness and Temperature Variance =
# =============================================
# ---- Richness vs Temporal SD Bottom Temp: Separate Panels ----
dev.new()
par(mfrow=auto.mfrow(length(rn)), oma=c(0,0,1,0))
rich_bt[!is.na(strat_mean_6yr_sd),plot(strat_mean_6yr_sd, richness, type="p", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("Richness vs Temporal SD (6 yrs) Temperature", side=3, line=-0.5, outer=TRUE)

dev.new()
par(mfrow=auto.mfrow(length(rn)), oma=c(0,0,1,0))
rich_bt[!is.na(strat_mean_9yr_sd),plot(strat_mean_9yr_sd, richness, type="p", col=pretty_col[una(reg)], main=pretty_reg[una(reg)]), by=c("reg")]
mtext("Richness vs Temporal SD (9 yrs) Temperature", side=3, line=-0.5, outer=TRUE)


# ---- Richness vs Temporal SD (6 yrs) Bottom Temp: All Together ----
dev.new()
par(oma=c(0,0,0,0))
rich_bt[!is.na(strat_mean_6yr_sd),plot(strat_mean_6yr_sd, richness, type="p", col=pretty_col[reg])]
rich_bt[!is.na(strat_mean_6yr_sd),j={
	v <- strat_mean_6yr_sd
	v2 <- strat_mean_6yr_sd^2
	mod <- lm(richness~v+v2)
	vals <- seq(min(strat_mean_6yr_sd), max(strat_mean_6yr_sd), by=0.01)
	pred <- predict.lm(mod, newdata=data.frame(v=vals, v2=vals^2))
	lwd <- ifelse(nrow(.SD)>=10, 2, 0.5)
	lines(vals, pred, col=pretty_col[reg], lwd=lwd)
}, by=c("reg")]
mtext("Richness vs Temporal SD (6 yrs) Temperature", side=3, line=1, outer=F)


# ---- Richness Anomaly vs Temporal SD (6 yrs) Bottom Temp: All Together ----
dev.new()
par(oma=c(0,0,0,0))
rich_bt[,richness_anomaly:=richness-mean(richness, na.rm=TRUE), by=c("reg")]
rich_bt[!is.na(strat_mean_6yr_sd),plot(strat_mean_6yr_sd, richness_anomaly, type="p", col=pretty_col[reg])]
rich_bt[!is.na(strat_mean_6yr_sd),j={
	v <- strat_mean_6yr_sd
	v2 <- strat_mean_6yr_sd^2
	mod <- lm(richness_anomaly~v+v2)
	vals <- seq(min(strat_mean_6yr_sd), max(strat_mean_6yr_sd), by=0.01)
	pred <- predict.lm(mod, newdata=data.frame(v=vals, v2=vals^2))
	lwd <- ifelse(nrow(.SD)>=10, 2, 0.5)
	lines(vals, pred, col=pretty_col[reg], lwd=lwd)
}, by=c("reg")]
mtext("Richness Anomaly vs Temporal SD (6 yrs) Temperature", side=3, line=1, outer=F)

# ---- Richness vs Temporal SD (9 yrs) Bottom Temp: All Together ----
dev.new()
par(oma=c(0,0,0,0))
rich_bt[!is.na(strat_mean_9yr_sd),plot(strat_mean_9yr_sd, richness, type="p", col=pretty_col[reg])]
rich_bt[!is.na(strat_mean_9yr_sd),j={
	v <- strat_mean_9yr_sd
	v2 <- strat_mean_9yr_sd^2
	mod <- lm(richness~v+v2)
	vals <- seq(min(strat_mean_9yr_sd), max(strat_mean_9yr_sd), by=0.01)
	pred <- predict.lm(mod, newdata=data.frame(v=vals, v2=vals^2))
	lwd <- ifelse(nrow(.SD)>=10, 2, 0.5)
	lines(vals, pred, col=pretty_col[reg], lwd=lwd)
}, by=c("reg")]
mtext("Richness vs Temporal SD (9 yrs) Temperature", side=3, line=1, outer=F)


# ---- Richness Anomaly vs Temporal SD (9 yrs) Bottom Temp: All Together ----
dev.new()
par(oma=c(0,0,0,0))
rich_bt[,richness_anomaly:=richness-mean(richness, na.rm=TRUE), by=c("reg")]
rich_bt[!is.na(strat_mean_9yr_sd),plot(strat_mean_9yr_sd, richness_anomaly, type="p", col=pretty_col[reg])]
rich_bt[!is.na(strat_mean_9yr_sd),j={
	v <- strat_mean_9yr_sd
	v2 <- strat_mean_9yr_sd^2
	mod <- lm(richness_anomaly~v+v2)
	vals <- seq(min(strat_mean_9yr_sd), max(strat_mean_9yr_sd), by=0.01)
	pred <- predict.lm(mod, newdata=data.frame(v=vals, v2=vals^2))
	# lwd <- ifelse(nrow(.SD)>=10, 2, 0.5)
	# f_stat <- summary(mod)$fst
	# lwd <- ifelse(pf(f_stat[1], f_stat[2], f_stat[3], lower.tail=FALSE) < 0.05, 2, 0.5)
	lwd=2
	lines(vals, pred, col=pretty_col[reg], lwd=lwd)
}, by=c("reg")]
mtext("Richness Anomaly vs Temporal SD (9 yrs) Temperature", side=3, line=1, outer=F)







