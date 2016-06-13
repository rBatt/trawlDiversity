library('trawlDiversity')
library('maps')

# ===========
# = Generic =
# ===========
regs <- comm_master[,una(reg)] #sapply(p, function(x)x$processed[,una(reg)])
pretty_reg <- c("ebs"="E. Bering Sea", "ai"="Aleutian Islands", "goa"="Gulf of Alaska", "wctri"="West Coast US", "gmex"="Gulf of Mexico", "sa"="Southeast US", "neus"="Northeast US", "shelf"="Scotian Shelf", "newf"="Newfoundland")
pretty_col <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','navy','#a65628','salmon','#999999')
names(pretty_col) <- names(pretty_reg)


# ===================================
# = Manuscript Main Text Candidates =
# ===================================
# ---- richness time series in 1 panel ----
# dev.new(width=3.5, height=5)
# pdf("~/Desktop/Figure1_richness_ts.pdf", width=3.5, height=5)
png("~/Desktop/Figure1_richness_ts.png", width=3.5, height=5, units='in', res=200)
par(mar=c(1.75,1.5,0.25,0.25),mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
comm_master[,plot(year, reg_rich, col=pretty_col[reg], xlab="Year", ylab="Estimated Richness", pch=20)]
comm_master[,lines(year, reg_rich, lwd=0.5, col=adjustcolor(pretty_col[reg], 0.5)),by='reg']
comm_master[,lines(year,fitted(lm(reg_rich~year))),by='reg']
comm_master[,legend("topleft",ncol=1,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.085,-0.01), bty='n')]
dev.off()

# ---- prevalence vs time to absence ----
# dev.new(width=3.5, height=3.5)
# pdf("~/Desktop/Figure2_prevalence_absenceTime.pdf", width=3.5, height=3.5)
png("~/Desktop/Figure2_prevalence_absenceTime.png", width=3.5, height=3.5, units='in', res=200)
avg_prev_abs <- spp_master[!is.na(ext_dist_sign) & ext_dist_sign!=0,list(prevalence=mean(propStrata)),by=c("reg","ext_dist_sign","stretch_type")]
par(mar=c(1.5,1.5,0.25,0.25), oma=c(0.1,0.1,0.1,0.1), cex=1, ps=8, mgp=c(0.75,0.1,0), tcl=-0.1, ylbias=0.35)
avg_prev_abs[,plot(ext_dist_sign, prevalence, col=adjustcolor(pretty_col[(reg)], 0.5), pch=16, xlab="Years until extinction", ylab="Prevalence")]
avg_prev_abs[,j={
	lines(ext_dist_sign, fitted(lm(prevalence~ext_dist_sign)), col=pretty_col[(reg[1])], lwd=1.5)
},by=c('reg',"stretch_type")]
comm_master[,legend("topleft",ncol=2,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.075, -0.02), bty='n')]
dev.off()

# ---- colonization rate map ----
# dev.new(height=3, width=7)
# pdf("~/Desktop/Figure3_cRate_map.pdf", width=7, height=3)
png("~/Desktop/Figure3_cRate_map.png", width=7, height=3, units='in', res=200)
map_layout <- trawl_layout()
par(mar=c(1.5,1.5,0.5,0.5), mgp=c(0.75,0.1,0), tcl=-0.1,ps=8, cex=1, oma=c(0.5,0.5,1,0.1))
layout(map_layout)
u_regs <- mapDat[,unique(reg)]
for(r in 1:lu(u_regs)){
	mapDat[reg==u_regs[r], plot_space(lon,lat, n_spp_col_weighted/avgRich, bty='l', xlab="", ylab="")]
	map(add=TRUE, fill=TRUE, col="white")
}
mtext(bquote(Specific~Colonization~Rate~(C~~y^-1~spp^-1)), side=3, outer=TRUE, font=2, line=-0.5)
mtext(bquote(Longitude~(phantom()*degree*E)), side=1, line=-0.4, outer=TRUE)
mtext(bquote(Latitude~(phantom()*degree*N)), side=2, line=-0.75, outer=TRUE)
dev.off()

# ---- richness vs detectability in 1 panel ----
# dev.new(width=3.5, height=3.5)
# pdf("~/Desktop/Figure4_rich_detect.pdf", width=3.5, height=3.5)
png("~/Desktop/Figure4_rich_detect.png", width=3.5, height=3.5, units='in', res=200)
par(mar=c(1.75,1.5,0.25,0.25),mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
comm_master[,plot(plogis(detect_mu_avg), reg_rich, col=adjustcolor(pretty_col[reg],0.5), xlab="Detectability", ylab="Estimated Richness", pch=16)]
comm_master[,lines(plogis(detect_mu_avg),fitted(lm(reg_rich~plogis(detect_mu_avg)))),by='reg']
comm_master[,legend("topright",ncol=2,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.02, -0.02), bty='n')]
dev.off()






# ---- prevalence vs time to absence (panel per region, non-average) ----
dev.new(width=3.5, height=3.5)
par(mfrow=c(3,3), mar=c(0.75,0.75,0.25,0.25), oma=c(0.75,0.75,0.1,0.1), cex=1, ps=8, mgp=c(0.05,0.01,0), tcl=-0.1, ylbias=0.35)
spp_master <- copy(trawlDiversity::spp_master)
spp_master[,ext_dist_sign2:=-ext_dist_sign]
spp_master[!is.na(stretch_type) & is.finite(ext_dist_sign2) & now_ext!=1,j={
	ylim <- .SD[,range(propStrata)]
	xlim <- .SD[,range(ext_dist_sign2)]
	propStrata_mods <- .SD[,j={
		if(.N>4 & lu(ext_dist_sign2)>=3){
			mod <- lm(propStrata~ext_dist_sign2)
			# print(coef(summary(mod)))
# 			if(nrow(coef(summary(mod)))==1){
# 				blah <<- mod
# 			}
			data.table(slope=coef(mod)[2], pval=coef(summary(mod))[2,4], propStrata=propStrata, propStrata_hat=fitted(mod), ext_dist_sign2=ext_dist_sign2[is.finite(propStrata)])
		}
	},by=c('spp','stretch_type')]
	# propStrata_mods[,pval_col:=adjustcolor(zCol(256, pval),0.5)]
	propStrata_mods[,pval_col:=ifelse(pval<=0.05, 'blue', 'black')]
	propStrata_mods[stretch_type=="pre_ext" & pval_col=='blue',pval_col:='red']
	setorder(propStrata_mods, -pval)
	propStrata_mods[,plot(ext_dist_sign2,propStrata, col=adjustcolor('black',0.05), pch=16, xlab="", ylab="")]
	mtext(pretty_reg[una(reg)], side=3, line=-0.75, adj=ifelse(reg[1]=="sa",0.1, 0.9), font=2, cex=1)
	propStrata_mods[,lines(ext_dist_sign2, propStrata_hat, col=pval_col),by=c('spp','stretch_type')]
},by='reg']
mtext("Years after extinction", side=1, line=0, outer=TRUE)
mtext("Prevalence", side=2, line=-0.05, outer=TRUE)
spp_master[,ext_dist_sign2:=NULL]
