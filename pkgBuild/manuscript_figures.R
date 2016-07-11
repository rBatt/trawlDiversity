
library('maps')
library('raster')
library('spatstat')
library('rbLib')
library('trawlDiversity')

# ===========
# = Generic =
# ===========
regs <- comm_master[,una(reg)] #sapply(p, function(x)x$processed[,una(reg)])
pretty_reg <- c("ebs"="E. Bering Sea", "ai"="Aleutian Islands", "goa"="Gulf of Alaska", "wctri"="West Coast US", "gmex"="Gulf of Mexico", "sa"="Southeast US", "neus"="Northeast US", "shelf"="Scotian Shelf", "newf"="Newfoundland")
pretty_col <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','navy','#a65628','salmon','#999999')
names(pretty_col) <- names(pretty_reg)


# ==================
# = Handy Function =
# ==================
qpng <- function(name, width=3.5, height=3.5, res=200, location="~/Desktop"){
	png(file.path(location,name), width=width, height=height, units='in', res=res)
}


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
# dev.new(width=3.5, height=6)
# pdf("~/Desktop/Figure2_prevalence_absenceTime.pdf", width=3.5, height=3.5)
png("~/Desktop/Figure2_prevalence_absenceTime.png", width=3.5, height=6, units='in', res=200)
avg_prev_abs <- spp_master[!is.na(ext_dist) & ext_dist!=0,list(prevalence=mean(propStrata)),by=c("reg","ext_dist","stretch_type")]
par(mfrow=c(2,1), mar=c(1.5,1.5,0.25,0.25), oma=c(0.1,0.1,0.1,0.1), cex=1, ps=8, mgp=c(0.75,0.1,0), tcl=-0.1, ylbias=0.35)
for(st in 1:2){
	t_st <- c("pre_ext","post_col")[st]
	t_xlab <- c("Years until extinction", "Years after colonization")[st]
	t_panel <- c("A","B")[st]
	leg_log <- c(TRUE, FALSE)[st]
	
	avg_prev_abs[stretch_type==t_st,plot(ext_dist, prevalence, col=adjustcolor(pretty_col[(reg)], 0.5), pch=16, xlab=t_xlab, ylab="Prevalence")]
	avg_prev_abs[stretch_type==t_st,j={
		lines(ext_dist, fitted(lm(prevalence~ext_dist)), col=pretty_col[(reg[1])], lwd=1.5)
	},by=c('reg')]
	
	legend("topleft", legend=t_panel, inset=c(-0.075, -0.03), bty='n', text.font=2)
	if(leg_log){
		comm_master[,legend("topright",ncol=1,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.02, -0.03), bty='n', x.intersp=1, y.intersp=0.65)]
	}	
}
dev.off()

# ---- colonization rate map ----
png("~/Desktop/Figure3_cRate_map.png", width=7, height=3, units='in', res=200)
map_layout <- trawl_layout()
par(mar=c(1.15,1.15,0.25,0.25), mgp=c(0.5,0.075,0), tcl=-0.1, ps=8, cex=1, oma=c(0.75,0.5,1,0.1))
layout(map_layout)

map_col <- grDevices::colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
toRast <- function(p){
	x = rep(p$xcol, length(p$yrow))
	y = rep(p$yrow, each=length(p$xcol))
	z = c(t(as.matrix(p)))
	raster::rasterFromXYZ(cbind(x, y, z))
}

u_regs <- mapDat[,unique(reg)]
rs <- X[,una(reg)]
nr <- length(rs)
mapPPP_col <- list()
for(r in 1:nr){
	td <- X[reg==rs[r]]
	mapPPP_col[[r]] <- spatstat::ppp(x=td[,lon], y=td[,lat], marks=td[,n_spp_col_weighted], window=mapOwin[[r]]) # /avgRich
	
	t_idw <- spatstat::Smooth(mapPPP_col[[r]], hmax=1)
	z <- toRast(t_idw)
	image(z, col=map_col, xlab="", ylab="")
	map(add=TRUE, fill=TRUE, col="slategray")
	
	zl <- range(values(z)*10, na.rm=TRUE)
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
	switch(rs[r],
		ebs = legend("topright", legend="A", bty='n', text.font=2, inset=c(-0.02,-0.15), cex=1.25, text.col='white'),
		ai = legend("topleft", legend="C", bty='n', text.font=2, inset=c(-0.065,-0.45), cex=1.25, xpd=T),
		goa = legend("topleft", legend="B", bty='n', text.font=2, inset=c(-0.065,-0.06), cex=1.25),
		wctri = legend("top", legend="E", bty='n', text.font=2, inset=c(0,0.15), cex=1.25, text.col='white'),
		gmex = legend("topleft", legend="G", bty='n', text.font=2, inset=c(-0.175,-0.12), cex=1.25, text.col='white'),
		sa = legend("topleft", legend="I", bty='n', text.font=2, inset=c(-0.15,-0.075), cex=1.25, text.col='white'),
		neus = legend("topleft", legend="H", bty='n', text.font=2, inset=c(-0.125,-0.05), cex=1.25, text.col='white'),
		shelf = legend("topleft", legend="D", bty='n', text.font=2, inset=c(-0.1,-0.125), cex=1.25, text.col='white'),
		newf = legend("topright", legend="F", bty='n', text.font=2, inset=c(-0.01,-0.05), cex=1.25)
	)
}
mtext(bquote(Colonization~Rate~(C[w]~~decade^-1)), side=3, outer=TRUE, font=2, line=-0.3)
mtext(bquote(Longitude~(phantom()*degree*E)), side=1, line=-0.15, outer=TRUE)
mtext(bquote(Latitude~(phantom()*degree*N)), side=2, line=-0.75, outer=TRUE)
dev.off()

# ---- richness vs detectability in 1 panel ----
# dev.new(width=3.5, height=3.5)
# pdf("~/Desktop/Figure4_rich_detect.pdf", width=3.5, height=3.5)
png("~/Desktop/Figure4_rich_detect.png", width=3.5, height=3.5, units='in', res=200)
par(mar=c(1.75,1.5,0.25,0.25),mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
comm_master[,plot((propTow_occ_avg), reg_rich, col=adjustcolor(pretty_col[reg],0.5), xlab="Within-site prevalence", ylab="Estimated richness", pch=16)]
comm_master[,lines((propTow_occ_avg),fitted(lm(reg_rich~propTow_occ_avg))),by='reg']
comm_master[,legend("topright",ncol=2,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.02, -0.02), bty='n')]
dev.off()



# ===================================
# = Manuscript Supplemental Figures =
# ===================================
# ---- Figure S1: Naive-MSOM Scatter ----
# png("~/Desktop/FigureS1_naive_msom_scatter.png", width=3.5, height=3.5, units='in', res=200)
qpng("FigureS1_naive_msom_scatter.png")
par(mfrow=c(3,3), mar=c(1.25,1.0,0.5,0.1), oma=c(0.35,0.5,0.1,0.1), mgp=c(0.25,0.1,0), tcl=-0.1, ps=8, cex=1)
comm_master[,j={
	plot(naive_rich, reg_rich, main=pretty_reg[una(reg)], type='p', xlab="", ylab="", cex.main=1)
	abline(a=0, b=1)
}, by='reg']
mtext("MSOM Richness", side=2, line=-0.2, outer=TRUE)
mtext("Naive Richness", side=1, line=-0.5, outer=TRUE)
dev.off()

# ---- Figure S2: Number of Species in Each C/E Category ----
categ_table <- t(spp_master[!duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ)])[c(4,1,2,3),]
colnames(categ_table) <- pretty_reg[colnames(categ_table)]
colnames(categ_table) <- gsub("^(.*) (.*)$", "\\1\n\\2", colnames(categ_table))
# dev.new(width=5, height=3.5)
qpng("FigureS2_categ_barplot.png", width=5)
par(cex=1, mar=c(3,2,1,0.1), ps=8)
bp <- barplot(categ_table, beside=T, legend=T, names.arg=rep("",ncol(categ_table)), args.legend=list(bty='n'))
text(colMeans(bp)-1, -11, labels=colnames(categ_table), srt=45, xpd=TRUE)
dev.off()

# ---- Figure S3: Colonization/ Extinction Time Series ----
# dev.new(width=3.5, height=3.5)
# png("~/Desktop/FigureS2_col_ext_ts.png", width=3.5, height=3.5, units='in', res=200)
qpng("FigureS3_col_ext_ts.png")
par(mfrow=c(3,3), mar=c(1.25,1.0,0.5,0.1), oma=c(0.35,0.5,0.1,0.1), mgp=c(0.25,0.1,0), tcl=-0.1, ps=8, cex=1)
comm_master[,j={
	ylim=range(c(n_col,n_ext));
	plot(year, n_col, main=pretty_reg[una(reg)], type='l', col='blue', ylim=ylim, xlab="", ylab="", cex.main=1)
	lines(year, n_ext, col='red')
}, by='reg']
mtext("Colonizations or Extinctions", side=2, line=-0.2, outer=TRUE)
mtext("Year", side=1, line=-0.5, outer=TRUE)
dev.off()

# ---- Figure S4: Map of Neighs and Local AC of Coloniz. Rate ----
png("~/Desktop/FigureS4_nb_moranI.png", width=7, height=3, units='in', res=200)
map_layout <- trawl_layout()
par(mar=c(0.25,0.25,0.25,0.25), mgp=c(0.25,0.075,0), tcl=-0.1, ps=8, cex=1, oma=c(0.1,0.1,0.1,0.1))
layout(map_layout)
pretty_reg <- c("ebs"="E. Bering Sea", "ai"="Aleutian Islands", "goa"="Gulf of Alaska", "wctri"="West\nCoast\nUS", "gmex"="Gulf of Mexico", "sa"="Southeast US", "neus"="Northeast US", "shelf"="Scotian Shelf", "newf"="Newfoundland")

rs <- mapDat[,una(reg)]
nr <- length(rs)
for(r in 1:nr){
	t_lac <- localAC[[rs[r]]]
	plot(mapOwin[[rs[r]]], coords=t_lac$I[,list(lon,lat)], add=FALSE, main="")
	box()
	if(rs[r]=='wctri'){
		mtext(pretty_reg[rs[r]], side=3, line=-3)
	}else{
		mtext(pretty_reg[rs[r]], side=3, line=-1)
	}
	plot(t_lac$nb, t_lac$I[,list(lon,lat)], add=TRUE, col=adjustcolor('black', 0.5), cex=0.8, lwd=0.5)
	sig_lac <- t_lac$I[,lI_pvalue]<0.05
	if(any(sig_lac)){
		zl <- range(t_lac$I[sig_lac,Ii], na.rm=TRUE)
		t_col <- rbLib::zCol(256,t_lac$I[sig_lac,Ii])
		map_col <- rbLib::zCol(6, 1:6)
		points(x=t_lac$I[sig_lac,lon], y=t_lac$I[sig_lac,lat], bg=t_col, pch=21, cex=1.1)
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
}
dev.off()




# ===============================
# = ***Select*** Backup Figures =
# ===============================
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
