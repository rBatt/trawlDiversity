
library('maps')
library('raster')
library('spatstat')
library("spdep")
library('rbLib')
library('trawlDiversity')

# ===========
# = Generic =
# ===========



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
richness_ts()
dev.off()

# ---- prevalence vs time to absence ----
# dev.new(width=3.5, height=6)
# pdf("~/Desktop/Figure2_prevalence_absenceTime.pdf", width=3.5, height=3.5)
png("~/Desktop/Figure2_prevalence_absenceTime.png", width=3.5, height=6, units='in', res=200)
prevalence_absenceTime()
dev.off()

# ---- colonization rate map ----
png("~/Desktop/Figure3_cRate_map.png", width=7, height=3, units='in', res=200)
cRate_map()
dev.off()

# ---- richness vs detectability in 1 panel ----
# dev.new(width=3.5, height=3.5)
# pdf("~/Desktop/Figure4_rich_detect.pdf", width=3.5, height=3.5)
png("~/Desktop/Figure4_rich_detect.png", width=3.5, height=3.5, units='in', res=200)
rich_detect()
dev.off()



# ===================================
# = Manuscript Supplemental Figures =
# ===================================
# ---- Figure S1: Naive-MSOM Scatter ----
# png("~/Desktop/FigureS1_naive_msom_scatter.png", width=3.5, height=3.5, units='in', res=200)
qpng("FigureS1_naive_msom_scatter.png")
naive_msom_scatter()
dev.off()

# ---- Figure S2: Number of Species in Each C/E Category ----
# dev.new(width=5, height=3.5)
qpng("FigureS2_categ_barplot.png", width=5)
categ_barplot()
dev.off()

# ---- Figure S3: Colonization/ Extinction Time Series ----
# dev.new(width=3.5, height=3.5)
# png("~/Desktop/FigureS2_col_ext_ts.png", width=3.5, height=3.5, units='in', res=200)
qpng("FigureS3_col_ext_ts.png")
col_ext_ts()
dev.off()

# ---- Figure S4: Map of Neighs and Local AC of Coloniz. Rate ----
png("~/Desktop/FigureS4_nb_moranI.png", width=7, height=3, units='in', res=200)
nb_moranI()
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
