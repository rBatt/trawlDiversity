
library("rstan")
library("trawlDiversity")
library("rbLib")
library("R2jags")
library("maps")
library("beanplot")
library("fields")

setwd("~/Documents/School&Work/pinskyPost/trawl")

load("trawlDiversity/pkgBuild/results/processedMsom/p.RData")

Figures <- list()


for(reg_num in 8:length(p)){
	

	t_prn <- p[[reg_num]]

	# ===========
	# = Figures =
	# ===========

	# ---- Figure 1 ----
	Figures <- plot_rich_bt_ts(t_prn, Figures)
	# dev.off()

	# ---- Figure 2 ----
	Figures <- plot_rich_bt_scatter(t_prn, Figures)
	# dev.off()

	# ---- Figure 3 ----
	Figures <- plot_btemp_map(t_prn, Figures)
	# dev.off()
	
	# ---- Figure 4 ----
	Figures <- plot_traceplot(t_prn, Figures)
	# dev.off()


	# ---- Figure 5 ----
	Figures <- plot_post_corr(t_prn, Figures, yr=1)
	# dev.off()


	# ---- Figure 6 ----
	Figures <- plot_col_vs_unobsSpp(t_prn, Figures)
	# dev.off()


	# ---- Figure 7 ----
	# ---- Number of Colonizations per Stratum ----	
	Figures <- plot_colExt_perStrat(t_prn, Figures)

	# ---- Figure 8 ----
	# ---- Plot Information and Identity of Colonizers, Leavers, etc ----
	Figures <- plot_ce_wrap(t_prn, Figures, spp_cat="col", width.max=12, max_spp_columns=12)
	Figures <- plot_ce_wrap(t_prn, Figures, spp_cat="ext", width.max=12, max_spp_columns=12)
	
	
	# ---- Figure 9:  ----
	# ---- Plot Number of Colonizers, Leavers, and Plot Temp Rank ----
	Figures <- plot_rank_temp(t_prn, Figures)
	
	

	# # ---- Figure 9 ----
	# # ---- Time Series of Species Detection Parameter (Beta) ----
	# fig9_name <- paste0("beta_spp_timeSeries", reg, ".png")
	# sppPar_mfr <- auto.mfrow(ab[,lu(spp)], tall=TRUE)
	# sP_w <- 8
	# sP_h <- sppPar_mfr[1] * sP_w / sppPar_mfr[2] * 1.1 # bonus height b/c of extra top margin space
	# fig9_dim <- c(width=sP_w, height=sP_h)
	#
	# dev.new(width=fig9_dim[1], height=fig9_dim[2])
	# par(mfrow=sppPar_mfr, mar=c(1,1,1,0.1), oma=c(0.1,0.1,1.5,0.1), ps=6, mgp=c(0.6,0.1,0), tcl=-0.1, cex=1)
	#
	# u_spp <- ab[,una(spp)]
	# for(s in 1:n_spp){
	# 	t_s <- u_spp[s]
	# 	t_dat <- ab[par=="beta" & (spp)==t_s]
	# 	plot_ab(X=t_dat, t_spp=t_s)
	# }
	# mtext(paste0(reg, " Beta[1] (detection intercept)"), outer=TRUE, side=3, font=2, line=0.25, cex=1.25)
	#
	# Figures[[reg]][['Figure9']][["figure"]] <- recordPlot()
	# Figures[[reg]][['Figure9']][["name"]] <- fig9_name
	# Figures[[reg]][['Figure9']][["dim"]] <- fig9_dim
	# dev.off()
	#
	#
	# # ---- Figure 10 ----
	# # ---- Time Series of Presence Parameters (alpha) ----
	# fig10_name <- paste0("alpha1_spp_timeSeries", reg, ".png")
	# fig10_dim <- c(width=sP_w, height=sP_h)
	#
	# dev.new(width=fig10_dim[1], height=fig10_dim[2])
	# par(mfrow=sppPar_mfr, mar=c(1,1,1,0.1), oma=c(0.1,0.1,1.5,0.1), ps=6, mgp=c(0.6,0.1,0), tcl=-0.1, cex=1)
	#
	# for(s in 1:n_spp){
	# 	t_s <- u_spp[s]
	# 	t_dat <- ab[par=="alpha" & ab_ind==1 & (spp)==t_s]
	# 	plot_ab(X=t_dat, t_spp=t_s)
	# }
	# mtext(paste0(reg, " Alpha[1] (presence intercept)"), outer=TRUE, side=3, font=2, line=0.25, cex=1.25)
	#
	# Figures[[reg]][['Figure10']][["figure"]] <- recordPlot()
	# Figures[[reg]][['Figure10']][["name"]] <- fig10_name
	# Figures[[reg]][['Figure10']][["dim"]] <- fig10_dim
	# dev.off()
	
	
	# ---- Figure 11:  ----

	graphics.off()
}

# ==========================
# = Plot Figures (or Save) =
# ==========================
plot_Figures <- function(x, FUN, ...){
	fig_fun <- match.fun(FUN)
	d <- x$dim

	if(FUN!="dev.new"){
		fig_fun(width=d[1], height=d[2], file=gsub("\\.png", paste0(".",FUN), x$name), ...)
		replayPlot(x$figure)
		dev.off()
	}else{
		fig_fun(width=d[1], height=d[2], ...)
		replayPlot(x$figure)
	}
}



for(i in 1:length(p)){
	plot_Figures(Figures[[i]][[1]], "dev.new")
}




