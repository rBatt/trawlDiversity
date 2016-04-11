
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


for(reg_num in 1:length(p)){
	

	t_prn <- p[[reg_num]]

	# ===========
	# = Figures =
	# ===========

	# ---- Figure 1 ----
	Figures <- plot_rich_bt_ts(t_prn, Figures, FUN="pdf")

	# ---- Figure 2 ----
	Figures <- plot_rich_bt_scatter(t_prn, Figures, FUN="pdf")

	# ---- Figure 3 ----
	Figures <- plot_btemp_map(t_prn, Figures, FUN="pdf")
	
	# ---- Figure 4 ----
	Figures <- plot_traceplot(t_prn, Figures, FUN="pdf")

	# ---- Figure 5 ----
	Figures <- plot_post_corr(t_prn, Figures, yr=1, FUN="pdf")

	# ---- Figure 6 ----
	# Figures <- plot_col_vs_unobsSpp(t_prn, Figures, FUN="pdf")

	# ---- Figure 7 ----
	# ---- Number of Colonizations per Stratum ----	
	Figures <- plot_colExt_perStrat(t_prn, Figures, FUN="pdf")

	# ---- Figure 8 ----
	# ---- Plot Information and Identity of Colonizers, Leavers, etc ----
	Figures <- plot_ce_wrap(t_prn, Figures, spp_cat="col", width.max=12, height.max=18, max_spp_columns=12, FUN="pdf")
	Figures <- plot_ce_wrap(t_prn, Figures, spp_cat="ext", width.max=12, height.max=18, max_spp_columns=12, FUN="pdf")
	Figures <- plot_ce_wrap(t_prn, Figures, spp_cat="both", width.max=12, height.max=18, max_spp_columns=12, FUN="pdf")
	
	# ---- Figure 9:  ----
	# ---- Plot Number of Colonizers, Leavers, and Plot Temp Rank ----
	Figures <- plot_rank_temp(t_prn, Figures, FUN="pdf")
	
}

# ==========================
# = Plot Figures (or Save) =
# ==========================
plot_Figures <- function(x, FUN, ...){
	fig_fun <- match.fun(FUN)
	d <- x$dim

	if(FUN!="dev.new"){
		fig_name <- gsub("\\.png", paste0(".",FUN), x$name)
		fig_fun(width=d[1], height=d[2], file=fig_name, ...)
		replayPlot(x$figure)
		dev_val <- dev.off()
	}else{
		fig_fun(width=d[1], height=d[2], ...)
		replayPlot(x$figure)
		fig_name <- NULL
	}
	return(fig_name)
}


save(Figures, "trawlDiversity/pkgBuild/figures/Figures.RData")

#
# for(i in 1:length(p)){
# 	plot_Figures(Figures[[i]][[1]], "dev.new")
# }


for(i in 1:length(p)){
	od <- getwd()
	
	fig_names <- c()
	td <- c()
	for(f in 1:length(Figures[[i]])){
		f_l <- Figures[[reg]][[f]][["fig_loc"]]
		if(!is.null(f_l)){
			fig_names[f] <- f_l
			td[f] <- od
		}else{
			td[f] <- tempdir()
			setwd(td[f])
			fig_names[f] <- plot_Figures(Figures[[i]][[f]], "pdf")
			setwd(od)
		}
	}
	ins <- file.path(td, fig_names) 
	i_reg <- p[[i]]$rd[,una(reg)]
	save_name <- paste0(od, "/trawlDiversity/pkgBuild/figures/processed_msom_figures_", i_reg, ".pdf")
	combine_pdf(shQuote(ins, "cmd"), shQuote(save_name, "cmd"))
	
	setwd(od)
}



