
library("rstan")
library("trawlDiversity")
library("rbLib")
library("R2jags")
library("maps")

setwd("~/Documents/School&Work/pinskyPost/trawl")
load("trawlDiversity/pkgBuild/results/processedMsom/p.RData")

Figures <- list()

# ---- Function to plot time series of species parameter ----
plot_ab <- function(X, t_spp){

	si <- sppImg(t_spp)
	if(!is.null(si)){
		par(new=T)
	}
	
	fin <- par("fin")[2]
	fac <- 0.01*fin^3/2

	plot(X[,year], X[,value], col=adjustcolor('gray', fac), cex=0.5, pch=21, bg=adjustcolor('white',fac))
	mu <- X[,list(mu=mean(value)),by="year"]
	mu[,lines(year, mu, lwd=2, col='gray')]
	mu[,lines(year, mu, lwd=1, col='white')]
	if(is.null(si)){
		common_name <- spp.key[spp==t_spp, una(common)]
		mtext(paste(t_spp, common_name, sep="\n"), side=3)
	}

}


for(reg_num in 1:length(p)){
	
	# ========
	# = Prep =
	# ========
	rd <- p[[reg_num]]$rd
	processed <- p[[reg_num]]$processed
	bt <- p[[reg_num]]$bt
	colonization <- p[[reg_num]]$colonization
	param_iters <- p[[reg_num]]$param_iters
	ab <- p[[reg_num]]$ab

	reg <- processed[,una(reg)]
	lang = "JAGS"

	if(lang == "Stan"){
		pars_trace <- c("Omega","alpha_mu[1]", "alpha_mu[2]", "alpha_mu[3]", "alpha_mu[4]", "alpha_mu[5]", "beta_mu[1]")
	}else{
		pars_trace <- c("Omega","alpha_mu[1]", "alpha_mu[2]", "alpha_mu[3]", "alpha_mu[4]", "alpha_mu[5]", "beta_mu")
	}

	naive_rich <- processed[,naive_rich, by='year']
	reg_rich <- processed[,reg_rich, by='year']
	bt_ann <- bt[,list(bt_ann=mean(bt)), by='year']

	n_pars <- length(pars_trace)
	n_yrs <- param_iters[,lu(year)]
	n_spp <- rd[,lu(spp)]


	# ===========
	# = Figures =
	# ===========

	# ---- Figure 1 ----
	fig1_name <- paste0("richness_bt_timeSeries_", reg, ".png")
	fig1_dim <- c(3.5, 6)
	
	dev.new(width=3.5, height=6)
	par(mfrow=c(3,1), mar=c(1.75,1.5,0.25,0.25), oma=c(0.1,0.1,0.75,0.1), mgp=c(0.75,0.1,0), tcl=-0.1, ps=8, cex=1)
	
	plot(naive_rich, type="o", ylab="Naive Region Richness", xlab="Year")
	plot(reg_rich, type="o", xlab="Year", ylab="MSOM Region Richness")
	plot(bt_ann, type="o", xlab="Year", ylab="Annual Mean Bottom Temperature")
	mtext(reg, side=3, line=0, outer=TRUE, font=2)
	
	Figures[[reg]][['Figure1']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure1']][["name"]] <- fig1_name
	Figures[[reg]][['Figure1']][["dim"]] <- fig1_dim
	dev.off()

	# ---- Figure 2 ----
	fig2_name <- paste0("richness_bt_scatter_", reg, ".png")
	fig2_dim <- c(3.5, 5)

	dev.new(width=3.5, height=5)
	par(mfrow=c(2,1), mar=c(1.75,1.5,0.25,0.25), oma=c(0.1,0.1,0.75,0.1), mgp=c(0.75,0.1,0), tcl=-0.1, ps=8, cex=1)
	
	plot(processed[,list(bt_ann,naive_rich)], type="p", ylab="Naive Region Richness", xlab="Annual Mean Bottom Temperature")
	plot(processed[,list(bt_ann,reg_rich)], type="p", ylab="MSOM Region Richness", xlab="Annual Mean Bottom Temperature")
	mtext(reg, side=3, line=0, outer=TRUE, font=2)
	
	Figures[[reg]][['Figure2']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure2']][["name"]] <- fig2_name
	Figures[[reg]][['Figure2']][["dim"]] <- fig2_dim
	dev.off()

	# ---- Figure 3 ----
	fig3_name <- paste0("btempMap_", reg, ".png")
	f3_mfrow <- auto.mfrow(n_yrs)
	f3_height <- 6
	f3_width <- f3_mfrow[2]*f3_height/f3_mfrow[1]
	fig3_dim <- c(f3_width, f3_height)
	
	dev.new(width=f3_width, height=f3_height)
	par(mfrow=f3_mfrow, oma=c(0.1,0.1, 1,0.1), mar=c(1,1,0.1,0.1), mgp=c(0.75,0.1,0), tcl=-0.15, cex=1, ps=8)
	
	bt[,j={
		plot(lon, lat, type="n")
		map(add=TRUE)
		points(lon, lat, col=bt_col, pch=20)
		mtext(unique(year), side=3, adj=0.1, line=-0.75, font=2)
	}, by="year"]
	mtext(paste(reg, "Bottom Temperature"), outer=TRUE, side=3, line=0.25, font=2)
	
	Figures[[reg]][['Figure3']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure3']][["name"]] <- fig3_name
	Figures[[reg]][['Figure3']][["dim"]] <- fig3_dim
	dev.off()
	
	# ---- Figure 4 ----
	fig4_name <- paste0("traceplot_", reg, ".png")
	f4_mfrow <- c(n_pars, n_yrs)
	f4_height <- 5
	f4_width <- f4_mfrow[2]*f4_height/f4_mfrow[1]
	fig4_dim <- c(f4_width, f4_height)

	dev.new(width=f4_width, height=f4_height)
	par(mfrow=f4_mfrow, oma=c(1,1, 1,0.1), mar=c(0.5,0.5,0.1,0.1), mgp=c(0.25,0.1,0), tcl=-0.1, cex=1, ps=6)
	
	for(h in 1:length(pars_trace)){
		for(i in 1:n_yrs){
			t_yr <- param_iters[,una(year)][i]
			t_iters <- param_iters[year==t_yr]
			mytrace(t_iters, pars=pars_trace[h], lang=lang, xaxt='n')
			if(i == 1){
				mtext(pars_trace[h], side=2, line=0.75)
			}
		}
	}
	
	Figures[[reg]][['Figure4']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure4']][["name"]] <- fig4_name
	Figures[[reg]][['Figure4']][["dim"]] <- fig4_dim
	dev.off()


	# ---- Figure 5 ----
	fig5_name <- paste0("posteriorCorrelation_", reg, ".png")
	fig5_dim <- c(7, 7)
	
	dev.new(fig5_dim[1], fig5_dim[2])
	
	pairs(param_iters[year==param_iters[,una(year)][1], eval(s2c(pars_trace))])
	
	Figures[[reg]][['Figure5']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure5']][["name"]] <- fig5_name
	Figures[[reg]][['Figure5']][["dim"]] <- fig5_dim
	dev.off()


	# ---- Figure 6 ----
	fig6_name <- paste0("Colonization_UnobsSpp_", reg, ".png")
	fig6_dim <- c(3.5, 3.5)
	
	dev.new(fig6_dim[1], fig6_dim[2])
	
	processed[,plot(unobs_rich[-length(unobs_rich)], n_col[-1], xlab="Unobserved species present last year", ylab="Species colonizing this year")]
	abline(a=0, b=1)
	
	Figures[[reg]][['Figure6']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure6']][["name"]] <- fig6_name
	Figures[[reg]][['Figure6']][["dim"]] <- fig6_dim
	dev.off()


	# ---- Figure 7 ----
	# ---- Number of Colonizations per Stratum ----
	r_lon <- colonization$n_spp_col_weighted_tot[,range(lon)]
	r_lat <- colonization$n_spp_col_weighted_tot[,range(lat)]
	if(diff(r_lon) > diff(r_lat)){
		fig7_mfr <- c(2,1)
		fig7_h <- 2.5
		fig7_w <- (diff(r_lon) * fig7_h / diff(r_lat)) * (fig7_mfr[2]/fig7_mfr[1])
	}else{
		fig7_mfr <- c(1,2)
		fig7_w <- 2.5
		fig7_h <- (diff(r_lat) * fig7_w / diff(r_lon)) / (fig7_mfr[2]/fig7_mfr[1])
	}

	fig7_name <- paste0("colonizations_per_stratum_", reg, ".png")
	fig7_dim <- c(fig7_w, fig7_h)
	dev.new(width=fig7_w, height=fig7_h)
	par(mfrow=fig7_mfr, mar=c(1.5,1.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=8, cex=1)

	# site-specific colonizations from data
	colonization$n_spp_col_weighted_tot[,plot_space(lon, lat, n_spp_col_weighted, TRUE, pch=19)]
	map(add=TRUE, fill=TRUE, col="white")
	# smoothed map for convex hull of site observations
	colonization$n_spp_col_weighted_tot[,plot_space(lon, lat, n_spp_col_weighted, pch=19)]
	map(add=TRUE, fill=TRUE, col="white")

	Figures[[reg]][['Figure7']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure7']][["name"]] <- fig7_name
	Figures[[reg]][['Figure7']][["dim"]] <- fig7_dim
	dev.off()
	
	
	# ---- Figure 8 ----
	# ---- Number of extinctions per stratum ----
	fig8_name <- paste0("extinctions_per_stratum_", reg, ".png")
	fig8_dim <- c(fig7_w, fig7_h)
	dev.new(width=fig7_w, height=fig7_h)
	par(mfrow=fig7_mfr, mar=c(1.5,1.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=8, cex=1)

	# site-specific extinctions from data
	colonization$n_spp_ext_weighted_tot[,plot_space(lon, lat, n_spp_ext_weighted, TRUE, pch=19)]
	map(add=TRUE, fill=TRUE, col="white")
	# smoothed map for convex hull of site observations
	colonization$n_spp_ext_weighted_tot[,plot_space(lon, lat, n_spp_ext_weighted, pch=19)]
	map(add=TRUE, fill=TRUE, col="white")

	Figures[[reg]][['Figure8']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure8']][["name"]] <- fig8_name
	Figures[[reg]][['Figure8']][["dim"]] <- fig8_dim
	dev.off()
	
	
	# ---- Figure 8.5 ----
	# ---- Colonizations relative to extinctions per stratum ----
	fig8.5_name <- paste0("ColRelExt_per_stratum_", reg, ".png")
	fig8.5_dim <- c(fig7_w, fig7_h)
	dev.new(width=fig7_w, height=fig7_h)
	par(mfrow=fig7_mfr, mar=c(1.5,1.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=8, cex=1)
	
	strat2lld <- function(x){
		s <- strsplit(x, split=" ")
		lon <- sapply(s, function(x)x[1])
		lat <- sapply(s, function(x)x[2])
		depth_interval <- sapply(s, function(x)x[3])
		data.table(lon=as.numeric(lon), lat=as.numeric(lat), depth_interval=as.numeric(depth_interval))
	} 

	# site-specific extinctions from data
	cre <- merge(colonization$n_spp_col_weighted_tot, colonization$n_spp_ext_weighted_tot, by=c("stratum","lon","lat","depth"), all=TRUE)
	# cre <- data.table(cre[,list(stratum, n_spp_ext_weighted, n_spp_col_weighted)], cre[,strat2lld(stratum)])
# 	cre[is.na(n_spp_ext_weighted), n_spp_ext_weighted:=0]
# 	cre[is.na(n_spp_col_weighted), n_spp_col_weighted:=0]
	cre[,rel_col_ext:=(n_spp_col_weighted - n_spp_ext_weighted) ]
	cre[,plot_space(lon, lat, rel_col_ext, TRUE, pch=19)]
	map(add=TRUE, fill=TRUE, col="white")
	# smoothed map for convex hull of site observations
	cre[,plot_space(lon, lat, rel_col_ext, pch=19)]
	map(add=TRUE, fill=TRUE, col="white")
	


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

plot_Figures(Figures[[1]][[1]], "dev.new")




