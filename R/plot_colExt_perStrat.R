#' Number of Colonizations and Extinctions in Each Stratum
#' 
#' Plots the number of colonizations, extinctions, and their difference as the time series total in each stratum
#' 
#' @param prn the p object (processed msom; output from \code{process_msomStatic})
#' @param Figures option list to which the figure and its information should be added
#' 
#' @details
#' Each of the 3 types of plots has two panels: The first panel is a map of the strata indicating the the values seen in the data. Each stratum is a point. The second panel is a thin plat spline of the first panel (see \code{fields::Tps}), which acts to smooth the values in the convex hull of the strata.
#' 
#' The units of the map are a bit complicated. They generally indicate the mean number of colonizations (or extinctions, or C-E) per year for each stratum. However, the number of 'colonizations' for a stratum in a given year is the number of species that colonized the region at that stratum, divided by the total number of strata that those same species colonized that year. Thus, species that colonize many strata do not contribute much to the colonization score of any single stratum, and the sum of each stratum's colonization score for a year is the number of species that colonized the whole region. It is this colonization score that is averaged over years (which can be different from the sum across years b/c some strata aren't sampled every year).
#' 
#' The function uses \code{unpack_p} to get much of the information it needs. The \code{Figures} object is a list whose first level is intended to be the region. The second level is specific to each figure. The third level has 3 elements: 'figure', 'name', and 'dim'. The 'figure' element is the result of using \code{\link{recordPlot}} on what is plotted. The name of the figure is, e.g., what the saved figure would be called. The dim is the width and height (in that order) in inches.
#' 
#' @return
#' Returns the Figure object
#' 
#' @seealso 
#' There are several functions that are run through the process_msom_figures script. Richness and temperature plots are \code{\link{plot_btemp_map}}, \code{\link{plot_rich_bt_scatter}}, and \code{\link{plot_rich_bt_ts}}. Figures for colonization, extinction, and the species and places associated with those processes are \code{\link{plot_ce_wrap}}, \code{\link{plot_col_vs_unobs}}, \code{\link{plot_colExt_perStrat}}, and \code{\link{plot_rank_temp}}. Figures for diagnostics are \code{\link{plot_traceplot}} and \code{\link{plot_post_corr}}.
#' 
#' @export
plot_colExt_perStrat <- function(prn, Figures){
	requireNamespace("fields", quietly=TRUE)
	unpack_p(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
	r_lon <- colonization$n_spp_col_weighted_tot[,range(lon)]
	r_lat <- colonization$n_spp_col_weighted_tot[,range(lat)]
	
	lay_logic <- diff(r_lon) > diff(r_lat)
	if(lay_logic){
		fig7_mfr <- c(3,2)
		fig7_h <- 3
		fig7_w <- (diff(r_lon) * fig7_h / diff(r_lat)) * (fig7_mfr[2]/fig7_mfr[1])
	}else{
		# fig7_mfr <- c(1*3,2)
# 		fig7_w <- 2.5#*3
# 		fig7_h <- (diff(r_lat) * fig7_w / diff(r_lon)) / (fig7_mfr[2]/fig7_mfr[1])
		fig7_mfr <- c(2,1*3)
		fig7_h <- 7#*3
		fig7_w <- (diff(r_lon) * fig7_h / diff(r_lat)) * (fig7_mfr[2]/fig7_mfr[1])
		
		# lay_logic <- !lay_logic
	}

	fig7_name <- paste0("colonizations_per_stratum_", reg, ".png")
	fig7_dim <- c(fig7_w, fig7_h)
	dev.new(width=fig7_w, height=fig7_h)
	if(lay_logic){
		par(mfrow=fig7_mfr)
	}else{
		par(mfcol=fig7_mfr)
	}
	par(mar=c(1.25,1.25,0.1,0.1), oma=c(0.1,1,1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=8, cex=1)

	# site-specific colonizations from data
	colonization$n_spp_col_weighted_tot[,plot_space(lon, lat, n_spp_col_weighted, TRUE, pch=19, ylab="", xlab="")]
	map(add=TRUE, fill=TRUE, col="white")
	mtext("Colonizations (C)", side=ifelse(lay_logic, 2, 3), line=ifelse(lay_logic, 1, 0.25))
	# smoothed map for convex hull of site observations
	colonization$n_spp_col_weighted_tot[,plot_space(lon, lat, n_spp_col_weighted, pch=19, ylab="", xlab="")]
	map(add=TRUE, fill=TRUE, col="white")
	
	
	# site-specific extinctions from data
	colonization$n_spp_ext_weighted_tot[,plot_space(lon, lat, n_spp_ext_weighted, TRUE, pch=19, ylab="", xlab="")]
	map(add=TRUE, fill=TRUE, col="white")
	mtext("Extinctions (E)", side=ifelse(lay_logic, 2, 3), line=ifelse(lay_logic, 1, 0.25))
	# smoothed map for convex hull of site observations
	colonization$n_spp_ext_weighted_tot[,plot_space(lon, lat, n_spp_ext_weighted, pch=19, ylab="", xlab="")]
	map(add=TRUE, fill=TRUE, col="white")
	
	# site-specific extinctions from data
	cre <- merge(colonization$n_spp_col_weighted_tot, colonization$n_spp_ext_weighted_tot, by=c("stratum","lon","lat","depth"), all=TRUE)
	# cre <- data.table(cre[,list(stratum, n_spp_ext_weighted, n_spp_col_weighted)], cre[,strat2lld(stratum)])
# 	cre[is.na(n_spp_ext_weighted), n_spp_ext_weighted:=0]
# 	cre[is.na(n_spp_col_weighted), n_spp_col_weighted:=0]
	cre[,rel_col_ext:=(n_spp_col_weighted - n_spp_ext_weighted) ]
	cre[,plot_space(lon, lat, rel_col_ext, TRUE, pch=19, ylab="", xlab="")]
	map(add=TRUE, fill=TRUE, col="white")
	mtext("C - E", side=ifelse(lay_logic, 2, 3), line=ifelse(lay_logic, 1, 0.25))
	# smoothed map for convex hull of site observations
	cre[,plot_space(lon, lat, rel_col_ext, pch=19, ylab="", xlab="")]
	map(add=TRUE, fill=TRUE, col="white")
	
	

	Figures[[reg]][['Figure7']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure7']][["name"]] <- fig7_name
	Figures[[reg]][['Figure7']][["dim"]] <- fig7_dim

	return(Figures)
	
	# ---- Figure 8 ----
	# ---- Number of extinctions per stratum ----
	# fig8_name <- paste0("extinctions_per_stratum_", reg, ".png")
# 	fig8_dim <- c(fig7_w, fig7_h)
# 	dev.new(width=fig7_w, height=fig7_h)
# 	par(mfrow=fig7_mfr, mar=c(1.5,1.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=8, cex=1)

	# # site-specific extinctions from data
# 	colonization$n_spp_ext_weighted_tot[,plot_space(lon, lat, n_spp_ext_weighted, TRUE, pch=19)]
# 	map(add=TRUE, fill=TRUE, col="white")
# 	# smoothed map for convex hull of site observations
# 	colonization$n_spp_ext_weighted_tot[,plot_space(lon, lat, n_spp_ext_weighted, pch=19)]
# 	map(add=TRUE, fill=TRUE, col="white")

	# Figures[[reg]][['Figure8']][["figure"]] <- recordPlot()
	# Figures[[reg]][['Figure8']][["name"]] <- fig8_name
	# Figures[[reg]][['Figure8']][["dim"]] <- fig8_dim
	# dev.off()
	
	
	# ---- Figure 8.5 ----
	# ---- Colonizations relative to extinctions per stratum ----
	# fig8.5_name <- paste0("ColRelExt_per_stratum_", reg, ".png")
	# fig8.5_dim <- c(fig7_w, fig7_h)
	# dev.new(width=fig7_w, height=fig7_h)
	# par(mfrow=fig7_mfr, mar=c(1.5,1.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=8, cex=1)
	#
	# strat2lld <- function(x){
	# 	s <- strsplit(x, split=" ")
	# 	lon <- sapply(s, function(x)x[1])
	# 	lat <- sapply(s, function(x)x[2])
	# 	depth_interval <- sapply(s, function(x)x[3])
	# 	data.table(lon=as.numeric(lon), lat=as.numeric(lat), depth_interval=as.numeric(depth_interval))
	# } 

	# # site-specific extinctions from data
# 	cre <- merge(colonization$n_spp_col_weighted_tot, colonization$n_spp_ext_weighted_tot, by=c("stratum","lon","lat","depth"), all=TRUE)
# 	# cre <- data.table(cre[,list(stratum, n_spp_ext_weighted, n_spp_col_weighted)], cre[,strat2lld(stratum)])
# # 	cre[is.na(n_spp_ext_weighted), n_spp_ext_weighted:=0]
# # 	cre[is.na(n_spp_col_weighted), n_spp_col_weighted:=0]
# 	cre[,rel_col_ext:=(n_spp_col_weighted - n_spp_ext_weighted) ]
# 	cre[,plot_space(lon, lat, rel_col_ext, TRUE, pch=19)]
# 	map(add=TRUE, fill=TRUE, col="white")
# 	# smoothed map for convex hull of site observations
# 	cre[,plot_space(lon, lat, rel_col_ext, pch=19)]
# 	map(add=TRUE, fill=TRUE, col="white")
}