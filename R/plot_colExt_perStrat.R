#' Number of Colonizations and Extinctions in Each Stratum
#' 
#' Plots the number of colonizations, extinctions, and their difference as the time series total in each stratum
#' 
#' @param prn the p object (processed msom; output from \code{process_msomStatic})
#' @param Figures option list to which the figure and its information should be added
#' @param FUN graphical device function
#' @param ... arguments passed to plot_device
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
#' There are several functions that are run through the process_msom_figures script. Richness and temperature plots are \code{\link{plot_btemp_map}}, \code{\link{plot_rich_bt_scatter}}, and \code{\link{plot_rich_bt_ts}}. Figures for colonization, extinction, and the species and places associated with those processes are \code{\link{plot_ce_wrap}}, \code{\link{plot_col_vs_unobsSpp}}, \code{\link{plot_colExt_perStrat}}, and \code{\link{plot_rank_temp}}. Figures for diagnostics are \code{\link{plot_traceplot}} and \code{\link{plot_post_corr}}.
#' 
#' @export
plot_colExt_perStrat <- function(prn, Figures, FUN="dev.new", ...){
	requireNamespace("maps", quietly=TRUE)
	pup <- unpack_p(prn)
mapply(x=names(pup), value=pup, function(x, value){assign(x, value, envir=parent.frame(n=2));invisible(NULL)})(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
	fig_num <- "Figure7"
	
	r_lon <- colonization$n_spp_col_weighted_tot[,range(lon)]
	r_lat <- colonization$n_spp_col_weighted_tot[,range(lat)]
	
	lay_logic <- diff(r_lon) > diff(r_lat)
	if(lay_logic){
		fig7_mfr <- c(3,2)
		fig7_h <- 3
		fig7_w <- (diff(r_lon) * fig7_h / diff(r_lat)) * (fig7_mfr[2]/fig7_mfr[1])
	}else{
		fig7_mfr <- c(2,1*3)
		fig7_h <- 7
		fig7_w <- (diff(r_lon) * fig7_h / diff(r_lat)) * (fig7_mfr[2]/fig7_mfr[1])
	}

	fig7_name <- paste0("colonizations_per_stratum_", reg, ".png")
	fig7_dim <- c(fig7_w, fig7_h)
	
	Figures[[reg]][[fig_num]][["figure"]] <- list()
	Figures[[reg]][[fig_num]][["name"]] <- fig7_name
	Figures[[reg]][[fig_num]][["dim"]] <- fig7_dim
	
	Figures[[reg]][[fig_num]] <- plot_device(Figures[[reg]][[fig_num]], FUN, ...)
	
	if(lay_logic){
		par(mfrow=fig7_mfr)
	}else{
		par(mfcol=fig7_mfr)
	}
	par(mar=c(1.25,1.25,0.1,0.1), oma=c(0.1,1,1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=8, cex=1)

	# site-specific colonizations from data
	colonization$n_spp_col_weighted_tot[,plot_space(lon, lat, n_spp_col_weighted, TRUE, pch=19, ylab="", xlab="")]
	maps::map(add=TRUE, fill=TRUE, col="white")
	mtext("Colonizations (C)", side=ifelse(lay_logic, 2, 3), line=ifelse(lay_logic, 1, 0.25))
	# smoothed map for convex hull of site observations
	colonization$n_spp_col_weighted_tot[,plot_space(lon, lat, n_spp_col_weighted, pch=19, ylab="", xlab="")]
	maps::map(add=TRUE, fill=TRUE, col="white")
	
	
	# site-specific extinctions from data
	colonization$n_spp_ext_weighted_tot[,plot_space(lon, lat, n_spp_ext_weighted, TRUE, pch=19, ylab="", xlab="")]
	maps::map(add=TRUE, fill=TRUE, col="white")
	mtext("Extinctions (E)", side=ifelse(lay_logic, 2, 3), line=ifelse(lay_logic, 1, 0.25))
	# smoothed map for convex hull of site observations
	colonization$n_spp_ext_weighted_tot[,plot_space(lon, lat, n_spp_ext_weighted, pch=19, ylab="", xlab="")]
	maps::map(add=TRUE, fill=TRUE, col="white")
	
	
	# site-specific extinctions from data
	cre <- merge(colonization$n_spp_col_weighted_tot, colonization$n_spp_ext_weighted_tot, by=c("stratum","lon","lat","depth"), all=TRUE)
	cre[,rel_col_ext:=(n_spp_col_weighted - n_spp_ext_weighted) ]
	cre[,plot_space(lon, lat, rel_col_ext, TRUE, pch=19, ylab="", xlab="")]
	maps::map(add=TRUE, fill=TRUE, col="white")
	mtext("C - E", side=ifelse(lay_logic, 2, 3), line=ifelse(lay_logic, 1, 0.25))
	# smoothed map for convex hull of site observations
	cre[,plot_space(lon, lat, rel_col_ext, pch=19, ylab="", xlab="")]
	maps::map(add=TRUE, fill=TRUE, col="white")
	
	
	
	if(is.null(Figures[[reg]][[fig_num]][["fig_loc"]])){
		Figures[[reg]][[fig_num]][["figure"]] <- recordPlot()
	}else{
		Figures[[reg]][[fig_num]][["figure"]] <- NULL
		dev.off()
	}
	
	return(Figures)
	
}