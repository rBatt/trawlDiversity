#' Plot Posterior Correlation
#' 
#' Plots the posterior correlation of hyperparameters from an MSOM output
#' 
#' @param prn the p object (processed msom; output from \code{process_msomStatic})
#' @param Figures option list to which the figure and its information should be added
#' @param yr integer indicating the index for the year to plot. Note that this is not the calendar year
#' 
#' @details
#' The function uses \code{unpack_p} to get much of the information it needs. The \code{Figures} object is a list whose first level is intended to be the region. The second level is specific to each figure. The third level has 3 elements: 'figure', 'name', and 'dim'. The 'figure' element is the result of using \code{\link{recordPlot}} on what is plotted. The name of the figure is, e.g., what the saved figure would be called. The dim is the width and height (in that order) in inches.
#' 
#' @return
#' Returns the Figure object
#' 
#' @seealso 
#' There are several functions that are run through the process_msom_figures script. Richness and temperature plots are \code{\link{plot_btemp_map}}, \code{\link{plot_rich_bt_scatter}}, and \code{\link{plot_rich_bt_ts}}. Figures for colonization, extinction, and the species and places associated with those processes are \code{\link{plot_ce_wrap}}, \code{\link{plot_col_vs_unobsSpp}}, \code{\link{plot_colExt_perStrat}}, and \code{\link{plot_rank_temp}}. Figures for diagnostics are \code{\link{plot_traceplot}} and \code{\link{plot_post_corr}}.
#' 
#' @export
plot_post_corr <- function(prn, Figures, yr=1){
	unpack_p(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
	fig5_name <- paste0("posteriorCorrelation_", reg, ".png")
	fig5_dim <- c(7, 7)
	
	dev.new(fig5_dim[1], fig5_dim[2])
	
	pairs(param_iters[year==param_iters[,una(year)][yr], eval(s2c(pars_trace))])
	
	Figures[[reg]][['Figure5']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure5']][["name"]] <- fig5_name
	Figures[[reg]][['Figure5']][["dim"]] <- fig5_dim
	
	return(Figures)
}