#' Plot Posterior Correlation
#' 
#' Plots the posterior correlation of hyperparameters from an MSOM output
#' 
#' @param prn the p object (processed msom; output from \code{process_msomStatic})
#' @param Figures option list to which the figure and its information should be added
#' @param yr integer indicating the index for the year to plot. Note that this is not the calendar year
#' @param FUN graphical device function
#' @param ... arguments passed to plot_device
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
plot_post_corr <- function(prn, Figures, yr=1, FUN="dev.new", ...){
	pup <- unpack_p(prn)
mapply(x=names(pup), value=pup, function(x, value){assign(x, value, envir=parent.frame(n=2));invisible(NULL)})(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
	fig_num <- "Figure5"
	
	fig5_name <- paste0("posteriorCorrelation_", reg, ".png")
	fig5_dim <- c(7, 7)
	
	Figures[[reg]][[fig_num]][["figure"]] <- list()
	Figures[[reg]][[fig_num]][["name"]] <- fig5_name
	Figures[[reg]][[fig_num]][["dim"]] <- fig5_dim
	
	Figures[[reg]][[fig_num]] <- plot_device(Figures[[reg]][[fig_num]], FUN, ...)
	
	pairs(param_iters[year==param_iters[,una(year)][yr], eval(s2c(pars_trace))])
	
	if(is.null(Figures[[reg]][[fig_num]][["fig_loc"]])){
		Figures[[reg]][[fig_num]][["figure"]] <- recordPlot()
	}else{
		Figures[[reg]][[fig_num]][["figure"]] <- NULL
		dev.off()
	}
	
	return(Figures)
}