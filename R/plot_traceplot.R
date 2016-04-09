#' Plot a Traceplot for MSOM Hyperparameters
#' 
#' Plots the traceplot for all years for the hyperparameters from MSOM model output
#' 
#' @param prn the p object (processed msom; output from \code{process_msomStatic})
#' @param Figures option list to which the figure and its information should be added
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
plot_traceplot <- function(prn, Figures){
	unpack_p(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
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
	
	return(Figures)
}