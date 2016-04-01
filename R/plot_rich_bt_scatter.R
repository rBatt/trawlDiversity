#' Scatter Plot of Richness and Bottom Temperature
#' 
#' Plot annual time series of richness vs annual mean of temperature
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
#' There are several functions that are run through the process_msom_figures script. Richness and temperature plots are \code{\link{plot_btemp_map}}, \code{\link{plot_rich_bt_scatter}}, and \code{\link{plot_rich_bt_ts}}. Figures for colonization, extinction, and the species and places associated with those processes are \code{\link{plot_ce_wrap}}, \code{\link{plot_col_vs_unobs}}, \code{\link{plot_colExt_perStrat}}, and \code{\link{plot_rank_temp}}. Figures for diagnostics are \code{\link{plot_traceplot}} and \code{\link{plot_post_corr}}.
#' 
#' @export
plot_rich_bt_scatter <- function(prn, Figures){
	unpack_p(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
	fig2_name <- paste0("richness_bt_scatter_", reg, ".png")
	fig2_dim <- c(3.5, 5)

	dev.new(width=fig2_dim[1], height=fig2_dim[2])
	par(mfrow=c(2,1), mar=c(1.75,1.5,0.25,0.25), oma=c(0.1,0.1,0.75,0.1), mgp=c(0.75,0.1,0), tcl=-0.1, ps=8, cex=1)
	
	plot(processed[,list(bt_ann,naive_rich)], type="p", ylab="Naive Region Richness", xlab="Annual Mean Bottom Temperature")
	plot(processed[,list(bt_ann,reg_rich)], type="p", ylab="MSOM Region Richness", xlab="Annual Mean Bottom Temperature")
	mtext(reg, side=3, line=0, outer=TRUE, font=2)
	
	Figures[[reg]][['Figure2']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure2']][["name"]] <- fig2_name
	Figures[[reg]][['Figure2']][["dim"]] <- fig2_dim
	
	return(Figures)
}