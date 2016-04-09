#' Plot colonizers and unobserved species
#' 
#' Colonizers from data this year versus MSOM estimates of last year's unobserved species
#' 
#' @param prn the p object (processed msom; output from \code{process_msomStatic})
#' @param Figures option list to which the figure and its information should be added
#' 
#' @details
#' The function uses \code{unpack_p} to get much of the information it needs. The \code{Figures} object is a list whose first level is intended to be the region. The second level is specific to each figure. The third level has 3 elements: 'figure', 'name', and 'dim'. The 'figure' element is the result of using \code{\link{recordPlot}} on what is plotted. The name of the figure is, e.g., what the saved figure would be called. The dim is the width and height (in that order) in inches.
#' 
#' I'm not sure how valuable this figure is. But the original idea was that if the survey didn't detect a species last year (but it was present), maybe the MSOM would have estimated it as a present but unobserved species. And if that unobserved species was observed this year, it could look like a colonization. In that case, observed colonizations might just be correlated with the number of unobserved species from last year. Obviously unobserved species from last year could have gone extinct, and, more importantly, the MSOM isn't perfect. So I'm hesitant to read too far into the results of this figure. However, I do think it is useful to note the relative magnitude of unobserved species in a given year, and the number of colonizers in a given year.
#' 
#' @return
#' Returns the Figure object
#' 
#' @seealso 
#' There are several functions that are run through the process_msom_figures script. Richness and temperature plots are \code{\link{plot_btemp_map}}, \code{\link{plot_rich_bt_scatter}}, and \code{\link{plot_rich_bt_ts}}. Figures for colonization, extinction, and the species and places associated with those processes are \code{\link{plot_ce_wrap}}, \code{\link{plot_col_vs_unobsSpp}}, \code{\link{plot_colExt_perStrat}}, and \code{\link{plot_rank_temp}}. Figures for diagnostics are \code{\link{plot_traceplot}} and \code{\link{plot_post_corr}}.
#' 
#' @export

plot_col_vs_unobsSpp <- function(prn, Figures){
	unpack_p(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
	fig6_name <- paste0("Colonization_UnobsSpp_", reg, ".png")
	fig6_dim <- c(3.5, 3.5)
	
	dev.new(fig6_dim[1], fig6_dim[2])
	par(mar=c(1.75,1.75,0.2,0.2), cex=1, ps=8, mgp=c(0.85,0.1,0), tcl=-0.1)
	processed[,plot(unobs_rich[-length(unobs_rich)], n_col[-1], xlab="Unobserved species present last year", ylab="Species colonizing this year")]
	mtext(reg, side=3)
	abline(a=0, b=1)
	
	Figures[[reg]][['Figure6']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure6']][["name"]] <- fig6_name
	Figures[[reg]][['Figure6']][["dim"]] <- fig6_dim
	
	return(Figures)
}