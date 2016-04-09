#' Plot the Temperature Rank of Colonizers and Leavers
#' 
#' For each category of species, plot the temperature rank and the number of species in that category
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
plot_rank_temp <- function(prn, Figures){
	requireNamespace("beanplot", quietly=TRUE)
	
	unpack_p(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
	
	tr <- rank_temp(rd)
	tr2 <- tr[[2]]
	tr2[spp%in%spp_col & !spp%in%spp_ext, status:="colonizer"]
	tr2[spp%in%spp_ext & !spp%in%spp_col, status:="leaver"]
	tr2[spp%in%spp_ext & spp%in%spp_col, status:="both"]
	tr2[!spp%in%spp_ext & !spp%in%spp_col, status:="neither"]
	
	f_dim <- c(3.5, 6)
	dev.new(width=f_dim[1], height=f_dim[2])
	par(mfrow=c(2,1), mar=c(2,2,1,0.1), mgp=c(1,0.1,0), tcl=-0.1, cex=1, ps=8)
	tr2[,j={
		nc <- sum(status=="colonizer")
		nl <- sum(status=="leaver")
		nb <- sum(status=="both")
		nn <- sum(status=="neither")
		n <- c(nc, nl, nb, nn)
		tot <- sum(n)
		barplot(n/tot, names.arg=c("colonizer", "leaver","both","neither"), main=reg, ylab="Proportion of Species")
	}]
	# tr2[,j={boxplot(bt_mean_rank~status, ylab="Species Temperature Rank", main=reg);NULL}]
	
	bLine <- rainbow(n=5, v=0.8, s=1)
	names(bLine) <- c("colonizer", "leaver","both","neither","blah")
	bFill <- rgb(t(col2rgb(bLine, alpha=TRUE)), alpha=40, maxColorValue=255)
	names(bFill) <- c("colonizer", "leaver","both","neither","blah")
	beanCol <- list(
		colonizer = c(bFill[1]),
		leaver = c(bFill[2]),
		both = c(bFill[3]),
		neither = c(bFill[4])
	)
	
	bFill <- bFill[names(bFill)%in%tr2[,una(status)]]
	beanCol <- beanCol[names(beanCol)%in%tr2[,una(status)]]
	bLine <- bLine[names(bLine)%in%tr2[,una(status)]]
	
	tr2[,j={beanplot(bt_mean_rank~status, ylab="Species Temperature Rank", main=reg, border=bLine, col=beanCol, ll=0.01, beanlinewd=1.5);NULL}] 
	
	
	Figures[[reg]][['Figure9']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure9']][["name"]] <- paste0("nCat_tRank", reg, ".png")
	Figures[[reg]][['Figure9']][["dim"]] <- f_dim
	
	return(Figures)
}