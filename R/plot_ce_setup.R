#' Setup Colonization & Extinction Plot
#' 
#' Creates layout and dimensions for c/e plot
#' 
#' @param spp2use character vector of species names
#' @param width.max the maximum width of the figure in inches
#' @param height.max the maximum height of the figure in inches
#' @param max_spp_columns the maximum number of columns to use when arranging panels
#' @param nPlots the number of plots per species; determined by \code{\link{plot_ce}}
#' 
#' @details
#' Shouldn't need to use this function directly. See \code{\link{plot_ce_wrap}}
#' 
#' @return
#' Returns a numeric vector of length 2 indicating the width and height of the output figure, in inches
#' 
#' @seealso 
#' There are several functions that are run through the process_msom_figures script. Richness and temperature plots are \code{\link{plot_btemp_map}}, \code{\link{plot_rich_bt_scatter}}, and \code{\link{plot_rich_bt_ts}}. Figures for colonization, extinction, and the species and places associated with those processes are \code{\link{plot_ce_wrap}}, \code{\link{plot_col_vs_unobsSpp}}, \code{\link{plot_colExt_perStrat}}, and \code{\link{plot_rank_temp}}. Figures for diagnostics are \code{\link{plot_traceplot}} and \code{\link{plot_post_corr}}.
#' 
#' @export
plot_ce_setup <- function(spp2use, width.max=12, height.max=7, max_spp_columns=15, nPlots=3){
	ncolspp <- max_spp_columns
	
	
	# if(length(spp2use) > ncolspp){
		nrowspp <- ceiling(length(spp2use)/ncolspp)
		lmat0 <- matrix(rep(1:(ncolspp*nrowspp),each=nPlots), nrow=nrowspp*nPlots, ncol=ncolspp)
		lmat <- lmat0
		lmat[] <- 1:prod(dim(lmat))
		lmat <- lmat[,apply(lmat < (length(spp2use)*nPlots), 2, function(x)any(x))]
	
		fwc_ldim <- dim(lmat)
		fwc_w <- min(max(1*fwc_ldim[2], 1.5), width.max) # width
		fwc_h <- min(max(fwc_w * fwc_ldim[1] / fwc_ldim[2] * 1.1, 3.0), height.max)
		fwc_dim <- c("width"=fwc_w, "height"=fwc_h)
		dev.new(width=fwc_dim[1], height=fwc_dim[2])
		layout(mat=lmat, heights=rep(c(1.5,1,1), nrowspp))
		par(mar=c(1,1,0.75,0.1) ,oma=c(0.1,0.1,0.5,0.1), ps=6, mgp=c(0.5, 0.01, 0), tcl=-0.01, cex=1)
	#
	# }else{
	# 	fwc_mfr <- c(nPlots, length(spp2use))
	# 	fwc_w <- min(max(1*fwc_mfr[2], 1.5), width.max) # width
	# 	fwc_h <- min(max(fwc_w * fwc_mfr[1] / fwc_mfr[2] * 1.1, 3.0), height.max)
	# 	fwc_dim <- c(width=fwc_w, height=fwc_h)
	# 	dev.new(width=fwc_dim[1], height=fwc_dim[2])
	# 	par(mfcol=fwc_mfr, mar=c(1,1,0.75,0.1) ,oma=c(0.1,0.1,0.5,0.1), ps=6, mgp=c(0.6, 0.1, 0), tcl=-0.1, cex=1)
	# }
	
	return(c(fig_dim=fwc_dim))

}