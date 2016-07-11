#' Colonization and Extinction Plotting Wrapper
#' 
#' Creates 3 plots for each species involved in some combination of colonization and extinction
#' 
#' @param prn the p object (processed msom; output from \code{process_msomStatic})
#' @param Figures option list to which the figure and its information should be added
#' @param spp_cat Category for colonization/ extinction to be chosen
#' @param width.max max width of figure in inches
#' @param height.max max height of figure in inches
#' @param max_spp_columns max number of columns
#' @param nPlots max number of plots per species
#' @param FUN the graphical device function
#' @param ... arguments passed to plot_device
#' 
#' @details
#' This figure is rather complicated, and is divided into 3 parts. This function is the piece that should be interacted with. The core piece of the function is \code{\link{plot_ce}}, which creates 3 core plots.
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
plot_ce_wrap <- function(prn, Figures, spp_cat=c("col","ext","both","neither"), width.max=12, height.max=7, max_spp_columns=15, nPlots=3, FUN="dev.new", ...){
	unpack_p(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
		
	spp_cat <- match.arg(spp_cat)
	fig_ind <- which(c("col","ext","both","neither")==spp_cat)
	
	fig_num <- paste0('Figure8.',fig_ind)
	
	fig_name <- c(
		"col" = paste0("who_colonized_", reg, ".png"),
		"ext" = paste0("who_left_", reg, ".png"),
		"both" = paste0("who_colonized_andLeft_", reg, ".png"),
		"neither" = paste0("who_no_col_nor_ext_", reg, ".png")
	)[spp_cat]
	
	spp2use <- list(
		"col" = spp_col_only,
		"ext" = spp_ext_only,
		"both" = spp_col_and_ext,
		"neither" = spp_neither
	)[[spp_cat]]
	
	main_title <- c(
		"col" = "Colonizing Species",
		"ext" = "Leaving Species",
		"both" = "Species that Colonize and Go Extinct",
		"neither" = "Species that Are Always Present"
	)[spp_cat]
	
	if(length(spp2use)==0){
		message("No species in this category")
		
		Figures[[reg]][[fig_num]][["figure"]] <- list()
		Figures[[reg]][[fig_num]][["name"]] <- fig_name
		Figures[[reg]][[fig_num]][["dim"]] <- c(3.5,3.5)
		
		Figures[[reg]][[fig_num]] <- plot_device(Figures[[reg]][[fig_num]], FUN, ...)
		
		par(mar=c(0.5,0.5,0.5,0.5), cex=1, ps=10)
		plot(1, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
		text(1,1, paste0("This Region Has 0\n", main_title))
		
		
	}else{
		ncolspp <- max_spp_columns
		nrowspp <- ceiling(length(spp2use)/ncolspp)
		lmat0 <- matrix(rep(1:(ncolspp*nrowspp),each=nPlots), nrow=nrowspp*nPlots, ncol=ncolspp)
		lmat <- lmat0
		lmat[] <- 1:prod(dim(lmat))
		lmat <- lmat[,apply(lmat < (length(spp2use)*nPlots), 2, function(x)any(x))]

		fwc_ldim <- dim(lmat)
		fwc_w <- min(max(1*fwc_ldim[2], 1.5), width.max) # width
		fwc_h <- min(max(fwc_w * fwc_ldim[1] / fwc_ldim[2] * 1.1, 3.0), height.max)
		fwc_dim <- c("width"=fwc_w, "height"=fwc_h)
		
		Figures[[reg]][[fig_num]][["figure"]] <- list()
		Figures[[reg]][[fig_num]][["name"]] <- fig_name
		Figures[[reg]][[fig_num]][["dim"]] <- fwc_dim
	
		Figures[[reg]][[fig_num]] <- plot_device(Figures[[reg]][[fig_num]], FUN, ...)
	
		layout(mat=lmat, heights=rep(c(1.5,1,1), nrowspp))
		par(mar=c(1,1,0.75,0.1) ,oma=c(0.1,0.1,0.5,0.1), ps=6, mgp=c(0.5, 0.01, 0), tcl=-0.01, cex=1)
	
		for(i in 1:length(spp2use)){
			t_sco <- spp2use[i]
			plot_ce(t_sco, pad_top_mar=ifelse(length(spp2use)>max_spp_columns, 2, 1), plt_pts=FALSE, use_ext=ifelse(spp_cat=="ext", TRUE, FALSE))
		}
		mtext(main_title, side=3, line=-0.5, outer=TRUE, font=2)
	
	
	}
	
	if(is.null(Figures[[reg]][[fig_num]][["fig_loc"]])){
		Figures[[reg]][[fig_num]][["figure"]] <- recordPlot()
	}else{
		Figures[[reg]][[fig_num]][["figure"]] <- NULL
		dev.off()
	}
	
	return(Figures)
	
}