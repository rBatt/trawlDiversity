#' Colonization and Extinction Plotting Wrapper
#' 
#' Creates 3 plots for each species involved in some combination of colonization and extinction
#' 
#' @param prn the p object (processed msom; output from \code{process_msomStatic})
#' @param Figures option list to which the figure and its information should be added
#' @param spp_cat Category for colonization/ extinction to be chosen
#' 
#' @details
#' This figure is rather complicated, and is divided into 3 parts. This function is the piece that should be interacted with. The core piece of the function is \code{\link{plot_ce}}, which creates 3 core plots. Setting up the figure layout and size can also be tricky, but is accomplished by \code{\link{plot_ce_setup}}.
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
plot_ce_wrap <- function(prn, Figures, spp_cat=c("col","ext","both","neither"), ...){
	unpack_p(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
	spp_cat <- match.arg(spp_cat)
	fig_ind <- which(c("col","ext","both","neither")==spp_cat)
	
	fig_name <- c(
		"col" = paste0("who_colonized", reg, ".png"),
		"ext" = paste0("who_left", reg, ".png"),
		"both" = paste0("who_colonized_andLeft", reg, ".png"),
		"neither" = paste0("who_no_col_nor_ext", reg, ".png")
	)[spp_cat]
	
	spp2use <- list(
		"col" = spp_col_only,
		"ext" = spp_ext_only,
		"both" = spp_col_and_ext,
		"neither" = spp_neither
	)[[spp_cat]]
	
	if(length(spp2use)==0){
		message("No species in this category")
		dev.new(width=3.5,height=3.5)
		plot(1, type="n")
		
		Figures[[reg]][[paste0('Figure8.',fig_ind)]][["figure"]] <- recordPlot()
		Figures[[reg]][[paste0('Figure8.',fig_ind)]][["name"]] <- fig_name
		Figures[[reg]][[paste0('Figure8.',fig_ind)]][["dim"]] <- c(3.5,3.5)
		
		return(Figures)
	}
	
	fd <- plot_ce_setup(spp2use, ...)
	for(i in 1:length(spp2use)){
		t_sco <- spp2use[i]
		plot_ce(t_sco, pad_top_mar=ifelse(length(spp2use)>15, 2, 0), plt_pts=FALSE, use_ext=ifelse(spp_cat=="ext", TRUE, FALSE))
	}
	
	
	Figures[[reg]][[paste0('Figure8.',fig_ind)]][["figure"]] <- recordPlot()
	Figures[[reg]][[paste0('Figure8.',fig_ind)]][["name"]] <- fig_name
	Figures[[reg]][[paste0('Figure8.',fig_ind)]][["dim"]] <- fd
	
	return(Figures)
	
}