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