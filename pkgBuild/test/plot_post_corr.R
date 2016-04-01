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