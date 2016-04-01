plot_col_vs_unobsSpp <- function(prn, Figures){
	unpack_p(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
	fig6_name <- paste0("Colonization_UnobsSpp_", reg, ".png")
	fig6_dim <- c(3.5, 3.5)
	
	dev.new(fig6_dim[1], fig6_dim[2])
	
	processed[,plot(unobs_rich[-length(unobs_rich)], n_col[-1], xlab="Unobserved species present last year", ylab="Species colonizing this year")]
	mtext(reg, side=3)
	abline(a=0, b=1)
	
	Figures[[reg]][['Figure6']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure6']][["name"]] <- fig6_name
	Figures[[reg]][['Figure6']][["dim"]] <- fig6_dim
	
	return(Figures)
}