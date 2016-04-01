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