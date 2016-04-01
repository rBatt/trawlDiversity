plot_rich_bt_ts <- function(prn, Figures){
	unpack_p(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
	fig1_name <- paste0("richness_bt_timeSeries_", reg, ".png")
	fig1_dim <- c(3.5, 6)
	
	dev.new(width=fig1_dim[1], height=fig1_dim[2])
	par(mfrow=c(3,1), mar=c(1.75,1.5,0.25,0.25), oma=c(0.1,0.1,0.75,0.1), mgp=c(0.75,0.1,0), tcl=-0.1, ps=8, cex=1)
	
	plot(naive_rich, type="o", ylab="Naive Region Richness", xlab="Year")
	plot(reg_rich, type="o", xlab="Year", ylab="MSOM Region Richness")
	plot(bt_ann, type="o", xlab="Year", ylab="Annual Mean Bottom Temperature")
	mtext(reg, side=3, line=0, outer=TRUE, font=2)
	
	Figures[[reg]][['Figure1']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure1']][["name"]] <- fig1_name
	Figures[[reg]][['Figure1']][["dim"]] <- fig1_dim
	
	return(Figures)
	
}