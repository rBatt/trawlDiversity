plot_btemp_map <- function(prn, Figures){
	unpack_p(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
	fig3_name <- paste0("btempMap_", reg, ".png")
	f3_mfrow <- auto.mfrow(n_yrs)
	f3_height <- 6
	f3_width <- f3_mfrow[2]*f3_height/f3_mfrow[1]
	fig3_dim <- c(f3_width, f3_height)
	
	dev.new(width=f3_width, height=f3_height)
	par(mfrow=f3_mfrow, oma=c(0.1,0.1, 1,0.1), mar=c(1,1,0.1,0.1), mgp=c(0.75,0.1,0), tcl=-0.15, cex=1, ps=8)
	
	bt[,j={
		plot(lon, lat, type="n")
		map(add=TRUE)
		points(lon, lat, col=bt_col, pch=20)
		mtext(unique(year), side=3, adj=0.1, line=-0.75, font=2)
	}, by="year"]
	mtext(paste(reg, "Bottom Temperature"), outer=TRUE, side=3, line=0.25, font=2)
	
	Figures[[reg]][['Figure3']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure3']][["name"]] <- fig3_name
	Figures[[reg]][['Figure3']][["dim"]] <- fig3_dim
	
	return(Figures)
}