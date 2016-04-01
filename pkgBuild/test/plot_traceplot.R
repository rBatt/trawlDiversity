plot_traceplot <- function(prn, Figures){
	unpack_p(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
	fig4_name <- paste0("traceplot_", reg, ".png")
	f4_mfrow <- c(n_pars, n_yrs)
	f4_height <- 5
	f4_width <- f4_mfrow[2]*f4_height/f4_mfrow[1]
	fig4_dim <- c(f4_width, f4_height)

	dev.new(width=f4_width, height=f4_height)
	par(mfrow=f4_mfrow, oma=c(1,1, 1,0.1), mar=c(0.5,0.5,0.1,0.1), mgp=c(0.25,0.1,0), tcl=-0.1, cex=1, ps=6)
	
	for(h in 1:length(pars_trace)){
		for(i in 1:n_yrs){
			t_yr <- param_iters[,una(year)][i]
			t_iters <- param_iters[year==t_yr]
			mytrace(t_iters, pars=pars_trace[h], lang=lang, xaxt='n')
			if(i == 1){
				mtext(pars_trace[h], side=2, line=0.75)
			}
		}
	}
	
	Figures[[reg]][['Figure4']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure4']][["name"]] <- fig4_name
	Figures[[reg]][['Figure4']][["dim"]] <- fig4_dim
	
	return(Figures)
}