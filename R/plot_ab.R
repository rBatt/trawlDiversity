#' Plot Alpha and Beta
#' 
#' Plots a species image and the parameter in a provided data.table
#' 
#' @param X a data.table with columns for year and value
#' @param t_spp character indicating a species name
#' @param plt_img logical indicating whether or not to plot a picture of the species in the background
#' @param plt_pts logical indicating whether or not to plot the posterior points, or just the mean of the posterior
#' @param ... arguments passed to both \code{sppImg} and \code{plot}
#' 
#' @details 
#' Used in process_msom_figures.R. X is assumed to have been subset to the specific species and the desired annual posterior for value, such as when value is a species-specific paraemter for alpha or beta.
#' 
#' The function tries to adjust the opacity of the points according to the size of the panel, with the intention of not completely obscuring the background image and so that opacity can be used as an indicator of posterior density. However, if the panel is too small, the transparency may be truncated to 0, so the posterior dots are just invisible.
#' 
#' @return
#' the output of \code{sppImg}, or if \code{plt_img} is \code{FALSE}, returns \code{NULL}
#' 
#' @export
plot_ab <- function(X, t_spp, plt_img=TRUE, plt_pts=TRUE, ...){

	if(plt_img){
		si <- sppImg(t_spp, ...)
	}else{
		si <- NULL
	}
	if(!is.null(si)){
		par(new=T)
	}
	
	fin <- par("fin")[2]
	fac <- 0.01*fin^3/2
	
	nX <- names(X)
	has_full_dist <- "value" %in% nX
	if(has_full_dist){
		mu <- X[,list(mu=mean(value)),by="year"]
	}else{
		mu <- X[,list(mu=mean(value_mu)),by="year"] # the aggregation shouldn't actually be necessary
	}
	if(plt_pts & has_full_dist){
		plot(X[,year], X[,value], col=adjustcolor('gray', fac), cex=0.5, pch=21, bg=adjustcolor('white',fac), ...)
		mu[,lines(year, mu, lwd=2, col='gray')]
	}else{
		mu[,plot(year, mu, lwd=2, col='gray', type="l", ...)]
	}
	mu[,lines(year, mu, lwd=1, col='white')]
	if(is.null(si) & plt_img){
		common_name <- spp.key[spp==t_spp, una(common)]
		mtext(paste(t_spp, common_name, sep="\n"), side=3)
	}
	
	return(si)
}