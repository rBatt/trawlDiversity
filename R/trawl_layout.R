#' Trawl Layout
#' 
#' Return a matrix appropriate for the \code{layout} of maps for 9 regions that are to scale
#' 
#' @details
#' The matrix returned is sized such that the longitude:latitude ratio is not distorted within a region (approximately). 
#' Furthermore, the sizes of the regions are to scale. 
#' One exception is that the Southeast US takes up a 3x3 grid instead of a 2x2, making is roughly double the size it should be (although lon-lat in this region are still in proportion). 
#' Also, the regions are arranged roughly according to their true geographic arrangement.
#' 
#' The functions expects that regions will be plotted in a specific order. That order is: \code{c('ebs','ai','goa','wctri','gmex','sa','neus','shelf','newf')}. 
#' Note that "wcann" is excluded. It cannot be simply substituted in place of "wctri", either, because the latitudinal extent is greater for "wcann" than for "wctri".
#' 
#' @examples
#' # subset data to unique lon-lat
#' d_ll <- data_all[,list(reg,lon,lat)]
#' setkey(d_ll, reg, lon, lat)
#' d_ll <- unique(d_ll)
#' 
#' # device dimension ~ layout matrix dimension
#' # i.e., 3:7 (height:width)
#' dev.new(height=3, width=7)
#' par(mar=c(0.5,0.5,0.5,0.5), mgp=c(0.5,0.1,0), tcl=-0.1, cex=1, ps=6)
#' layout(trawl_layout())
#' for(r in 1:9){
#' 	tr <- c('ebs','ai','goa','wctri','gmex','sa','neus','shelf','newf')[r]
#' 	d_ll[reg==tr,plot(lon,lat,main=reg[1])]
#' }
#' 
#' @export
trawl_layout <- function(){
	lay_grid <- matrix(1:84, nrow=6, ncol=14)
	squares <- list(
		ebs = c(1,26),
		ai = c(6,42),
		goa = c(3,41),
		wctri = c(45,48),
		gmex = c(53,66),
		sa = c(70,84), # should really be c(70,77), ... just made it bigger to fill in gap and b/c it's a square
		neus = c(67,81),
		shelf = c(31,44),
		newf = c(49,64)
	)
	squares_ind <- lapply(squares, function(x)arrayInd(which(lay_grid%in%x),.dim=dim(lay_grid)))
	map_layout <- array(NA, dim(lay_grid))
	for(i in 1:length(squares_ind)){
		map_layout[squares_ind[[i]]] <- i
	}
	map_layout[is.na(map_layout)] <- 0
	
	return(map_layout)
}