#' Make owin from Region Outline
#' 
#' Make an owin object for spatstat package
#' 
#' @param X a data.table of mapDat, e.g.
#' @param outlines region outlines, e.g., the outlines object
#' 
#' @examples
#' trawlDiversity::make_owin(mapDat, outlines)
#' 
#' @export
make_owin <- function(X, outlines){
	requireNamespace("spatstat", quietly = TRUE)
	
	rs <- X[,una(reg)]
	nr <- length(rs)
	X_owin <- list()
	for(r in 1:nr){
		td <- X[reg==rs[r]]
		o <- outlines[reg==rs[r], list(x=lonP, y=latP)]
		o_r <- outlines[reg==rs[r], list(x=rev(lonP), y=rev(latP))]
	
		xr <- td[,range(c(lon,o[,x]), na.rm=TRUE)]
		yr <- td[,range(c(lat,o[,y]), na.rm=TRUE)]
		ow_exp1 <- bquote(spatstat::owin(xr, yr, poly=o)) # use when I correctly traced outline counterclockwise
		ow_exp2 <- bquote(spatstat::owin(xr, yr, poly=o_r)) # reversed (use if traced outline in wrong direction)
	
		X_owin[[rs[r]]] <- tryCatch(td[, eval(ow_exp1)], error=function(cond)td[, eval(ow_exp2)])
	}
	return(X_owin)
}