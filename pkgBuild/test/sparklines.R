#' Add Sparklines to a Figure
#' 
#' Adds lines and (optionally) axes positioned relative to reference coordinates and scaled to the size of a character
#' 
#' @param x,y numeric vector of values for sparklines
#' @param x_pt,y_pt numeric length 1 indicating coordinate to use as a reference for positioning sparkline
#' @param y_align character indicating the portion of \code{y} to use for vertical alignment relative to reference coordinate
#' @param x_align, same as \code{y_align}, except for \code{x} and horizontal alignment; fewer options than \code{y_align}
#' @param scale_xy Logical; if TRUE (default) \code{x} and \code{y} are rescaled to be the size of a character
#' @param x_cex,y_cex numeric, scaling factor to adjust the size of the sparklines in either direction; does nothing if \code{scale_xy} is FALSE
#' @param ax_sides Integer value(s) in 1:4 indicating the sides of the sparkline to be bordered by a solid 'axis' line; if \code{NULL} (default), no lines are drawn
#' @param col color of the sparkline
#' @param lwd width of the sparkline
#' @param acol color of the axis lines
#' @param awd width of the axis lines
#' 
#' @details
#' Created as a way of associating a mini time series with a particular point in a primary time series. For example, if plotting monthly averages, could add a sparkline for each month showing daily values.
#' 
#' @return Returns NULL invisibly
#' 
#' @examples
#' # fake data
#' x <- 1:20
#' y <- cumsum(rnorm(20))
#' 
#' # plot squared values
#' # to show that sparkline x-y 
#' # can be on a different scale
#' plot(x^2, y^2, type='l')
#' 
#' # add sparklines
#' # showing a mini version of the time series
#' # leading into and out of each point
#' # of the big time series
#' for(i in 1:length(x)){
#' 	sparklines(x, y, y_pt=y[i]^2, x_pt=x[i]^2, y_align="right", x_align="right", ax_sides=c(1,2), col="red", awd=0.5)
#' 	sparklines(x, y, y_pt=y[i]^2, x_pt=x[i]^2, y_align="left", x_align="left", ax_sides=c(1,4), col="blue", awd=0.5)
#' }
#' 
#' @export
sparklines <- function(x, y, x_pt, y_pt, y_align=c("bot","top","mid", "left", "right"), x_align=c("left","right","mid"), scale_xy=TRUE, x_cex=1, y_cex=1, ax_sides=NULL, col="black", lwd=1, acol="black", awd=1){
	
	stopifnot(is.null(ax_sides) || all(ax_sides%in%(1:4)))
	y_align <- match.arg(y_align)
	x_align <- match.arg(x_align)
	
	
	scale_fac <- function(type=c("x","y")){
		type <- match.arg(type)
		c1 <- mean(par("cin"))/par("pin")
		switch(type, x=c1[1]*diff(par()$usr[1:2])*x_cex, y=c1[2]*diff(par()$usr[3:4])*y_cex)
	}
	
	set_xy <- function(type=c("x","y")){
		type <- match.arg(type)
		v <- switch(type, x=x, y=y)
		v_align <- switch(type, x=x_align, y=y_align)
		v_pt <- switch(type, x=x_pt, y=y_pt)
		
		if(scale_xy){
			vn <- (v-min(v))/max(abs((v-min(v)))) * scale_fac(type)
		}else{
			vn <- v
		}
		
		va_val <- switch(v_align,
			left = vn[1],
			bot = min(vn),
			right = vn[length(vn)],
			top = max(vn),
			mid = mean(range(vn))
		)
		
		vn <- vn - (va_val - v_pt)
		return(vn)
	}
	
	xs <- set_xy("x")
	ys <- set_xy("y")
	
	if(!is.null(ax_sides)){
		xadj <- (diff(range(xs))/1E2)*c(-1,-1,1,1) # adjust 'axis' lines to be slightly outside the range of values being plotted
		yadj <- (diff(range(ys))/1E2)*c(-1,1,1,-1)
		xy_corners <- expand.grid(y=range(ys), x=range(xs))[c(1,2,4,3),2:1] # rows are: bottom-left, top-left, top-right, bottom-right
		xyc_adj <- xy_corners + cbind((xadj),(yadj)) # corners, slightly extended outward (adjusted)
		adf <- data.frame(x0=xyc_adj[c(1,1,2,3),1], x1=xyc_adj[c(4,2,3,4),1], y0=xyc_adj[c(1,1,2,3),2], y1=xyc_adj[c(4,2,3,4),2])
		
		for(s in 1:length(ax_sides)){
			do.call('segments', c(adf[ax_sides[s],],list(col=acol, lwd=awd)))
		}
	}
	
	lines(xs, ys, col=col, lwd=lwd)
	
	invisible(NULL)
}









