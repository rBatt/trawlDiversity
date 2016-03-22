#' Plot Spatial Data
#' 
#' Plot spatial data consistent of X, Y, and Z components
#' 
#' @param x x values (independent)
#' @param y y values (independent)
#' @param z z values (dependent value)
#' @param scatter logical, if true, plots a scatter plot of raw values; else, smooths data with Krigging
#' @param ... arguments passed to plot via colScat
#' 
#' @export
plot_space <- function(x, y, z, scatter=FALSE, ...){
	
	requireNamespace("fields", quietly=TRUE)
	requireNamespace("rbLib", quietly=TRUE)
	requireNamespace("reshape2", quietly=TRUE)
	
	colScat <- function(x, y, z, ...){
		zols <- rbLib::zCol(256, z)
		d.ord <- order(z)
		plot(x[d.ord], y[d.ord], col=zols[d.ord], ...)
	}
	
	if(scatter){
		colScat(x, y, z, ...)
	}else{
	
		smooth_hat_fit <- fields::Tps(matrix(c(x,y), ncol=2), Y=z, lon.lat=TRUE, scale.type="unscaled")
		# smooth_hat_fit <- n_spp_col_weighted_tot[,Krig(matrix(c(lon,lat), ncol=2), Y=n_spp_col_weighted)]
		smooth_hat <- fields::predictSurface(smooth_hat_fit)
		smooth_hat_arr <- array(smooth_hat$z, dim=dim(smooth_hat$z), dimnames=list(lon=smooth_hat$x, lat=smooth_hat$y))
		smooth_hat_melt <- data.table(reshape2::melt(smooth_hat_arr))
		smooth_hat_melt[,colScat(lon, lat, value, ...)]
	}

	
	invisible(NULL)
}