#' Plot Device
#' 
#' Prepare the graphical device for use with processed msom figures
#' 
#' @param x a Figures list
#' @param FUN a character indicating a graphical device function
#' @param path folder where the figure will be saved, if FUN is not \code{dev.new}
#' @param ... other arguments passed to \code{FUN}
#' 
#' @details
#' Used in the functions called in the plotting functions used in process_msom_figures.R
plot_device <- function(x, FUN, path="trawlDiversity/pkgBuild/figures/pmfPieces", ...){
	fig_fun <- match.fun(FUN)
	d <- x$dim

	if(FUN!="dev.new"){
		fig_name <- gsub("\\.png", paste0(".",FUN), x$name)
		fig_loc <- file.path(path, fig_name)
		fig_fun(width=d[1], height=d[2], file=fig_loc, ...)
		x$fig_loc <- fig_loc
	}else{
		fig_fun(width=d[1], height=d[2], ...)
		x$fig_loc <- NULL
	}
	return(x)
}