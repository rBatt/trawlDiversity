#' Update MSOM results from JAGS
#' 
#' Uses autojags to update msom reuslts until converged
#' 
#' @param out output of a jags model
#' @param arguments to pass to \code{autojags}
#' 
update_msom <- function(out, ...){
	requireNamespace("R2jags", quietly = TRUE)
	
	R2jags::recompile(out)
	updated_out <- autojags(out, ...)
	return(updated_out)
}

#' @describeIn update_msom
get_update_msom <- function(file, r, n.iter, n.thin ...){
	requireNamespace("R2jags", quietly = TRUE)
	
	load(file) # lods rm_out
	
	if(missing(n.iter)){
		n.iter <- 1E3
	}
	if(missing(n.thin)){
		n.thin <- max(1, floor((n.iter/2)/200))
	}
	
	for(i in 1:length(rm_out)){
		rm_out[[r]][[i]]$out <- update_msom(rm_out[[r]][[i]]$out, n.iter=n.iter, n.thin=n.thin, ...)
	}
	
}