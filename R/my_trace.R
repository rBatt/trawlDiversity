#' Posterior Traceplots
#' 
#' Generate traceplots from either JAGS or Stan output
#' 
#' @param x model output, or a data.table result from the output of \code{\link{get_iters}}
#' @param pars names of parameters to be plotted. Uses \code{get_iters} defaults for getting parameters.
#' @param lang model language (of \code{x})
#' @param ... arguments passed to plot
#' 
#' @details
#' Easier to handle than other traceplots (if, e.g., wanted to add trace plot from multiple models to same figure, hard to do w/ Stan function). Works for JAGS or Stan output (via R2jags and rstan).
#' 
#' @export
mytrace <- function(x, pars, lang, ...){
	
	if(!any(class(x)=="data.table")){
		sims <- get_iters(x, pars, lang=lang)
	}else{
		sims <- x
	}
	
	cn <- pars #colnames(sims)
	for(h in 1:length(cn)){
		ylim <- sims[,range(eval(s2c(cn[h]))[[1]])]
		n_chains <- sims[,lu(chain)]
		cols <- adjustcolor(1:n_chains, 0.5)
		for(i in 1:n_chains){
			
			if(i == 1){
				sims[chain==i, plot(eval(s2c(cn[h]))[[1]], ylim=ylim, type="l", col=cols[i], xlab="",ylab=cn[h], ...)]
			}else{
				sims[chain==i, lines(eval(s2c(cn[h]))[[1]], col=cols[i])]
			}
		}
	}
}