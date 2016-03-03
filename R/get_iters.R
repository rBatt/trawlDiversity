#' Get sampled iterations of posterior
#' 
#' Given model output from either Stan or JAGS, return a data.table of posterior samples for specified parameters
#' 
#' @param X an object containing model output, e.g., as returned by \code{stan} or \code{jags}
#' @param pars Character vector of parameter names. Note that naming coventions differ slightly for Stan and JAGS
#' @param lang character vector indicating the language that generated X (either "JAGS" or "Stan")
#' 
#' @details important to remember that models in the two languages might have different parameters/ parameter names. For example, if there is only one \code{beta_mu} parameter, it is \code{beta_mu} in JAGS and \code{beta_mu[1]} in Stan. Also, \code{logit_psi} and \code{logit_theta} are tracked in Stan but not JAGS, and \code{Z} parameters are only in JAGS. Both the way the models are written and the languages themselves contribute to these differences.
#' 
#' @return
#' A data.table with columns for each of \code{pars} and a column for \code{chain}, the latter of which indicates which of the chains the sample is from.
#' 
#' @export
get_iters <- function(X, pars, lang=c("JAGS","Stan")){
	
	lang <- match.arg(lang)
	
	if(lang=="JAGS"){
		# sims <- x$BUGSoutput$sims.list[pars]
		sims0 <- lapply(apply(X$BUGSoutput$sims.array[,,pars], 2, function(x)list(as.data.table(x))), function(x)x[[1]])
		sims <- vector("list", length(sims0))
		for(i in 1:length(sims)){
			# sims[[i]] <- as.data.table(cbind(sims0[[i]], chain=rep(i, nrow(sims0[[i]]))))
			sims[[i]] <- sims0[[i]][,chain:=rep(i, nrow(sims0[[i]]))]
			sims[[i]] <- setnames(sims[[i]], old=colnames(sims[[i]]), new=c(pars, "chain"))
		} 
	}else if(lang == "Stan"){
		sims <- lapply(X@sim$samples, function(x)x[pars])
		for(i in 1:length(sims)){
			sims[[i]]$chain <- rep(i, length(sims[[i]][[1]]))
		}
	}
	
	sims <- rbindlist(sims)
	
	return(sims)
	
}