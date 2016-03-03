#' Get sampled iterations of posterior
#' 
#' 

get_iters <- function(X, pars, lang){
	
	if(lang=="JAGS"){
		# sims <- x$BUGSoutput$sims.list[pars]
		sims0 <- lapply(apply(X$BUGSoutput$sims.array[,,pars], 2, function(x)list(as.matrix(x))), function(x)x[[1]])
		sims <- vector("list", length(sims0))
		for(i in 1:length(sims)){
			sims[[i]] <- as.data.table(cbind(sims0[[i]], chain=rep(i, nrow(sims0[[i]]))))
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