#' Resample Observed Local Richness
#' 
#' The number of tows affects richness, this function calculates site-level richness after rarefying tow intensity
#' 
#' @param X a data.table that is a subset of \code{\link{trawlDiversity::data_all}}. Subset to a region-year.
#' 
#' @details 
#' Value returned is the rarefied (to 1 tow) cross-site average richness. Each site in the region is expected, on average, to have this many species observed in it if 1 tow is made per site.
#' 
#' @return returns a numeric value of length 1, the average local richness 
#' 
#' @seealso \code{\link{range_sample}}
#' 
#' @export
local_rich_sample <- function(X){
	format_local_rich <- function(Xdat){
		uspp <- data.table(spp=Xdat[,unique(spp)])
		site_tow_spp <- Xdat[,CJ(K=K, spp=spp, unique=TRUE), by='stratum']
		setkey(site_tow_spp)
		X <- Xdat[site_tow_spp, list(stratum,K,spp,abund), on=c("stratum","K","spp")]
		X[is.na(abund), abund:=0L]
		return(X)
	}
	
	X <- format_local_rich(Xdat=X)
	X[,abund:=as.integer(abund)]
	setkey(X, stratum, spp)
	
	nStrat <- X[,length(unique(stratum))]
	
	innerF_calc <- function(x, nT=1){
		x[,j={
			nk <- nrow(.SD)
			p <- sum(abund)/nk
			prob <- 1 - (1-p)^nT
			list(prob=prob)
		},by=c('stratum','spp')]
	}
	
	probs <- innerF_calc(X)
	site_rich <- probs[,list(local_rich_samp=sum(prob)),by=c("stratum")]
	local_rich_samp <- site_rich[,mean(local_rich_samp)]
	return(local_rich_samp)
}