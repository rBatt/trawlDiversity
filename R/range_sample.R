#' Format data to calculate range
#' Reformat data for use in resampling to estimate range size
#' @param Xdat a data.table that is a subset of \code{trawlDiversity::data_all}; subset to one region-year
#' @param sppName the name of a species; corresponds to the "spp" column in \code{trawlDiversity::data_all}
format_range <- function(Xdat, sppName){
	site_tow_combos <- Xdat[,list(K=unique(K)), by='stratum']
	X <- Xdat[(spp)==sppName][site_tow_combos, list(stratum,K,abund) ,on=c("stratum","K")] # could be faster if I didn't do the species subsetting this way; faster to get formatted for all species, then re-write range_sample to automatically loop through species. Bleh.
	X[is.na(abund), abund:=0]
	return(X)
}

#' Resample the observed range of species occurrences
#' 
#' Resample tows among sites to determine range size while accounting for variability in tow frequency
#' 
#' @param X A matrix with rows being sites (indicating range size) and columns as tows within a site. Elements of matrix are taken to be binary: 0 for absent, 1 for present
#' @param sppName character of length 1, corresponding to "spp" in \code{trawlDiversity::data_all}
#' @param tow_target an integer of length 1 indicating the number of tows that should be represented in the estimated range size. I.e., the range size 
#' @param calc logical, if TRUE (default), the range size is calculated directly, rather than using the resampling technique that achieves the equivalent result (when nPerm is large)
#' @param nPerm integer of length 1, number of permutations when resampling tows. I.e., simulated sample size behind returned value (and each point in rarefaction curve)
#' @param rarefy logical, default TRUE. Compute the rarefaction curve?
#' @param tow_max integer of length 1 indicating the maximum extent of the rarefaction curve when \code{rarefy} is TRUE; when NULL (default), \code{tow_max} is set to the maximum number of tows for any site in \code{X}
#' 
#' @details 
#' Resampling is done with replacement. Therefore, the rarefaction curve will not asymptote to the observed range size at the maximum number of tows (default value when \code{tow_max} is NULL). The asymptote will always be at the observed range size.
#' 
#' @return  
#' If \code{rarefy} is FALSE, returns a numeric of length 1 indicating the estimated range size as a fraction of the total number of rows (sites) in \code{X}.
#' Otherwise, returns same as above but with the attribute "rc", the Rarefaction Curve, that is a vector of \code{length(tow_max)} indicating the estimated range size for a number of tows equal to the index \code{seq_len(tow_max)}.
#' 
#' @seealso \code{\link{local_rich_sample}}
#' 
#' @examples
#' # ---- range size at 1 tow ----
#' Xdat <- data_all[reg=="shelf" & year==1985, list(stratum,spp,K,Kmax,abund)]
#' range_sample(X=Xdat, sppName="Leucoraja erinacea")
#' 
#' # ---- calculate and plot rarefaction curve ----
#' rc <- range_sample(X=Xdat, sppName="Leucoraja erinacea", rarefy=TRUE)
#' plot(attributes(rc)$rc, type='o')
#' 
#' # ---- calculate for all species in a region for 2 years ----
#' # calculate resampled range
#' test <- data_all[reg=='shelf' & year%in%(1989:1990),j={
#' 	uspp <- unique(spp)
#' 	tf <- function(x){
#' 		range_sample(X=.SD, sppName=x, nPerm=3)
#' 		# Note: low nPerm for fast example
#' 	}
#' 	rs <- sapply(uspp, tf)
#' 	data.table(spp=names(rs), range_size_samp=unname(rs))
#' },by=c("reg","year")]
#' 
#' # calculate observed range
#' obs <- data_all[reg=='shelf' & year%in%(1989:1990), j={
#' 	nS <- length(unique(stratum))
#' 	.SD[,j={
#' 		list(obs_rs=length(unique(stratum))/nS)
#' 	},by="spp"]
#' },by=c("reg","year")]
#' 
#' # merge and plot
#' obs_test <- merge(obs, test, by=c("reg","year","spp"))
#' obs_test[,plot(obs_rs, range_size_samp)]
#' abline(a=0, b=1)
#' 
#' @export
range_sample <- function(X, sppName, tow_target=1, calc=TRUE, nPerm=100, rarefy=FALSE, tow_max=NULL){
	
	X <- format_range(Xdat=X, sppName=sppName)
	X[,abund:=as.integer(abund)]
	setkey(X, stratum)
	
	if(is.null(tow_max) & isTRUE(rarefy)){tow_max <- X[,max(K)]}
	
	nStrat <- X[,length(unique(stratum))]
	
	if(calc){
		innerF_calc <- function(x, nT){
			x[,j={
				nk <- nrow(.SD)
				p <- sum(abund)/nk
				prob <- 1 - (1-p)^nT
			},by=c('stratum')][,sum(V1)/(nStrat)]
		}
		
		if(!isTRUE(rarefy)){
			innerF_calc(X, nT=tow_target)
		}else{
			rc <- rep(NA_real_, tow_max)
			for(nT in 1:tow_max){
				rc[nT] <- innerF_calc(X, nT=nT)
			}
			out <- rc[tow_target]
			attr(out, "rc") <- rc
			return(out)
		}
		
		
	}else{
	
		innerF <- function(x, nT){
			x[,max(sample(abund, size=nT, replace=TRUE)),by="stratum"][,sum(V1)/(nStrat)]
		}
	
		if(!isTRUE(rarefy)){
			sum(replicate(nPerm, innerF(X, nT=tow_target)))/nPerm
		}else{
			rc <- rep(NA_real_, tow_max)
			for(nT in 1:tow_max){
				rc[nT] <- sum(replicate(nPerm, innerF(X, nT=nT)))/nPerm
			}
			out <- rc[tow_target]
			attr(out, "rc") <- rc
			# plot(attributes(out)$rc, type='o')
			return(out)
		}
	}
	
}


