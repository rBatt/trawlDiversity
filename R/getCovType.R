#' Get Covariate Type
#' 
#' Given the output of \code{\link{msomData}}, get the type of covariates (constant or mean of random variable, or sd of random variable).
#' 
#' @param staticDataUV the output of msomData, containing at least U or V element
#' @param staticDataUV_sd the output of msomData, containing at least U or V element, where the covariates have the "_sd" suffix attached to their names, and these covariates correspond to similarly named (sans-suffix) elements of U and V in staticDataUV
#' @param UV character of length 1, indicating whether to treat
#' @param type type of covariate to find
#' 
#' @export
# ---- Split Covariates into Constants and Random Variables ----
getCovType <- function(staticDataUV, staticDataUV_sd, UV=c("U","V"), type=c("constant","mu","sd")){
	
	UV <- match.arg(UV)
	
	# Dimension Sizes, Number, and Names
	dims <- dim(staticDataUV[[UV]])
	nD <- length(dims)
	dn <- dimnames(staticDataUV[[UV]])[[nD]]
	dn_sd <- dimnames(staticDataUV_sd[[UV]])[[nD]]
	
	# Get Names of Desired Covariates
	covNames <- switch(type,
		constant = dn[!dn%in%gsub("_sd","",dn_sd)],
		mu = dn[dn%in%gsub("_sd","",dn_sd)],
		sd = dn_sd
	)
	
	# Get the Desired Covariates
	covOut <- switch(type,
		constant = staticDataUV[[UV]][,,covNames],
		mu = staticDataUV[[UV]][,,covNames],
		sd = staticDataUV_sd[[UV]][,,covNames]
	)
	
	# Convert to Array, Get Dimension Information
	covOut <- as.array(covOut)
	cO_dims <- dim(covOut)
	cO_nD <- length(cO_dims)
	
	# Make Sure covOut Dimensions match Full Cov Dimensions
	if(cO_nD<nD){
		# if here, it's b/c the length of one of the dims is 1
		# (or if covNames is 1, which can be diff from corresp dim b/c of constant vs rv covs) 
		# and in that case, I want an array of the orignal major dims to be returned
		stopifnot(length(covNames)==1 | any(dims==1))
		
		# get a new dims, changing the last dim to be the number of covNames
		dims2 <- dims
		dims2[nD] <- length(covNames)
		
		# add back the original dimensions
		cO_newDim <- c(cO_dims,rep(1,nD-cO_nD)) 
		dimnames2 <- dimnames(staticDataUV[[UV]])
		dimnames2[[nD]] <- covNames
		dn_1 <- dimnames2[dims2==1]
		covOut <- array(covOut, dim=cO_newDim, dimnames=c(dimnames(covOut),dn_1))
		
		# put dimensions back in right order
		new_order <- order(c(which(dims2!=1), which(dims2==1)))
		covOut <- aperm(covOut, new_order)
		
	}
	
	# Return
	return(covOut)
}
