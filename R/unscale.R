#' Unscale Parameters
#' 
#' Unscale parameters from a polynomial multiple linear regression
#' 
#' @param a1 intercept term
#' @param a2 coefficient for first predictor
#' @param a3 coefficient for first predictor squared
#' @param a4 cofficient for second predictor
#' @param a5 coefficient for second predictor squared
#' @param mu1 mean used to scale first predictor
#' @param mu2 mean used to scale second predictor
#' @param sigma1 standard deviation used to scale first predictor
#' @param sigma2 standard deviation used to scale second predictor
#' 
#' @details
#' Assume a model of the general form \code{y ~ a1 + a2*x1 + a3*x1^2 + a4*x2 + a5*x2^2}. Scaling is assumed to be \code{scale(x1)}, e.g. All arguments can be vectors, for which each element should refer to a different model/ regression.
#' 
#' Created to unscale parameters used in an MSOM. In this case, each iteration on each chain effectively referes to a different set of parameters. So if the posterior was 100 iterations sampled, the 'a' arguments would be of length 100, and mu and sigma arguments could be length 1 (and recycled), or have the same value repeated 100 times. Alterntively, if the 'a' arguments refer to models for which the predictors were scaled differently, just make sure that the mu and sigma values are repeated to match to the corresponding models. No checking is done, but this should be fairly intuitive: the coefficients should match up with the scaling constants used to scale the predictors used in the model that estimated the coefficients.
#' 
#' @return
#' Either a vector (if a1 is of length 1) or a matrix (if a1 is longer than 1) of unscaled parameters.
#' 
#' @export
unscale <- function(a1, a2, a3, a4, a5, mu1, mu2, sigma1, sigma2){
	# i2 <- intercept + sum(-(betaLinear*mu/sigma) + betaQuad*mu^2/sigma^2)
	i2 <- a1 - a2*mu1/sigma1 - a4*mu2/sigma2 + a3*mu1^2/sigma1^2 + a5*mu2^2/sigma2^2
	# bl2 <- betaLinear/sigma - 2*betaQuad*mu/sigma^2
	a22 <- a2/sigma1 - 2*a3*mu1/sigma1^2
	a42 <- a4/sigma2 - 2*a5*mu2/sigma2^2
	# bq2 <- betaQuad/sigma^2
	a32 <- a3/sigma1^2
	a52 <- a5/sigma2^2

	if(length(a1)==1){
		return(c(a1=i2, a2=a22, a3=a32, a4=a42, a5=a52))
	}else{
		return(cbind(a1=i2, a2=a22, a3=a32, a4=a42, a5=a52))
	}
}
