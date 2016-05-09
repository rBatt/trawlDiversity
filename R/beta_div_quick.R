#' Calculate Beta Diversity Quickly
#' 
#' Calculates beta diversity using a modified version of Legendre & Caceres's \code{beta.div}
#' 
#' @param Y a matrix of presence-absence data; species are columns, sites are rows
#' @param method the method for calculating beta diversity
#' 
#' @details in Ecology Letters (2013) 16: 951-963 "Beta diversity as the variance of community data: dissimilarity coefficients and partitioning" a function \code{beta.div} is provided in the supplement. This function is an abbreviated version that can be used for 4 out of the original's 20 options, and only returns the total beta diversity (original returned many other values). This function was written to be much faster than the original, as the original can be quite slow if being run for 10's of thousands of matrices (think Bayesian posterior of MSOM results).
#' 
#' @return Returns a single numeric value that is the total beta diversity of the matrix (SS_total/(n-1))
#' 
#' @seealso \link{process_msomStatic}
#' 
#' @examples
#' Y <- matrix(rbinom(18, 1, 0.5), ncol=3)
#' bd1 <- beta_div_quick(Y, "hellinger")
#' bd2 <- beta_div_quick(Y, "jaccard")
#' 
#' @export
beta_div_quick <- function(Y, method=c("hellinger", "jaccard","sorensen","ochiai")){
	method <- match.arg(method)
	
	BD.group1 <- function(Y, method){
		# if(method=="hellinger") Y = decostand(Y, "hellinger")
		Y <- vegan::decostand(Y, "hellinger")
		s <- scale(Y, center=TRUE, scale=FALSE)^2   # eq. 1
		SStotal <- sum(s)          # eq. 2
		BDtotal <- SStotal/(n-1)   # eq. 3
		BDtotal
	}

	BD.group2 <- function(Y, method){
		D <- switch(method,
			jaccard = ade4::dist.binary(Y, method=1),
			sorensen = ade4::dist.binary(Y, method=5),
			ochiai = ade4::dist.binary(Y, method=7)
		)
		SStotal <- sum(D^2)/n      # eq. 8
		BDtotal <- SStotal/(n-1)   # eq. 3
		BDtotal
	}
	
	n <- nrow(Y)
	
	if(method %in% c("jaccard","sorensen","ochiai")){
		requireNamespace("ade4", quietly = TRUE)
		res <- BD.group2(Y, method)
	}else{
		requireNamespace("vegan", quietly = TRUE)
		res <- BD.group1(Y, method)
	}

	return(res)
}