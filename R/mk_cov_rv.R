#' Make Covariate (and their powers) a RV
#' 
#' Get the mean and SD of a covariate
#' 
#' @param X data.table
#' @param cName name of covariate
#' @param across RV is computed across this dimension
#' @param by RV is computed separately for each of these dimensions
#' @param pow Raise to this power. Default is 2.
#' 
#' @details
#' Adds columns for the covarate, it's power, and both their sd's. Columns for SD's are named by appending "_sd", column for power is named by appending the numeric value of \code{pow}.
#' 
#' @return
#' NULL is retured, but the data.table is modified, with appropriate columns. Original covariate column will be modified.
#' 
#' @examples
#' X <- data.table(a=rnorm(1000), b=sample(letters, 1000, replace=TRUE), key="b")
#' X[,K:=seq_along(a),by="b"]
#' mk_cov_rv_pow(X, "a", across="K", by="b", pow=2)
#' 
#' @import trawlData
#' @export
mk_cov_rv <- function(X, cName, across, by){
	cn_sd <- paste0(cName, "_sd")
	sd2 <- function(x, ...){
		if(sum(!is.na(x))==1){return(0)}else{return(sd(x, ...))}
	}
	X[,c(cn_sd):=sd2(.SD[,list(mean(eval(s2c(cName))[[1]],na.rm=T)),by=c(across)][,V1], na.rm=TRUE), by=c(by)]
	X[,c(cName):=meanna(eval(s2c(cName))[[1]]), by=c(by)]
	
	invisible(NULL)
}
 
#' @describeIn mk_cov_rv Get the mean and SD of a covariate, and the same for the covariate^pow
#' @export
mk_cov_rv_pow <- function(X, cName, across, by, pow=2){
	mk_cov_rv(X, cName, across, by)
	cn_sd <- paste0(cName, "_sd")
	cn_pow <- paste0(cName, pow)
	cn_pow_sd <- paste0(cn_pow, "_sd")
	
	#' Propagate standard deviation from raising to power
	#'
	#' Calculates the standard deviation of a random variable created by multiplying another RV by a constant and raising it to a power
	#'
	#' @param A The random variables being operated on
	#' @param sigma_A Standard deviations of A
	#' @param b The exponent
	#' @param a the multiplication coefficient
	#'
	#' @details
	#' Of the form f = a*A^b. Sigma must be provided. Elements of A are treated as potentially separate RV's.
	prop_sd_pow <- function(A, sigma_A, b=2, a=1){
		f <- a*A^b
		sigma_f <- (f*b*sigma_A)
	}
	
	X[,c(cn_pow):=eval(s2c(cName))[[1]]^pow]
	X[,c(cn_pow_sd):=prop_sd_pow(A=eval(s2c(cName))[[1]], sigma_A=eval(s2c(cn_sd))[[1]], b=pow)]
	
	invisible(NULL)
	
}