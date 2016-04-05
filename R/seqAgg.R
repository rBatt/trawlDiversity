#' Sequential Aggregation
#' 
#' Sequentially aggregate columns in a data.table, varying the function and 'by' values
#' 
#' @param X a data.table
#' @param val character vector indicating columns in X that should be sequuentially aggregated
#' @param FUN a list of functions
#' @param by a list of character vectors; should be of same length as FUN
#' @param ... arguments to be passed to each element of \code{FUN}
#' 
#' @details
#' The lengths of \code{FUN} and \code{by} should be of the same length, but the length of \code{val} is independent.
#' 
#' Useful when a summary statistic needs to be computed sequentially (possibly taking the mean over some dimensions before computing the standard deviation, e.g.)
#' 
#' @examples
#' neus_dat <- data_all[reg=="neus"&year==2010]
#' a1 <- c("K","stratum")
#' a2 <- c("stratum")
#' a3 <- c("")
#' a_list <- list(a1, a2, a3)
#' f_list <- list(mean, mean, sd)
#' 
#' neus_dat[,sd(btemp,na.rm=TRUE)] # naive
#' seqAgg(neus_dat, val="btemp", FUN=f_list, by=a_list, na.rm=TRUE) # sequential
#' seqAgg(neus_dat, val=c("btemp","depth"), FUN=f_list, by=a_list, na.rm=TRUE) # 2 columns
#' 
#' @export
seqAgg <- function(X, val, FUN, by, ...){
	for(f in 1:length(FUN)){
		t_f <- FUN[[f]]
		t_b <- by[[f]]
		X <- X[,lapply(eval(s2c(val)), t_f, ...), by=t_b]
		setnames(X, paste0("V",seq_along(val)), val)
	}
	return(X)
}

