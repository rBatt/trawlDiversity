#' Time Slope
#' 
#' Calculate the temporal trendline (quickly)
#' 
#' @param x temporal predictor variable
#' @param y response variable
#' @param ... others arguments; not used
#' 
#' @details
#' Requires at least 3 non-na values. 
#' 
#' For benchmarking, see: https://github.com/rBatt/trawl/commit/a04d8e166d7093538c698368d48542a0e60bd16c
#' 
#' @export
timeSlope<- function(x, y, ...){
	stopifnot(all(!is.na(x)))
	nona <- !is.na(y)
	
	if(sum(nona)<3){
		NA
	}else{
		ly <- length(y)
		if(ly>=4E3){
			tvec <- x[nona]
			y2 <- y[nona]
			cor(y2, tvec)*(sd(y2)/sd(tvec))

		}else{
			tvec <- cbind(1,x)[nona,] # the 1 is a dummy vector for the intercept
			ttvec <- t(tvec)
			(solve(ttvec%*%tvec)%*%ttvec%*%y[nona])[2]
		}
	}
}