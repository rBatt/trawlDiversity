#' Streak Length
#' 
#' Finds streaks in a vector containing consecutively-repeated values, and in each streak indicates the number of values into the streak so far
#' 
#' @param x vector of numbers
#' @param streak_val value in \code{x} considered as part of streak
#' @param fill_val value in output to indicate indices that are part of non-streak
#' 
#' @details
#' Functions like \code{\link{rle}}, but within each run (streak) indicates the current position of each value where the first value in the streak is 1 and the last value in the streak is the length of the streak
#' 
#' @return Returns a vector the same length as \code{x} that contains streak lengths and 0's
#' 
#' @examples
#' streak_length(c(1,1,1,0,0,1,0,1,0,0,1,1))
#' 
#' @export
streak_length <- function(x, streak_val=1, fill_val=0){
	rl <- rle(x)
	rlv <- rl$values
	streaks <- sapply(rl$lengths[rlv==streak_val & is.finite(rlv)], seq_len)	
	
	finite_x <- is.finite(x)
	x[x!=streak_val | !finite_x] <- fill_val
	x[x==streak_val & finite_x] <- Reduce(c, streaks)
	
	return(x)
}