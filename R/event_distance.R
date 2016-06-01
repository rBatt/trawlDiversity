#' Event Distance
#' 
#' Distance from each event to nearest non-event
#' 
#' @param x vector of values indicating events (and arbitrary values indicating non-events)
#' @param positions indicating locations separating distance between values in \code{x}; see example use in \code{\link{event_stretches}} and its use of \code{year}; default assumes event spacing
#' @param event_value value in \code{x} inidicating an event
#' @param keep_sign Logical; should distances retain their sign? See Note
#' 
#' @details Finds whether each value in a streak (run) is closer to the beginning or the end of the streak, and then finds the distance to whichever is nearest.
#' 
#' @note When \code{keep_sign} is TRUE, the sign of values equidistant from 2 non-events may be misleading: even though the event is not closer to either, the sign will be negative.
#' 
#' @export
event_distance <- function(x, positions=seq_along(x), event_value=1, keep_sign=FALSE){
	# distance of each event to nearest non-event
	
	event_index <- which(x == event_value)
	nonEvent_index <- seq_along(x)[-event_index]
	
	if(length(event_index) == 0 | length(nonEvent_index) == 0){
		return(as(rep(NA, length(x)), class(x)))
	}
	
	event_position <- positions[event_index]
	nonEvent_position <- positions[nonEvent_index]
	
	# dists <- abs(outer(event_position, nonEvent_position, "-"))
	# event_dists <- rep(0, length(x))
	# event_dists[event_index] <- apply(dists, 1, min)
	
	dists_s <- outer(nonEvent_position, event_position, "-")
	dists <- abs(dists_s)
	event_dists <- rep(0, length(x))
	signed_dists <- apply(dists_s, 2, function(x)x[which.min(abs(x))])
	if(keep_sign){
		event_dists[event_index] <- signed_dists
	}else{
		event_dists[event_index] <- abs(signed_dists)
	}
	
	return(as(event_dists, class(x)))
}