#' Find and Classify Unbroken Stretches of an Event
#' 
#' Like, run length encoding (see \code{\link{rle}}), this function finds stretches of a vector that contain consecutive repeats of a particular value. This function then procedes by classifying each event in a stretch as being closer to the beginning or end of the stretch.
#' 
#' @param X a data.table with columns for \code{present}, \code{year}, \code{now_ext}, \code{col}, and \code{ext_dist_sign} (see Details)
#' 
#' @details The contents of the columns in \code{X} are as follows:
#' \tabular{ll}{
#' 	\code{present} \tab 1 for event, 0 for non-event. "Streaks" refer to consecutive 1's \cr
#' 	\code{year} \tab calendar year, or any time index indicating distance (space, time, etc) between events. \cr
#' 	\code{now_ext} \tab "now extinct"; integer vector of 1 or 0. 1 indicates that \code{present} contains the first 0 following a streak of 1's \cr
#' 	\code{col} \tab "colonization"; integer vector of 1 or 0. 1 indicates that \code{present} contains the first 1 at the beginning of a streak of 1's \cr
#' 	\code{ext_dist_sign} \tab "signed distance to extinction"; output from \code{\link{event_distance}}; magntiude is distance in \code{year}s to non-event, with negative numbers indicating that closest non-event is in past (post-colonization) \cr
#' }
#' 
#' @return Returns a named list of length 2: \code{stretch_id} and \code{hybrid_part}. 
#' 
#' The values in \code{stretch_id} indicate the following:  
#' \tabular{ll}{
#' 	-1 \tab the stretch is an extinction-only stretch, which can only happen for a stretch at the beginning of a vector b/c the first element in \code{present} is 1, so there is no initial 0-to-1 transition \cr
#' 	-2 \tab colonization-only stretch; only occurs at end of vector where last value in \code{present} is 1 (no extinction) \cr
#' 	1,2,... \tab a hybrid stretch, occurs for a stretch in \code{present} that is begun and terminated by a 0; this is the most generic case; the integer value serves to distinguish among multiple hybrid stretches \cr
#' }
#' 
#' When the stretch is classified as a "hybrid", some of the 1's will be closer to the initial 0 than to the terminating 0. Those closer to the inital 0 are considered "post-colonization", whereas those closer to the final 0 are considered "pre-extinction". Proximiting is defined using \code{year}. Ties are ~arbitrarily attributed to post-coloniztion.
#' 
#' The values in \code{hybrid_part} indicate the following:
#' \tabular{ll}{
#' 	1 \tab post-colonization \cr
#' 	2 \tab pre-extinction \cr
#' }
#' 
#' @export
event_stretches <- function(X){
	# must have columns for present, year, now_ext, col, ext_dist_sign
	# returns columns of stretch_id, hybrid_part
	
	present <- X[, present]
	year <- X[, year]
	now_ext <- X[, now_ext]
	col <- X[, col]
	ext_dist_sign <- X[, ext_dist_sign]

	# For convenience, define rle of present
	rls0 <- rle(present)
	
	# Output components, empty
	# stretch_id = integer vector identifying the type and ID of stretch to which each element belongs
	stretch_id <- rep(0, length(present)) # 0 means no stretch
	hybrid_part <- rep(0, length(present))

	# For any stretches to exist, the full series will have all(c(0,1)%in%una(present)) 
	has_stretches <- all(c(0,1)%in%una(present))
	
	# If the current species does not have both presences and absenences, just reutrn 0's
	if(!has_stretches){
		return(list(stretch_id=stretch_id, hybrid_part=hybrid_part))
	}

	# If it has stretches, the series will be comprised of stretches taking on one or more of the following forms:
		# extinction-only
		# colonization-only
		# hybrid stretch

	# Define start and stop of e-only and c-only stretches
	has_col <- any(col==1)
	has_now_ext <- any(now_ext==1)
	e_only_start <- (year==min(year) & present == 1) # extinction-only
	e_only_end <- (now_ext==1 & year==min(year[now_ext==1])) # extinction-only
	c_only_start <- (col==1 & year==max(year[col==1])) # colonization-only
	c_only_end <- (year==max(year) & present==1) # colonization-only

	# A e-only or c-only stretch exists if (any(start_logic) & any(end_logic))
	has_eo <- any(e_only_start) & any(e_only_end)
	has_co <- any(c_only_start) & any(c_only_end)

	# The existance and number of hybrid stretches in a series is defined by:
	n_hybrid <- sum(head(rls0$values[-1],-1)==1)

	# ID for extinction-only stretch
	if(has_eo){
		ind <- seq_along(present)
		eo_ind <- ind[e_only_start] : ind[e_only_end]
		stretch_id[eo_ind] <- -1
	}

	# ID for colonization-only stretch
	if(has_co){
		ind <- seq_along(present)
		co_ind <- ind[c_only_start] : ind[c_only_end]
		stretch_id[co_ind] <- -2
	}

	# ID for hybrid stretches
	if(n_hybrid>0){
		rls <- rls0
		rls$values[c(1,length(rls$values))] <- 0
		rls$values <- rls$values * cumsum(rls$values)
		after_hybrid <- c(0, diff(rls$values)) < 0
		rls$lengths <- rls$lengths - as.integer(after_hybrid) # removes the now_ext==1 from non-hybrid classification
		rls$lengths <- rls$lengths + as.integer(rls$values > 0) # adds the now_ext==1 to the hybrid classification
		hybrid_id <- inverse.rle(rls)
		stretch_id[hybrid_id!=0] <- hybrid_id[hybrid_id!=0]
	
		# A single hybrid stretch can be further subdivided into two parts:
			# post-colonization: an initial portion of the stretch that is closer to its start than end
			# pre-extinction: a latter portion of the stretch that is closer to its end than start
		# Categorization as post-c or pre-e corresponds to whether the nearest absence is in the past or future:
			# post-c: ext_dist_sign<0
			# pre-e: ext_dist_sign>=0
		hybrid_part[stretch_id>0 & ext_dist_sign<0] <- 1 # post-colonization (first part)
		hybrid_part[stretch_id>0 & ext_dist_sign>=0] <- 2 # pre-extinction (second part)
	}
	out <- list(stretch_id=stretch_id, hybrid_part=hybrid_part)
	
	
	return(out)
	
}
