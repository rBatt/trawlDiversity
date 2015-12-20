#' Strip the not-indexed equivalents of the elements whose names imply indexing
#' 
#' An index is implied if a character vector contains the '['. If \code{x} contains both 'w' and 'w[1]' as elements, will return only 'w[1]'. Purpose is for avoiding redundancy when an object name can stand as a shortcut for an object and all of its indices. Care is to be taken, though, because if 'w' is length 2, 'w[2]' would be omitted in the above example, even though by specifying 'w' one may have intended to imply 'w[2]'.
#' 
#' @param x A character vector
#' 
#' @details
#' This function may behave strangely if \code{x} contains the '[' character that does not refer to an index, or is not followed by ']', which together surround regular R-style indexing syntax. See Examples.
#' 
#' @return
#' If \code{x} has positive length, returns a character vector of at least length 1 , but no longer than \code{length(x)}
#' 
#' @examples
#' # Typical/ intended use
#' ifInd_strip_noInd(c("theta","w","w[1]"))
#'
#' # Odd use
#' # Note that the last element is valid R syntax
#' # Yet Z is not removed
#' ifInd_strip_noInd(c("Z", "Z[", "do.call('[', list(Z,1))"))
#' 
#' # However, if normal and odd indexing exists,
#' # only the standard "Z" is stripped out
#' ifInd_strip_noInd(c("Z", "Z[1]", "Z[", "do.call('[', list(Z,1))"))
#' 
#' @export
ifInd_strip_noInd <- function(x){
	strip_ind <- unique(gsub("\\[[0-9,]{1,}\\]$", "", x))
	has_ind_ind <- x[grepl("\\[", x)]
	has_ind_noInd <- unique(gsub("\\[[0-9,]{1,}\\]$", "", 	x[grepl("\\[", x)]))

	x[x%in%strip_ind[!strip_ind%in%has_ind_noInd] | x%in%has_ind_ind]
}