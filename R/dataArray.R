#' Array Format for MSOM Data
#' 
#' Given a data.table of trawl data , turn biological response and covariates into acceptable array format for MSOM analysis
#' 
#' @param X A data.table of trawl data
#' @param bioName A character vector of length one specifying a column in \code{X} to be used as the biological response variable
#' @param covNames Names of columns in \code{X} to be used as covariates
#' 
#' @return
#' A list of length \code{length(c(bioName, covNames))} whose elements are named accordingly, and the contents of each element contain the corresponding arrays. Array dimensions for bioName are \code{stratum~K~spp~year}. Array dimensions for covNames are \code{stratum~K~year}. Grand dimension names for the bioName array are defaults for \code{trawlData::trawlCast}; names for covNames arrays are as in formula.
#' 
#' @export
dataArray <- function(X, bioName="abund", covNames=c("btemp", "nAgg2")){
	requireNamespace("trawlData", quietly = TRUE)
	
	X <- copy(X)

	ebs.c <- trawlCast(X, stratum~K~spp~year, valueName=bioName)
	
	covCast <- function(valueName){
		trawlData::trawlCast(X, 
				stratum~K~year, 
				valueName=valueName, 
				fixAbsent=FALSE, 
				fun.aggregate=meanna, 
				valFill=NA_real_, 
				grandNamesOut=c("stratum","K","year")
			)
	}
	
	covList <- lapply(covNames, covCast)
	
	
	out <- c(list(ebs.c), covList)
	names(out) <- c(bioName, covNames)
	return(out)
	
}