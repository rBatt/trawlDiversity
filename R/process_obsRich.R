#' Process Observed Data for Richness
#' 
#' Used observations, as opposed to model output, to summarize elements related to regional richness
#' 
#' @param X a data.table containing data for a region. Should \code{\link{data_all}} (1 year of which used as input to \code{\link{run_msom}})
#' @param msom_yrs Optional vector of integers indicating calendar years analyzed by the msom, and that will be used to subset \code{X}
#' 
#' @details 
#' The analysis performed by this function was previously contained in process_msomStatic.
#' 
#' @export
process_obsRich <- function(X, msom_yrs){
	requireNamespace("rbLib", quietly=TRUE)
	
	# rd
	if(missing(msom_yrs)){
		rd <- X
	}else{
		rd <- X[year%in%msom_yrs]
	}
	
	# rd_yr
	rd_yr <- rd[,sort(una(year))]
	
	# naive_rich
	naive_rich <- rd[,list(naive_rich=lu(spp)), keyby=c("year")]
	
	# colonization
	colonization <- get_colonizers(d=rd)
	
	# bt
	vars <- c("lon","lat","btemp","depth")
	funs <- list(mean,mean)
	bys <- list(c("haulid","stratum","year"),c("stratum","year"))
	bt_depth <- seqAgg(rd, val=vars, FUN=funs, by=bys, na.rm=TRUE)
	setnames(bt_depth, c("btemp"), c("bt"))
	bt_depth[,c("bt_col","depth_col"):=list(rbLib::zCol(256, bt), rbLib::zCol(256, depth))]
	bt <- bt_depth
	
	
	output <- list(rd_yr=rd_yr, rd=rd, naive_rich=naive_rich, colonization=colonization, bt=bt)
	
	return(output)
	
}
