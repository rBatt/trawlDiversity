#' Get Data from All Regions
#' 
#' Get data from all regions, trimmed and ready for MSOM
#' 
#' @param region Character indicating region abbreviation
#' 
#' @return a data.table
#' 
#' @export
get_data_all <- function(region){
	poss_regs <- c("ebs", "ai", "goa", "wctri", "wcann", "gmex", "sa", "neus", "shelf", "newf")
	if(missing(region)){region <- poss_regs}
	regs <- match.arg(region, poss_regs, several.ok=TRUE)
	pretty_reg <- c("ebs"="E. Bering Sea", "ai"="Aleutian Islands", "goa"="Gulf of Alaska", "wc"="West Coast US", "gmex"="Gulf of Mexico", "sa"="Southeast US", "neus"="Northeast US", "shelf"="Scotian Shelf", "newf"="Newfoundland")
	pr <- names(pretty_reg)

	# names(pretty_col) <- pr


	# ====================================
	# = Get Trimmed Data for Each Region =
	# ====================================
	data_in_regs <- list()
	for(r in 1:length(regs)){
		t_reg <- regs[r]
		t_trimmed <- trim_msom(t_reg, gridSize=0.5, depthStratum=100, tolFraction=0.15, grid_stratum=TRUE, plot=FALSE)
		setkey(t_trimmed, year, stratum, haulid, spp)
		data_in_regs[[regs[r]]] <- t_trimmed
	}


	data_all <- rbindlist(data_in_regs, fill=TRUE)

	return(data_all)
}