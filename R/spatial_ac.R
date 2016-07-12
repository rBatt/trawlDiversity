#' Compute Local Spatial Autocorrelation
#' 
#' Sites don't need to be gridded or random; defines neighbors and calculates local spatial autocorrelation as Moran's I
#' 
#' @param lon longitude
#' @param lat latitude
#' @param value value of interest for which to calculate AC
#' 
#' @details
#' The neighborhood is defined in a 2-step process. First, the nearest neighbor for each location is found, and the maximum distance among nearest neighbors recorded. Second, the neighborhood for a focal point is then defined as all points within that distance. So it's defined as all points within a raidus with a size that is the minimum distance required for all points to have at least 1 neighbor.
#' 
#' Local Moran's I is calculated for each site. P value is adjusted to number of sites in each "region" (not the same 'region' used in the trawlDiversity package); see \code{spdep::localmoran}.
#' 
#' @return
#' A named list of length 3. The first element of this list, \code{max2NDist} is the minimum distance and is a length 1 numeric vector. The second element, \code{nb} is the neighborhood object returned by \code{spdep::knn2nb}. The third element, \code{I}, is a data.table that contains columns as returned by \code{spdep::localmoran} with a few adjustments: the column for the pvalue is changed to \code{lI_pvalue} and columns for \code{lon} and \code{lat} of the sites are added.
#' 
#' @examples
#' # calculate local I on random subset of volcano
#' data(volcano)
#' lon <- rep(seq_len(ncol(volcano)), each=nrow(volcano))
#' lat <- rep(seq_len(nrow(volcano)), ncol(volcano))
#' rind <- sample(seq_along(lon), 100) # random subset
#' value <- c(volcano)
#' spac <- spatial_ac(lon[rind], lat[rind], value[rind])
#' 
#' # plot neighborhood
#' plot(spac[["nb"]], cbind(lon[rind], lat[rind]))
#' 
#' # plot significant I's
#' sig_ind <- spac[["I"]][,lI_pvalue] < 0.05
#' plot(spac[["nb"]], cbind(lon[rind], lat[rind]))
#' spac[["I"]][sig_ind, points(lon, lat, col=zCol(256, Ii), pch=19)]
#' mapLegend(zlim=spac[["I"]][sig_ind,range(Ii)], cols=rbLib::zCol(6, 1:6))
#' 
#' @export
spatial_ac <- function(lon, lat, value){
	requireNamespace("spdep", quietly = TRUE)
	
	localAC <- list()
	X <- data.table(lon, lat, value)
	locs_ll <- as.matrix(X[,list(lon,lat)])
	locs <- simplify2array(ll2km(X[,lon], X[,lat]))[,2:1]
	
	nn1 <- spdep::knn2nb(spdep::knearneigh(locs, k=1))
	max2NDist <- max(unlist(spdep::nbdists(nn1, locs)))
	localAC$max2NDist <- max2NDist
	localAC$nb <- spdep::dnearneigh(locs, d1=0, d2=max2NDist) # graph2nb(gabrielneigh(locs))
	ml_out <- spdep::localmoran(X[,value], spdep::nb2listw(localAC$nb), p.adjust.method="BH")
	localAC$I <- data.table(locs_ll, ml_out)
	setnames(localAC$I, "Pr(z > 0)", "lI_pvalue")
	
	return(localAC)
}