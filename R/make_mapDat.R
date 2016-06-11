#' Make data for maps related to diversity analysis
#' 
#' For each region, get the data appropriate for maps of richness and environmental metrics.
#' 
#' @param p The \code{p} data object from the process_msom_script.R file.
#' 
#' @details
#' The data.table contains summary statistics to collapse time (long-term mean, e.g.), and retains stratum-specific information.
#' 
#' @return a data.table
#' 
#' @export
make_mapDat <- function(p){
	mapDat <- list()
	for(r in 1:length(p)){
		pt1 <- trawlAgg(
			p[[r]]$rd,
			bioFun=meanna, envFun=meanna,
			envCols=c("btemp","depth","stemp","lon","lat"),
			bio_lvl="spp",time_lvl="year",space_lvl="stratum",
			metaCols=c("reg"),meta.action="unique1"	
		)
		pt2 <- pt1[,j={
			avgRich <- .SD[,lu(spp),by="time_lvl"][,meanna(V1)]
			sdRich <- .SD[,lu(spp),by="time_lvl"][,sd(V1, na.rm=TRUE)]
			avgBtemp <- .SD[,meanna(btemp),by="time_lvl"][,meanna(V1)]
			sdBtemp <- .SD[,meanna(btemp),by="time_lvl"][,sd(V1, na.rm=TRUE)]
			data.table(avgRich=avgRich, sdRich=sdRich, avgBtemp=avgBtemp, sdBtemp=sdBtemp)
		},by=c("reg","stratum")]
		to_merge <- c(p[[r]]$colonization[c("n_spp_col_weighted_tot","n_spp_ext_weighted_tot","n_cep")])
		mapDat[[r]] <- merge(to_merge[["n_spp_col_weighted_tot"]], to_merge[["n_spp_ext_weighted_tot"]], by=c("stratum","lon","lat","depth"),all=TRUE)
		mapDat[[r]] <- merge(mapDat[[r]], pt2, by="stratum",all=TRUE)
	}
	mapDat <- rbindlist(mapDat)[reg!="wcann"]
	# mapDim <- mapDat[,list(r_lon=diff(range(lon)),r_lat=diff(range(lat))),by="reg"]
	# mapDim[,c("lon_scale","lat_scale"):=list(r_lon/min(r_lon), r_lat/min(r_lat))]
	# mapDim[,ll_ratio:=r_lon/r_lat]
	
	return(mapDat)
}