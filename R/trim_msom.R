#' Trim trawl data for msom
#' 
#' Performs the \code{trawlTrim} function the clean region data, then trims down to species (not genera) that have been observed at least 10 times. Does aggregating within a haul/ spp using \code{trawlAgg}. Adds an abundance column, which is actually just 0 if wtcpue == 0, and 1 if wtcpue > 0.
#' 
#' @param reg Character name of region or data.table to be passed to \code{\link{trawlTrim}}.
#' @param gridSize numeric grid size to be passed to \code{\link{ll2strat}}
#' @param grid_stratum logical, default TRUE; whether or not the strata should be defined on a grid, or used original stratum definition
#' @param depthStratum The depth range to use to define stratum; the default does not use depth at all. If using, a reasonable value might be 100.
#' @param tolFraction The fraction of years that a stratum can not be sampled yet still be included in output. Default is 1/3, indicating that a stratum will be included in output so long as it is sampled for at least 2/3 of years. Actual number of years tolerated is rounded up. This value determines the strat_tol argument in \code{\link{check_strat}}
#' @param plot Logical, default FALSE; for \code{check_strat}, should the stratum tolerance be plotted?
#' 
#' @return
#' A data.table 
#' 
#' @export
trim_msom <- function(reg, gridSize=1, grid_stratum=TRUE, depthStratum=NULL, tolFraction=1/3, plot=FALSE){
	
	# use default trimming
	X.t <- trawlTrim(reg)
	
	
	if(grid_stratum){
		# define stratum on gridSizeº grid
		X.t[,stratum:=ll2strat(lon, lat, gridSize=gridSize)]
		if(!is.null(depthStratum)){
			X.t[,stratum:=paste(stratum, roundGrid(depth, depthStratum))]
		}
	}
	
	
	# drop strata that weren't sampled every year
	if(reg == "gmex" | reg == "neus"){
		X.t <- X.t[(year)!=2015,]
	}
	
	strat_tol <- X.t[,ceiling(lu(year)*tolFraction)]
	check_strat(X.t, reg, gridSize=gridSize, strat_tol=strat_tol, append_keep_strat=TRUE, plot=plot)
	
	X.t <- X.t[(keep_strat)]
	X.t[,keep_strat:=NULL]
	
	# get the names of the species that were observed at least 10 times
	lots_spp <- X.t[,apply(table(spp, stratum,year)>1, 1, sum)]>=10
	names_lots_spp <- names(lots_spp[lots_spp])
	
	# drop taxa that aren't species or weren't observed at least 10 times
	X.t <- X.t[(spp%in%names_lots_spp) & (taxLvl=="species" | taxLvl=="subspecies")]
	# X.t[,lu(spp)]

	# aggregate to make sure a haul
	# doesn't have duplicates of individuals
	# this aggregation uses the biological sum
	
	if(reg== 'neus'){
		X.t <- trawlAgg(
			X = X.t,
			bioFun = una,
			envFun = meanna,
			bio_lvl = "ref",
			space_lvl = "haulid",
			time_lvl = "haulid",
			bioCols= c("wtcpue"),
			envCols = c("stemp","btemp","depth","lon","lat"),
			metaCols = c("reg","datetime","season","year","lon","lat","stratum","common","spp","species"),
			meta.action = "unique1"
		)
	}
	
	X.agg1 <- trawlAgg(
		X = X.t,
		bioFun = sumna,
		envFun = meanna,
		bio_lvl = "spp",
		space_lvl = "haulid",
		time_lvl = "haulid",
		bioCols= c("wtcpue"),
		envCols = c("stemp","btemp","depth","lon","lat"),
		metaCols = c("reg","datetime","season","year","lon","lat","stratum","common","species"),
		meta.action = "unique1"
	)
	# X.agg1[nAgg!=1,j={par(mar=c(0,0,3,0),mfrow=rbLib::auto.mfrow(lu(spp))); for(i in 1:lu(spp)){trawlData::sppImg(unique(spp)[i],unique(common)[i])}}]

	# drop new "time" column
	# change name of agg to indicate step
	# add a 1 and 1/10 degree grid
	X.agg1[,time_lvl:=NULL]
	setnames(X.agg1, "nAgg", "nAgg1")
	# X.agg1[,aggStrat:=ll2strat(lon, lat, gridSize=0.1)]

	# ========================
	# = Give Abundance Value =
	# ========================
	X.agg1[,abund:=as.integer(wtcpue>0)]
	
	
	# ============
	# = Define K =
	# ============
	X.agg1[,K:=as.integer(as.factor(haulid)), by=c("reg","year","stratum")]
	X.agg1[,Kmax:=max(K), by=c("reg","year","stratum")]
	
	return(X.agg1)
	
}