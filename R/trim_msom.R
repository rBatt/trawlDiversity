#' Trim trawl data for msom
#' 
#' Performs the \code{trawlTrim} function the clean region data, then trims down to species (not genera) that have been observed at least 10 times. Does aggregating within a haul/ spp. Adds an abundance column, which is actually just 0 if wtcpue == 0, and 1 if wtcpue > 0.
#' 
#' @param reg name of region to be passed to \code{\link{trawlTrim}}
#' @param gridSize grid size to be passed to \code{\link{ll2strat}}
#' @param plot Logical, default FALSE; for \code{check_strat}, should the stratum tolerance be plotted?
#' 
#' @return
#' A data.table 
#' 
#' @export
trim_msom <- function(reg, gridSize=1, plot=FALSE){
	
	# use default trimming
	X.t <- trawlTrim(reg)
	
	# define stratum on 1º grid
	X.t[,stratum:=ll2strat(lon, lat, gridSize=gridSize)]
	
	
	# drop strata that weren't sampled every year
	check_strat(X.t, reg, gridSize=gridSize, strat_tol=0, append_keep_strat=TRUE, plot=plot)
	X.t <- X.t[(keep_strat)]
	X.t[,keep_strat:=NULL]
	
	# get the names of the species that were observed at least 10 times
	lots_spp <- X.t[,rowSums(table(spp,haulid))>=10]
	names_lots_spp <- names(lots_spp[lots_spp])
	
	# drop taxa that aren't species or weren't observed at least 10 times
	X.t <- X.t[(spp%in%names_lots_spp) & (taxLvl=="species" | taxLvl=="subspecies")]
	X.t[,lu(spp)]

	# aggregate to make sure a haul
	# doesn't have duplicates of individuals
	# this aggregation uses the biological sum
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
	X.agg1[nAgg!=1,j={par(mar=c(0,0,3,0),mfrow=rbLib::auto.mfrow(lu(spp))); for(i in 1:lu(spp)){trawlData::sppImg(unique(spp)[i],unique(common)[i])}}]

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