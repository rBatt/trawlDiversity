#' Get Colonizers
#' 
#' Identify species colonizing, going extinct, and persisting between years
#' 
#' @param d An input data.table with columns spp, year (may have others, but they aren't used)
#' 
#' @return
#' A list of length 5:
#' \tabular{ll}{
	#'\code{col_dt} \tab colonization events for each year-stratum-species. \cr
	#' \code{col_ext_dt} \tab colonization (1), extinction (-1), or no change for each species-year \cr
	#' \code{n_cep} \tab the number of colonizations, extinctions, or no-change events in the whole region for each year \cr
	#' \code{n_spp_col_weighted} \tab appropriate aggregation of the column by this name in \code{col_dt} \cr
	#'\code{n_spp_col_weighted_tot} \tab like \code{n_spp_col_weighted}, but summed up across years (and divided by \code{yrs_sampled}) \cr
#'	}
#' 
#' The columns in \code{col_dt} are:
#' \tabular{ll}{
#' \code{col_logic} \tab indicating whehterh or not the colonization or extinction event occurred, is specific to year-stratum-spp \cr
#' \code{n_spp_col} \tab number of species that colonized a stratum in a year \cr
#' \code{n_strat_col} \tab number of strata colonized by this specie sin this year \cr
#' \code{n_spp_col_weighted} \tab Like \code{n_strat_col}, but divided by \code{n_strat_col} (e.g. if a species colonizes 5 strata, it only contributes 0.2 to each stratum). \cr
#' \code{yrs_sampled} \tab The number of years a stratum was sampled. \cr
#' }
#' 
#' @export
get_colonizers <- function(d){
	d <- copy(d)
	
	# ---- baseline steps for getting colonization/ extinction ----
	tbl <- d[,table(spp, year)] # counts (table) of spp in each year
	col_ext <- apply(tbl>0, 1, diff) # for each species, get the presence/absence (>0 logic) sequential difference (col/ext)
	col_ext_dt <- data.table(reshape2::melt(col_ext, varnames=c("year","spp")), key=c("year","spp")) # turn the table into a data.table
	col_ext_dt[,spp:=as.character(spp)] # make spp a character, not factor
	
	# ---- make extinction events occur in the year before the species is absent ----
	# done so that extinctions can be associated with a (last known) place (stratum)
	years <- d[,sort(una(year))] # all the unique years, in order
	ext_year <- years[match(col_ext_dt[value==-1, year], years) - 1] # the year vector for extinctions
	ext_spp <- col_ext_dt[value==-1, spp] # the species vector for extinctions
	col_ext_dt[value==-1, value:=0] # swap the old extinction year to 0
	col_ext_dt <- rbind(data.table(year=years[1], spp=d[,sort(una(spp))], value=0), col_ext_dt) # add 1st yr b/c new ext timing
	col_ext_dt[paste0(spp,year)%in%paste0(ext_spp,ext_year), value:=-1] # change in the extinctions to previous year
	
	# ---- make lon, lat, and depth consistent w/in a stratum-year, instead of just w/in a haulid ----
	# doing this makes things consistent in the long term
	# but without resorting to the gridded lon-lat in the stratum definition, 
	# which would cause strata w/ same lon-lat but different depths to be plotted perfectly on top of each other
	# but at the same time gives a consistent (but unique!) spatial definition for events in x-y coordinates
	d_lld <- d[!duplicated(haulid), list(lon, lat, depth), by=c('stratum', 'year')] # now each haul represented once
	d_lld <- d_lld[, list(lon=meanna(lon), lat=meanna(lat), depth=meanna(depth)), by=c('stratum','year')] # strat-year avg; hauls equal
	d_lld <- d_lld[,list(lon=meanna(lon), lat=meanna(lat), depth=meanna(depth)), by=c('stratum')] # strat avg; years equal
	d[,c('lon','lat','depth'):=NULL]
	d <- merge(d, d_lld, by=c('stratum'))
	
	# ---- Tally up annual colonizations, extinctions, and no-changes ----
	n_cep <- col_ext_dt[,list(n_col=sum(value==1), n_ext=sum(value==-1), n_pers=sum(value==0)), by="year"]
	
	d_new <- merge(d, col_ext_dt, by=c("year","spp"), all=TRUE)
	setnames(d_new, "value", "col_ext")

	get_ce <- function(cet, ce=c('col','ext')){
		ce <- c(col=1, ext=-1)[match.arg(ce)]
		cet <- copy(cet)

		ce_tbl <- cet[,table(year, stratum,spp,col_ext)>0][,,,as.character(ce)]
		ce_dt <- data.table(reshape2::melt(ce_tbl, value.name="col_logic"))
		ce_dt[,c("stratum","spp"):=list(as.character(stratum),as.character(spp))]
		setkey(ce_dt, year, stratum, spp)
		ce_dt[,n_spp_col:=sum(col_logic),by=c("year","stratum")]
		ce_dt[,n_strat_col:=sum(col_logic),by=c("year","spp")]
		ce_dt[n_strat_col>0, n_spp_col_weighted:=sum(as.integer(col_logic)/n_strat_col), by=c("year","stratum")]
		ce_dt[n_strat_col<=0, n_spp_col_weighted:=0]
		
		yrs_sampled <- reshape2::melt(d[,apply(table(stratum,year)>0, 1, sum)], value.name="yrs_sampled")
		yrs_sampled <- data.table(stratum=rownames(yrs_sampled), yrs_sampled=yrs_sampled[,1])
		ce_dt <- merge(ce_dt, yrs_sampled, by=c("stratum"), all=TRUE)
	
		n_spp_ce_weighted <- ce_dt[, list(n_spp_col_weighted=una(n_spp_col_weighted)),by=c("year","stratum","yrs_sampled")]
		n_spp_ce_weighted <- merge(n_spp_ce_weighted, unique(d[,list(lon, lat, depth),keyby=c("year","stratum")]), by=c("year","stratum"))
	
		n_spp_ce_weighted_tot <- n_spp_ce_weighted[,list(lon=mean(lon), lat=mean(lat), depth=mean(depth), n_spp_col_weighted=sum(n_spp_col_weighted)/una(yrs_sampled)),by=c("stratum")]
		
		if(ce == -1){
			setnames(ce_dt, c("n_spp_col", "n_strat_col", "n_spp_col_weighted"),  c("n_spp_ext", "n_strat_ext", "n_spp_ext_weighted"))
			setnames(n_spp_ce_weighted, "n_spp_col_weighted", "n_spp_ext_weighted")
			setnames(n_spp_ce_weighted_tot, "n_spp_col_weighted", "n_spp_ext_weighted")
			
			return(list(ext_dt=ce_dt, n_spp_ext_weighted=n_spp_ce_weighted, n_spp_ext_weighted_tot=n_spp_ce_weighted_tot))
		}else{
			return(list(col_dt=ce_dt, n_spp_col_weighted=n_spp_ce_weighted, n_spp_col_weighted_tot=n_spp_ce_weighted_tot))
		}
	}


	cols <- get_ce(d_new, "col")
	exts <- get_ce(d_new, "ext")
	ce_out <- c(list(col_ext_dt=col_ext_dt, n_cep=n_cep), cols, exts)
	return(ce_out)

}