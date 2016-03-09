#' Get Colonizers
#' 
#' Identify species colonizing, going extinct, and persisting between years
#' 
#' @param d An input data.table with columns spp, year (may have others, but they aren't used)
#' 
#' @return
#' A list of length 5:
#' \tabular{ll}{
	#'\code{col_dt} \tab colonization events for each year-stratum-species. \code{col_logic} is specific to year-stratum-spp, \code{n_spp_col} is specific to year-stratum (number of species colonizing), and \code{n_strat_col} is specific to year-spp (number of strata occupied by that species in that year [if it colonized the region that year]). Like the weighted version, \code{n_spp_col_weighted} indicates the number of species colonizing a stratum in a year, but if a species colonizes 5 strata, it only contributes 0.2 to each, instead of 1. The number of years a stratum was sampled is \code{yrs_sampled}.
	#' \code{col_ext_dt} \tab colonization (1), extinction (-1), or no change for each species-year \cr
	#' \code{n_cep} \tab the number of colonizations, extinctions, or no-change events in the whole region for each year \cr
	#' \code{n_spp_col_weighted} \tab appropriate aggregation of the column by this name in \code{col_dt} \cr
	#'\code{n_spp_col_weighted_tot} \tab like \code{n_spp_col_weighted}, but summed up across years (and divided by \code{yrs_sampled}) \cr
#'	}
#' 
#' @export
get_colonizers <- function(d){
	
	tbl <- d[,table(spp, year)]
	col_ext <- apply(tbl>0, 1, diff)
	col_ext_dt <- data.table(reshape2::melt(col_ext, varnames=c("year","spp")), key=c("year","spp"))
	
	n_cep <- col_ext_dt[,list(n_col=.SD[value==1, sum(value)], n_ext=.SD[value==-1, abs(sum(value))], n_pers=.SD[value==0, lu(spp)]), by="year"]
	
	d_new <- merge(d, col_ext_dt, by=c("year","spp"), all=TRUE)
	setnames(d_new, "value", "col_ext")
	
	col_tbl <- d_new[col_ext==1,table(year, stratum,spp)>0]
	col_dt <- data.table(reshape2::melt(col_tbl, value.name="col_logic"))
	col_dt[,n_spp_col:=sum(col_logic),by=c("year","stratum")]
	col_dt[,n_strat_col:=sum(col_logic),by=c("year","spp")]
	col_dt[n_strat_col>0, n_spp_col_weighted:=sum(as.integer(col_logic)/n_strat_col), by=c("year","stratum")]
	col_dt[n_strat_col<=0, n_spp_col_weighted:=0]
	
	yrs_sampled <- reshape2::melt(d[,apply(table(stratum,year)>0, 1, sum)], value.name="yrs_sampled")
	yrs_sampled <- data.table(stratum=rownames(yrs_sampled), yrs_sampled=yrs_sampled[,1])
	col_dt <- merge(col_dt, yrs_sampled, by=c("stratum"), all=TRUE)
	
	n_spp_col_weighted <- col_dt[, list(n_spp_col_weighted=una(n_spp_col_weighted)),by=c("year","stratum","yrs_sampled")]
	n_spp_col_weighted <- merge(n_spp_col_weighted, unique(d[,list(lon, lat, depth),keyby=c("year","stratum")]), by=c("year","stratum"))
	
	n_spp_col_weighted_tot <- n_spp_col_weighted[,list(lon=mean(lon), lat=mean(lat), depth=mean(depth), n_spp_col_weighted=sum(n_spp_col_weighted)/una(yrs_sampled)),by=c("stratum")]
		
	
	return(list(col_dt=col_dt, col_ext_dt=col_ext_dt, n_cep=n_cep, n_spp_col_weighted=n_spp_col_weighted, n_spp_col_weighted_tot=n_spp_col_weighted_tot))

}