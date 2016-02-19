#' Get species that show up and remain
#' 
#' A crude method for determining which species show up in the data set due to enhanced detection
#' 
#' @param X A data.table containing trawl data
#' @param n_yrs_start after this many years, species seen for the first time "show up"; acts as a burn-in
#' @param n_yrs_remain how many years must be remaining in the time series after a species "shows up" to evaluate it as "sticking around"
#' @param n_yrs_miss_after how many years the species can be missing during remaining (post show-up) years and still be considered to "stick around"
#' 
#' @details
#' Function looks for species that appear after the second year, but before the final 5 years, and that in the remaining years (of which there must be at least 5) the species is observed in at least every year but 1. So a species that shows up after many years of sampling, and then is seen every year thereafter, will be returned by this function. Those species might be showing up and sticking around b/c of changing detection by the survey (I'd expect species to come on the scene more slowly, such that after it first shows up, there should still be some years it isn't observed)
#' 
#' @return
#' A data.table containing the region, year of appearance, species, stats for how long and how consistently that species showed up.
#' 
#' @export
get_show_up_spp <- function(X, n_yrs_start=2, n_yrs_remain=5, n_yrs_miss_after=1){
	added_spp <- X[, j={
		old <- .SD[year==una(year)[1],una(spp)]
		new <- list()
		uyr <- .SD[,una(year)]
		for(y in 2:lu(year)){
			td <- .SD[(year)==uyr[y]]
			t_spp <- td[,una(spp)]
			spp_added <- c(NA, setdiff(t_spp, old))
			old <- una(c(old,t_spp))
			new[[y]] <- data.table(year=td[,una(year)], spp=spp_added)

		}

		rbindlist(new)

	}, by="reg"]

	n_added_spp <- added_spp[,list(n_new_spp=lu(spp)-1), by=c("reg","year")]

	t_X <- copy(X)
	setkey(t_X, reg, spp, stratum, haulid)
	t_added_spp <- copy(added_spp[!is.na(spp)])
	setkey(t_added_spp, reg, spp)
	t_added_spp[,year_added:=year]
	t_added_spp[,year:=NULL]

	X_added <- merge(t_X, t_added_spp, by=c("reg","spp"))
	spp_added_freq <- X_added[,j={
		yrs <- sort(una(year))
		uya <- una(year_added)
		.SD[,list(
			year_added=una(year_added),
			total_years = lu(yrs)+1,
			years_elapsed=sum(yrs <= una(year_added))+1,
			years_remaining=sum(yrs > una(year_added)),
			years_seen_after=lu(year)-1
		), by="spp"]

	}, by=c("reg")]

	show_up_spp <- spp_added_freq[years_elapsed>=n_yrs_start & ((years_remaining - years_seen_after)<=n_yrs_miss_after) & years_remaining >=n_yrs_remain]
	
	return(show_up_spp)
	
}