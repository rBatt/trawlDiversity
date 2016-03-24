#' Rank Species by Temperature
#' 
#' Rank species by the tempeartures at which they are found in the data
#' 
#' @param data A data.table with columns for spp, year, stratum, btemp
#' 
#' @details
#' Defining the thermal preference of species according to the temperatures for which they are observed can be problematic either due to censoring or due to confounding factors. Toward the latter point, if a region is warming and a species is only present in the second half of the time series, it may appear to have a higher temperature preference than a species present in all years. However, the species may have shown up in the data set because of some other variable (e.g., it was an invasive that was recently transported to the region) or because of changes in detectability.
#' 
#' It is the point of detectability that motivates this function. By making the assumption that detectability is constant within a year (within a region), it then becomes 'fair' to compare thermal preferences among species for that year. It is then only problematic to use data for all years. To make use of data for all years we need to reduce the influence of annual differences in temperature on our summary statistic. This function tries to achieve that goal by ranking each species thermal preference by year, then summarizing that rank (as a percentage) among years.
#' 
#' Thus, the summary statistic provided by this function answers the question "When this species is present, what proportion of species were found in colder locations in each year, and what is the average of that proportion among years?"
#' 
#' A problem with this new statistic is that 'thermal rank' is subject to change if the number of warm-water to cold-water species changes among years, which of course can change as regional temperatures change over time. The severity of this problem depends on the breadth of species' thermal tolerances, the spatial variability in temperature, and the temporal variability in temperature. If thermal tolerances are broad and/or if spatial variation is large, and if in addition to either of those temporal variability is relatively low, then this potential bias is minimized.
#' 
#' @return
#' A list of 2 data.table s
#' 
#' @export
rank_temp <- function(data){
	temp_rank_annual <- data[, j={
		syh <- .SD[,list(
			bt_min=min(btemp, na.rm=TRUE), 
			bt_mean=mean(btemp, na.rm=TRUE),
			bt_max=max(btemp, na.rm=TRUE)
		),by=c("spp","year","stratum")]
			
		sy <- .SD[,list(
			bt_min=min(btemp, na.rm=TRUE), 
			bt_mean=mean(btemp, na.rm=TRUE),
			bt_max=max(btemp, na.rm=TRUE)
		),by=c("spp","year")]
	
	
		t_rank <- sy[,list(
			spp = spp,
			bt_min = bt_min,
			bt_mean = bt_mean,
			bt_max = bt_max,
			bt_min_rank = rank(bt_min)/max(rank(bt_min)),
			bt_mean_rank = rank(bt_mean)/max(rank(bt_mean)),
			bt_max_rank = rank(bt_max)/max(rank(bt_max))
		),by="year"]
	
		t_rank
	
	}]

	sr <- data.table:::melt.data.table(temp_rank_annual, id.vars=c("year","spp"))
	sr[!is.finite(value), value:=NA]
	temp_rank <- dcast.data.table(sr[,mean(value, na.rm=TRUE), keyby=c("spp","variable")], formula=spp~variable, value.var="V1")
	
	list(temp_rank_annual=temp_rank_annual, temp_rank=temp_rank)
}