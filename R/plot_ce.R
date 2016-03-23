#' Plot a Species that Colonizes and/or goes Extinct
#' 
#' Plots 3 panels for a given species: its image and time series of detectability, the location of colonizations, and time series of biomass
#' 
#' @param t_sco a character indicating a species name
#' 
#' @details Used in process_msom_figures.R
#' 
#' @returns nothing
#' 
#' @export
plot_ce <- function(t_sco, ...){
	
	omar <- par("mar")
	t_detect <- ab[par=="beta" & (spp)==t_sco]
	par(mar=omar+c(0,0,2,0))
	plot_ab(X=t_detect, t_spp=t_sco, TRUE, ylab="beta", xlab="", ...)
	
	par(mar=omar)
	t_wtcpue <- rd[spp==t_sco, list(wtcpue=mean(wtcpue)), by=c("stratum","year")][,list(wtcpue=mean(wtcpue)), keyby="year"]
	t_wtcpue <- merge(data.table(year=rd[,una(year)]), t_wtcpue, by="year", all=TRUE)
	t_wtcpue[is.na(wtcpue), wtcpue:=0]
	t_wtcpue[,plot(year, wtcpue, type="l", xlab="", ...)]
	t_wtcpue[,points(year, wtcpue, pch=20, col=c("blue","red")[(1+as.integer(wtcpue>0))])]
	
	par(mar=omar)
	t_col_strat <- colonization$col_dt[spp==t_sco & (col_logic), list(year, stratum, col_logic)]
	col_yrs <- t_col_strat[,una(year, na.rm=TRUE)]
	t_col_strat <- t_col_strat[,list(col_logic=sum(col_logic)), by=c("stratum")]
	t_sll <- colonization$n_spp_col_weighted_tot[,list(stratum,lon, lat)]
	t_col_strat <- merge(t_sll, t_col_strat, by=c("stratum"), all=TRUE)
	t_col_strat[is.na(col_logic), col_logic:=0]
	t_col_strat[,plot_space(lon, lat, as.integer(col_logic), TRUE, pch=20, ylab="", xlab="", ...)]
	map(add=TRUE, fill=TRUE, col="white")
	mtext(paste(col_yrs, collapse=" "), side=3, line=0, adj=0.01, font=2)
}