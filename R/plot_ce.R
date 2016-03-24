#' Plot a Species that Colonizes and/or goes Extinct
#' 
#' Plots 3 panels for a given species: its image and time series of detectability, the location of colonizations, and time series of biomass
#' 
#' @param t_sco a character indicating a species name
#' @param pad_top_mar numeric indicating the amount to increase the size of the top margin for the first of the 3 panels
#' @param plt_pts logical indicating whether or not to plot points for first figure (see \code{plot_ab})
#' @param ... arguments passed to \code{plot_ab}, \code{plot}, and \code{plot_space}
#' 
#' @details Used in process_msom_figures.R
#' 
#' @returns nothing
#' 
#' @export
plot_ce <- function(t_sco, pad_top_mar=2, plt_pts=FALSE, ...){
	
	omar <- par("mar")
	t_detect <- ab[par=="beta" & (spp)==t_sco]
	par(mar=omar+c(0,0,pad_top_mar,0))
	plot_ab(X=t_detect, t_spp=t_sco, plt_img=TRUE, plt_pts=plt_pts, ylab="beta", xlab="",  ...)
	
	par(mar=omar)
	t_wtcpue <- rd[spp==t_sco, list(wtcpue=mean(wtcpue)), by=c("stratum","year")][,list(wtcpue=mean(wtcpue)), keyby="year"]
	t_wtcpue <- merge(data.table(year=rd[,una(year)]), t_wtcpue, by="year", all=TRUE)
	t_wtcpue[is.na(wtcpue), wtcpue:=0]
	t_wtcpue[,plot(year, wtcpue, type="l", xlab="", ...)]
	t_wtcpue[,points(year, wtcpue, pch=20, col=c("blue","red")[(1+as.integer(wtcpue>0))])]
	
	par(mar=omar)
	t_col_strat <- colonization$col_dt[spp==t_sco & (col_logic), list(year, stratum, col_logic)]
	col_yrs <- t_col_strat[,sort(una(year, na.rm=TRUE))]
	t_col_strat <- t_col_strat[,list(col_logic=sum(col_logic)), by=c("stratum")]
	t_sll <- colonization$n_spp_col_weighted_tot[,list(stratum,lon, lat)]
	t_col_strat <- merge(t_sll, t_col_strat, by=c("stratum"), all=TRUE)
	t_col_strat[is.na(col_logic), col_logic:=0]
	t_col_strat[,plot_space(lon, lat, as.integer(col_logic), TRUE, pch=20, ylab="", xlab="", ...)]
	map(add=TRUE, fill=TRUE, col="white")
	mtext(paste(col_yrs, collapse=" "), side=3, line=0, adj=0.01, font=2)
}