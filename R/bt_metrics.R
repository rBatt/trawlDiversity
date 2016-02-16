#' Bottom Temperature Metrics
#' 
#' Computer bottom temperature and its mean and standard deviation in space and time
#' 
#' @param regs Character vector of region names
#' 
#' @export

bt_metrics <- function(regs){
	
	requireNamespace("zoo", quietly = TRUE)
	
	bt_all0 <- list()
	for(r in 1:length(regs)){
	
		t_reg <- regs[r]
	
		data_in_all0 <- trim_msom(t_reg, gridSize=0.5, depthStratum=100, tolFraction=0.15, grid_stratum=TRUE, plot=FALSE)
	
		nSite_min <- data_in_all0[,lu(stratum), by="year"][,min(V1)]

		data_in_all <- data_in_all0
		setkey(data_in_all, year, stratum, haulid, spp)
	
		haul_bt <- data_in_all[,list(haul_bt=mean(btemp, na.rm=TRUE)), by=c("reg","year","stratum","haulid")]
		strat_bt <- haul_bt[,list(strat_bt=mean(haul_bt, na.rm=TRUE)), by=c("reg","year","stratum")]
		reg_bt <- strat_bt[,list(strat_mean=mean(strat_bt, na.rm=TRUE), strat_sd=sd(strat_bt, na.rm=TRUE)), by=c("reg","year")]
	
		reg_bt <- reg_bt[CJ(reg=una(reg), year=min(year):max(year))]
	
		roll_win <- function(x, stat=c("mean","sd"), width=c("6","9"), ...){
			stat <- match.fun(stat)
			width <- as.integer(match.arg(width))
		
			zoo::rollapplyr(x, by=1, width=width, fill=NA, FUN=stat, na.rm=TRUE)
		}
	
		reg_bt[,strat_mean_6yr_mean:=roll_win(strat_mean, stat="mean", width="6")]
		reg_bt[,strat_mean_9yr_mean:=roll_win(strat_mean, stat="mean", width="9")]
	
		reg_bt[,strat_mean_6yr_sd:=roll_win(strat_mean, stat="sd", width="6")]
		reg_bt[,strat_mean_9yr_sd:=roll_win(strat_mean, stat="sd", width="9")]
	
		bt_all0[[r]] <- reg_bt
	
	}

	bt_all <- rbindlist(bt_all0)
	bt_all <- bt_all[!is.na(strat_mean)]
	
	return(bt_all)
}

