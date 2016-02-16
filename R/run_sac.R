#' Run Species Accumulation Curve Estimate of Richness
#' 
#' Calculate richness using a species accumulation curve
#' 
#' @param regs Character vector of region names
#' @param data_regs Option list of data.tables; list should be named according to \code{regs
#' 
#' @export

run_sac <- function(regs, data_regs){

	requireNamespace("reshape2", quietly = TRUE)
	requireNamespace("vegan", quietly = TRUE)
	
	sac_fun <- function(td, nSite_rich=20){
		pres <- reshape2::acast(t_data, stratum~spp, value.var='wtcpue', fun.aggregate=length)
		w <- t_data[,list(nK=lu(haulid)), by="stratum"][,nK]
		sac <- vegan::specaccum(pres, w=w, method="random")
	
		sac_rich <- predict(sac, nSite_rich)
	
		return(list(sac=sac, sac_rich=sac_rich))
	}

	sac_out <- vector("list", length(regs))

	for(r in 1:length(regs)){
	
		t_reg <- regs[r]
	
		if(missing(data_regs)){
			data_in_all0 <- trim_msom(t_reg, gridSize=0.5, depthStratum=100, tolFraction=0.15, grid_stratum=TRUE, plot=FALSE)
	
			nSite_min <- data_in_all0[,lu(stratum), by="year"][,min(V1)]

			data_in_all <- data_in_all0
			setkey(data_in_all, year, stratum, haulid, spp)
		}else{
			data_in_all <- data_regs[[regs[r]]]
		}
		
		u_yrs <- data_in_all[,unique(year)]
		n_spp <- data_in_all[,list(n_spp=lu(spp)), by="year"]
		
		sac_out[[r]] <- vector("list", length(u_yrs))
	
		for(i in 1:length(u_yrs)){
			t_data <- data_in_all[year==u_yrs[i]]
	
			# msg_reg <- toupper(t_data[,unique(reg)])
# 			msg_yr_id <- paste0("Year = ",u_yrs[i])
# 			msg_yr_cnt <-paste0("(", i, " of ", length(u_yrs), ")")
# 			msg_progress <- paste(msg_reg, msg_yr_id, msg_yr_cnt)
# 			cat(paste("\n\n\n", msg_progress, "\n"))
# 			print(Sys.time())
# 			flush.console()
	
			sac_out[[r]][[i]] <- sac_fun(t_data, nSite_rich=nSite_min)
			sac_out[[r]][[i]]$u_yr <- u_yrs[i]
		
			# cat("\n\n")
		
		}
			#
		# append_r <- paste0("_r", r, ".RData")
		# save_name <- gsub("\\.RData", append_r, rm_out[[r]][[i]][[3]]["save_path"])
		# save_name <- gsub("^\\./", "./trawlDiversity/pkgBuild/results/", save_name)
		# save(rm_out, file=save_name, compress="xz")
	
		# cat(paste0("\n\n\n\n\n",  paste(rep("=", 50), collapse=""), "\n\n\n"))
	
	}


	rich0 <- list()
	for(r in 1:length(regs)){
		t_rich <- sapply(sac_out[[r]], function(x)x$sac_rich)
		t_yrs <- sapply(sac_out[[r]], function(x)x$u_yr)
	
		t_dt <- data.table(reg = regs[r], year = t_yrs, richness = t_rich)
		setkey(t_dt, reg, year)
	
		t_dt <- t_dt[CJ(reg=una(reg), year=min(year):max(year))]
	
		rich0[[r]] <- t_dt
	}
	rich <- rbindlist(rich0)
	rich <- rich[!is.na(richness)]
	
	return(rich)
	
}
#
#
#
# # ---- Get Bottom Temperature for all Regions ----
# bt_all0 <- list()
# for(r in 1:10){
#
# 	t_reg <- regs[r]
#
# 	data_in_all0 <- trim_msom(t_reg, gridSize=0.5, depthStratum=100, tolFraction=0.15, grid_stratum=TRUE, plot=FALSE)
#
# 	nSite_min <- data_in_all0[,lu(stratum), by="year"][,min(V1)]
#
# 	data_in_all <- data_in_all0
# 	setkey(data_in_all, year, stratum, haulid, spp)
#
# 	haul_bt <- data_in_all[,list(haul_bt=mean(btemp, na.rm=TRUE)), by=c("reg","year","stratum","haulid")]
# 	strat_bt <- haul_bt[,list(strat_bt=mean(haul_bt, na.rm=TRUE)), by=c("reg","year","stratum")]
# 	reg_bt <- strat_bt[,list(strat_mean=mean(strat_bt, na.rm=TRUE), strat_sd=sd(strat_bt, na.rm=TRUE)), by=c("reg","year")]
#
# 	reg_bt <- reg_bt[CJ(reg=una(reg), year=min(year):max(year))]
#
# 	roll_win <- function(x, stat=c("mean","sd"), width=c("6","9"), ...){
# 		stat <- match.fun(stat)
# 		width <- as.integer(match.arg(width))
#
# 		rollapplyr(x, by=1, width=width, fill=NA, FUN=stat, na.rm=TRUE)
# 	}
#
# 	reg_bt[,strat_mean_6yr_mean:=roll_win(strat_mean, stat="mean", width="6")]
# 	reg_bt[,strat_mean_9yr_mean:=roll_win(strat_mean, stat="mean", width="9")]
#
# 	reg_bt[,strat_mean_6yr_sd:=roll_win(strat_mean, stat="sd", width="6")]
# 	reg_bt[,strat_mean_9yr_sd:=roll_win(strat_mean, stat="sd", width="9")]
#
# 	bt_all0[[r]] <- reg_bt
#
# }
#
# bt_all <- rbindlist(bt_all0)
# bt_all <- bt_all[!is.na(strat_mean)]
#
#





# ---- save results ----


# ---- Time Series of Regional Richness ----
# par(mfrow=c(3,3), mar=c(2,2,0.5,0.5), ps=8, mgp=c(0.5,0.1,0), tcl=-0.1)
# for(r in 1:length(regs)){
# 	if(r == 5){
# 		next
# 	}
#
# 	t_rich <- sapply(sac_out[[r]], function(x)x$sac_rich)
# 	t_yrs <- sapply(sac_out[[r]], function(x)x$u_yr)
#
# 	col <- c("black")
#
# 	if(r == 4){
# 		tri_length <- length(t_rich)
# 		t_rich <- c(t_rich, sapply(sac_out[[r+1]], function(x)x$sac_rich))
# 		t_yrs <- c(t_yrs, sapply(sac_out[[r+1]], function(x)x$u_yr))
# 		col <- rep('red', length(t_rich))
# 		col[1:tri_length] <- "black"
# 	}
#
#
# 	plot(t_yrs, t_rich, type="o", main=regs[r], col=col)
#
# }

