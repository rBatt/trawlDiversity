#' Run Species Accumulation Curve Estimate of Richness
#' 
#' Calculate richness using a species accumulation curve
#' 
#' @param regs Character vector of region names
#' @param data_regs Option list of data.tables; list should be named according to \code{regs}
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
	
			sac_out[[r]][[i]] <- sac_fun(t_data, nSite_rich=nSite_min)
			sac_out[[r]][[i]]$u_yr <- u_yrs[i]
			
		}
	
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
