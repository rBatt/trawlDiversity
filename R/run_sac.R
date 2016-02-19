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
			d0 <- trim_msom(t_reg, gridSize=0.5, depthStratum=100, tolFraction=0.15, grid_stratum=TRUE, plot=FALSE)
			d <- d0
		}else{
			d <- data_regs[(reg) == regs[r]]
		}
		setkey(d, year, stratum, haulid, spp)
		
		nSite_min <- d[,lu(stratum), by="year"][,min(V1)]
		u_yrs <- d[,unique(year)]
		n_spp <- d[,list(n_spp=lu(spp)), by="year"]
		
		sac_out[[r]] <- vector("list", length(u_yrs))
	
		for(i in 1:length(u_yrs)){
			t_data <- d[year==u_yrs[i]]
	
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
