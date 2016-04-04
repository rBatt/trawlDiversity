# rich_contrib <- function(){}
	


library("trawlDiversity")
library("rbLib")
library("maps")

regs <- c("ebs", "ai", "goa", "wctri", "wcann", "gmex", "sa", "neus", "shelf", "newf")
pretty_reg <- c("ebs"="E. Bering Sea", "ai"="Aleutian Islands", "goa"="Gulf of Alaska", "wc"="West Coast US", "gmex"="Gulf of Mexico", "sa"="Southeast US", "neus"="Northeast US", "shelf"="Scotian Shelf", "newf"="Newfoundland")
pr <- names(pretty_reg)

pretty_col <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','navy','#a65628','salmon','#999999')
names(pretty_col) <- pr


# ====================================
# = Get Trimmed Data for Each Region =
# ====================================
data_in_regs <- list()
for(r in 1:length(regs)){
	t_reg <- regs[r]
	t_trimmed <- trim_msom(t_reg, gridSize=0.5, depthStratum=100, tolFraction=0.15, grid_stratum=TRUE, plot=FALSE)
	setkey(t_trimmed, year, stratum, haulid, spp)
	data_in_regs[[regs[r]]] <- t_trimmed
}


data_all <- rbindlist(data_in_regs, fill=TRUE)
# data_all <- data_all[reg!="wcann"]
# data_all[reg=="wctri", reg:="wc"]

# data_all <- data_all[abund>0 & !is.na(abund)]


# =======================================
# = Get Bottom Temperature and Richness =
# =======================================
bt_all <- bt_metrics(pr, data_regs=data_all)
rich <- run_sac(pr, data_regs=data_all)
rich_bt <- merge(bt_all, rich, by=c("reg","year"), all=TRUE)

rich_bt <- rich_bt[reg!="wcann"]
rich_bt[reg=="wctri", reg:="wc"]

rich_bt[,reg_order:=as.numeric(factor(reg, levels=pr, ordered=TRUE))]
setorder(rich_bt, reg_order)
rich_bt[,reg_order:=NULL]




# =================
# = StART PLAYING =
# =================

d <- data_all[reg=="wc"]

d <- d[,j={
	uyr <- .SD[,una(year)]
		
	old <- .SD[year==una(year)[1],una(spp)]
	d_out <- list()
	for(y in 2:lu(year)){
		td <- .SD[(year)==uyr[y]]
		t_spp <- td[,una(spp)]
		colonizers <- c(setdiff(t_spp, old))
		old <- una(c(t_spp))
		
		td[,colonize:=(spp%in%colonizers)]
		d_out[[y]] <- td
	}

	rbindlist(d_out)
	
}]


d <- d[,j={
	uyr <- .SD[,una(year)]
		
	d_out <- list()
	for(y in 1:(lu(year)-1)){
		
		td <- .SD[(year)==uyr[y]]
		t_spp <- td[,una(spp)]
		
		nxt_yr <- .SD[year==una(year)[y+1],una(spp)]
		
		extinct <- c(setdiff(nxt_yr, t_spp))		
		td[,extinct:=(spp%in%extinct)]
		d_out[[y]] <- td
	}

	rbindlist(d_out)
	
}]



d[,j={
	t_tab <- d[,table(stratum, spp, year)]
	n_strat <- nrow(t_tab)
	apply(t_tab>0, c(2,3), function(x)sum(x)/n_strat)
}]