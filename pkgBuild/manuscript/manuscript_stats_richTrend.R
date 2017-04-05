library('trawlDiversity')
library('rbLib')
library('lme4')
library('car')
library('multcomp')

setwd("~/Documents/School&Work/pinskyPost/trawlDiversity")
load("pkgBuild/results/processedMsom/p.RData")


# ---- long-term trends ----
# rich_trend_kendall <- comm_master[,cor.test(reg_rich, year, method="kendall")[c("estimate","p.value")], by='reg']
rich_naive_trend_kendall <- comm_master[,tsTau(year, naive_rich), by='reg'] # comm_master[,Kendall::Kendall(year, reg_rich), by='reg']


# ---- wrapper function to avoid repetition between reg and local rich calcs (and possibly others in future) ----
iter_trend_wrapper <- function(iter_obj, valName){
	cast_val_iter <- as.matrix(dcast(iter_obj, iter~year, value.var=valName)[,iter:=NULL])
	val_year <- as.integer(colnames(cast_val_iter))
	c(list(reg=iter_obj[,unique(reg)]), as.list(post_trend(val_year, cast_val_iter, nSamp=1E4)))
}


# ---- regional richness trend ----
# script originally written to only do this one (regional, not local etc)
tau_list <- list()
pb <- txtProgressBar(min=0, max=length(p))
for(i in 1:length(p)){
	setTxtProgressBar(pb, i)
	tau_list[[i]] <- iter_trend_wrapper(iter_obj=p[[i]]$reg_rich_iter, valName="reg_rich")
}
rich_trend_kendall <- rbindlist(tau_list)


# ---- local richness trend ----
tau_list <- list()
pb <- txtProgressBar(min=0, max=length(p))
for(i in 1:length(p)){
	setTxtProgressBar(pb, i)
	tau_list[[i]] <- iter_trend_wrapper(iter_obj=p[[i]]$local_rich_iter, valName="local_rich")
}
local_rich_trend_kendall <- rbindlist(tau_list)


# ---- save results ----
save(rich_trend_kendall, file="pkgBuild/results/rich_trend_kendall.RData")
save(rich_naive_trend_kendall, file="pkgBuild/results/rich_naive_trend_kendall.RData")
save(local_rich_trend_kendall, file="pkgBuild/results/local_rich_trend_kendall.RData")


