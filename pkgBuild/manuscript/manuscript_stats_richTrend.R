library('trawlDiversity')
library('rbLib')
library('lme4')
library('car')
library('multcomp')

setwd("~/Documents/School&Work/pinskyPost/trawlDiversity")

# ---- long-term trends ----
# rich_trend_kendall <- comm_master[,cor.test(reg_rich, year, method="kendall")[c("estimate","p.value")], by='reg']
rich_naive_trend_kendall <- comm_master[,tsTau(year, naive_rich), by='reg'] # comm_master[,Kendall::Kendall(year, reg_rich), by='reg']

load("pkgBuild/results/processedMsom/p.RData")
tau_list <- list()
pb <- txtProgressBar(min=0, max=length(p))
for(i in 1:length(p)){
	setTxtProgressBar(pb, i)
	cast_rich_iter <- as.matrix(dcast(p[[i]]$reg_rich_iter, iter~year, value.var="reg_rich")[,iter:=NULL])
	rich_year <- as.integer(colnames(cast_rich_iter))
	tau_list[[i]] <- c(list(reg=p[[i]]$reg_rich_iter[,unique(reg)]), as.list(post_trend(rich_year, cast_rich_iter, nSamp=1E4)))
}
rich_trend_kendall <- rbindlist(tau_list)

save(rich_trend_kendall, file="pkgBuild/results/rich_trend_kendall.RData")
save(rich_naive_trend_kendall, file="pkgBuild/results/rich_naive_trend_kendall.RData")


