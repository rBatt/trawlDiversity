library(trawlDiversity)
library(trawlData)

ll_mesh <- c(0.25, 0.5, 1.0)
d_mesh <- c(50, 100, 300, 500, 800)
regs <- c("ai", "ebs", "gmex", "goa", "neus", "newf", "sa", "shelf", "wctri")
bin_sites <- replicate(length(regs), matrix(NA, nrow=length(ll_mesh), ncol=length(d_mesh), dimnames=list(lonlat=ll_mesh, depth=d_mesh)), simplify=FALSE)
names(bin_sites) <- regs


ctr <- 0
for(r in 1:length(regs)){
	for(ll in 1:length(ll_mesh)){
		for(d in 1:length(d_mesh)){
			ctr <- ctr + 1
			print(ctr)
			bin_sites[[r]][ll,d] <- tryCatch(trim_msom(reg=regs[r], gridSize=ll_mesh[ll], grid_stratum=TRUE, depthStratum=d_mesh[d], tolFraction=0.15, plot=FALSE, cull_show_up=FALSE)[,lu(stratum)], error=function(cond)NA)
		}
	}
}

setwd("~/Documents/School&Work/pinskyPost/trawlDiversity/pkgBuild/manuscript")
save(ll_mesh, d_mesh, regs, bin_sites, ctr, file="manuscript_binSize.RData")
# save.image("manuscript_binSize.RData")