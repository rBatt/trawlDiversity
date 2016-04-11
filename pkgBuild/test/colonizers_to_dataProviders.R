
library("rstan")
library("trawlDiversity")
library("rbLib")
library("R2jags")
library("maps")
library("beanplot")
library("fields")

setwd("~/Documents/School&Work/pinskyPost/trawl")


# ==================
# = Get Colonizers =
# ==================
col_info <- list()

regs <- data_all[,una(reg)]

for(r in 1:length(regs)){
	clnz <- get_colonizers(data_all[reg==(regs[r])])
	spp_yr <- clnz$col_ext_dt[col==1 | ext==1, 
		list(year_firstColonization=min(year[col==1]), year_lastExtinction=max(year[ext==1])), 
		keyby="spp"
	]
	cmn <- spp.key[spp%in%spp_yr[,spp], list(common=una(common)), keyby="spp"]
	col_info[[r]] <- data.table(reg=regs[r], merge(spp_yr, cmn, by="spp"))
}

col_info　<- rbindlist(col_info)
col_info[!is.finite(year_firstColonization), year_firstColonization:=NA]
col_info[!is.finite(year_lastExtinction), year_lastExtinction:=NA]

col_only <- col_info[is.finite(year_firstColonization) & is.na(year_lastExtinction)]
col_only[,year_lastExtinction:=NULL]

setkey(col_only, reg, year_firstColonization, spp)

col_only[,lu(spp), by="reg"]


# ========
# = Save =
# ========
save(col_only, file="trawlDiversity/pkgBuild/spp_check/col_only.RData")

col_only[,j={
	fn <- paste0("trawlDiversity/pkgBuild/spp_check/", una(reg), "_colonizers_check.csv")
	write.csv(.SD, file=fn, row.names=FALSE)
},by="reg"]
