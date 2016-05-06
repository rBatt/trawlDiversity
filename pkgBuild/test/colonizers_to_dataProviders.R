
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

full_ci <- list()

regs <- col_only[,una(reg)]
colonize_refs <- list()
for(r in 1:length(regs)){
	t_reg <- regs[r]
	t_clean <- switch(t_reg,
		ebs = clean.ebs,
		ai = clean.ai,
		goa = clean.goa,
		wctri = clean.wctri,
		wcann = clean.wcann,
		gmex = clean.gmex,
		sa = clean.sa,
		neus = clean.neus, 
		shelf = clean.shelf,
		newf = clean.newf
	)
	
	t_spp <- col_info[reg==t_reg, una(spp)]
	n_yr <- t_clean[,lu(year)]
	bad_flag <- t_clean[,flag%in%c("badJWM","bad","avoid","egg")&!is.na(flag)]
	t_ci <- t_clean[spp%in%t_spp & !bad_flag, list(reg=una(reg), ref=paste(una(ref), collapse=", "), full_firstYear=min(year), full_lastYear=max(year), full_propYear=lu(year)/n_yr), by=c("spp")]
	
	full_ci[[r]] <- t_ci
	
	t_cr <- t_clean[spp%in%t_spp & !bad_flag, list(spp, ref, flag, reg)]
	setkey(t_cr, ref, spp, flag, reg)
	t_cr <- unique(t_cr)
	colonize_refs[[r]] <- t_cr
	
}
colonize_refs <- rbindlist(colonize_refs)
colonize_refs[,reg:=NULL]
setkey(colonize_refs, ref, spp, flag)
colonize_refs <- unique(colonize_refs)

coliz_extinc_full <- rbindlist(full_ci)
setkey(coliz_extinc_full, reg, spp)

# ---- colonization only from full data set ----
setkey(col_only, reg, spp)
col_only_full <- coliz_extinc_full[col_only]
col_only_full[,full_propYear:=round(full_propYear, 2)]
setnames(col_only_full, c("ref","year_firstColonization"), c("orig_tax_id", "subset_firstYear"))
setcolorder(col_only_full, c("reg","spp","orig_tax_id", "common", "subset_firstYear", "full_firstYear","full_lastYear","full_propYear"))


# ---- colonization or extinction from full data set ----
setkey(col_info, reg, spp)
col_ext_full <- coliz_extinc_full[col_info]
col_ext_full[,full_propYear:=round(full_propYear, 2)]
setnames(col_ext_full, c("ref","year_firstColonization","year_lastExtinction"), c("orig_tax_id", "subset_firstYear","subset_lastYear"))
setcolorder(col_ext_full, c("reg","spp","orig_tax_id", "common", "subset_firstYear", "subset_lastYear", "full_firstYear","full_lastYear","full_propYear"))




# ========
# = Save =
# ========
save(colonize_refs, file="trawlDiversity/pkgBuild/spp_check/colonize_refs.RData")

col_ext_full[,j={
	fn <- paste0("trawlDiversity/pkgBuild/spp_check/", una(reg), "_full_col_or_ext_check.csv")
	write.csv(.SD[,list(spp, orig_tax_id, common, full_firstYear, full_lastYear, full_propYear)], file=fn, row.names=FALSE)
},by="reg"]

col_only_full[,j={
	fn <- paste0("trawlDiversity/pkgBuild/spp_check/", una(reg), "_full_colonizers_check.csv")
	write.csv(.SD[,list(spp, orig_tax_id, common, full_firstYear, full_lastYear, full_propYear)], file=fn, row.names=FALSE)
},by="reg"]

save(col_only, file="trawlDiversity/pkgBuild/spp_check/col_only.RData")

col_only[,j={
	fn <- paste0("trawlDiversity/pkgBuild/spp_check/", una(reg), "_colonizers_check.csv")
	write.csv(.SD, file=fn, row.names=FALSE)
},by="reg"]
