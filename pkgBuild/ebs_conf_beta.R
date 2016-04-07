
library(trawlDiversity)

setwd("~/Documents/School&Work/pinskyPost/trawl/")
load("trawlDiversity/pkgBuild/results/processedMsom/p_new.RData")

ebs_detect_fish <- fread("Detectability/Ebs Spp ID Confidence Fish.csv", header=TRUE)
ebs_detect_inverts <- fread("Detectability/Ebs Spp ID Confidence Inverts.csv", header=TRUE)


fish_rs <- reshape2::melt(ebs_detect_fish, id.vars="Taxon", value.name="conf", variable.name="year")
inverts_rs <- reshape2::melt(ebs_detect_inverts, id.vars="Taxon", value.name="conf", variable.name="year")

det <- rbind(fish_rs, inverts_rs)
det[,year:=as.integer(as.character(year))]
det[,conf:=as.factor(conf)]
setnames(det, "Taxon", "spp")

spp_ind <- det[,spp] %in% spp.key[,una(spp)]
conf <- det[spp_ind]


p_beta <- p[[1]]$ab[par=="beta",list(beta=mean(value)), by=c("spp","year")]

beta_conf <- merge(p_beta, conf, all=FALSE, by=c("spp","year"))


# detability predicting colonizers?
cdt <- p[[1]]$colonization$col_dt
edt <- p[[1]]$colonization$ext_dt

detect_spp <- conf[year>1983,lu(conf)>1, by="spp"][(V1),una(spp)]
col_spp <- cdt[(col_logic)&year>1983,una(spp)]
ext_spp <- edt[(col_logic)&year>1983,una(spp)]

bad_col <- col_spp[col_spp%in%detect_spp]
bad_ext <- ext_spp[ext_spp%in%detect_spp]

# col_yr <- cdt[(col_logic)&year>1983, list(col=lu(spp)),keyby="year"]
# ext_yr <- edt[(col_logic)&year>1983, list(ext=lu(spp)),keyby="year"]
# r_yr <- merge(col_yr,ext_yr,all=TRUE)
# r_yr[is.na(r_yr)] <- 0
# r_yr[,r:=(col-ext)]
# r_yr[,plot(year, r, type="o")]
#
# col_yr_fix <- cdt[(col_logic)&(!spp%in%bad_col)&year>1983, list(col=lu(spp)),keyby="year"]
# ext_yr_fix <- edt[(col_logic)&(!spp%in%bad_ext)&year>1983, list(ext=lu(spp)),keyby="year"]
# r_yr_fix <- merge(col_yr_fix,ext_yr_fix,all=TRUE)
# r_yr_fix[is.na(r_yr_fix)] <- 0
# r_yr_fix[,r:=(col-ext)]
# r_yr_fix[,plot(year, r, type="o")]

plot(data_all[reg=="ebs", lu(spp), keyby="year"], type="o")
plot(data_all[reg=="ebs"&!spp%in%detect_spp, lu(spp), keyby="year"], type="o")



# Aggregate boxplot
beta_conf[,j={boxplot(beta~conf, xlab="Reported Uncertainty", ylab="MSOM Detectability"); NULL}]

beta_conf[,summary(lm(beta~conf))]

beta_conf[,j={
	summary(lmer(beta~as.factor(conf)+(1|spp)))
}]

lmer(beta~conf|spp, data=beta_conf)


beta_conf[,n_conf:=trawlData::lu(conf),by="spp"]

dev.new()
par(mar=c(2.5,2.5,0.5,0.5), cex=1, ps=10, mgp=c(1.25,0.5,0), tcl=-0.1)
beta_conf[n_conf>1,plot(1,1, ylim=range(beta), xlim=range(as.integer(conf)), xlab="Reported Uncertainty", ylab="MSOM Detectability", type="n")]
beta_conf[n_conf>1,j={
	i_conf <- as.integer(conf)
	print(i_conf)
	tmod <- lm(beta~i_conf)
	print(tmod)
	abline(lm(beta~i_conf))
	# ord <- order(i_conf)
	# lines(i_conf[ord], beta[ord])
}, by=c("spp")]

beta_conf[,beta_range:=diff(range(beta)),by="spp"]


