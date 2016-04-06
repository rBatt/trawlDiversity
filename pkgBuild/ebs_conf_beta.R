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


