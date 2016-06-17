

library('trawlDiversity')
library('rbLib')
library('lme4')
library('car')
library('multcomp')

# ========================
# = Richness Time Series =
# ========================
# ---- mean richness ----
# smallest and largest long-term averages
comm_master[,mean(reg_rich), by='reg'][reg%in%c("gmex","shelf")]

# ---- long-term variability ----
comm_master[,sd(reg_rich),by='reg']#[,mean(V1)]
comm_master[,sd(reg_rich),by='reg'][,mean(V1)]

# ---- long-term trends ----
rich_trend_kendall <- comm_master[,cor.test(reg_rich, year, method="kendall")[c("estimate","p.value")], by='reg']
rich_naive_trend_kendall <- comm_master[,cor.test(naive_rich, year, method="kendall")[c("estimate","p.value")], by='reg']

rich_trend_kendall[,BH:=p.adjust(p.value, method="BH")]
setcolorder(rich_trend_kendall, c("reg","estimate","BH","p.value"))
rich_naive_trend_kendall[,BH:=p.adjust(p.value, method='BH')[,"BH"]]
setcolorder(rich_naive_trend_kendall, c("reg","estimate","BH","p.value"))

print(rich_naive_trend_kendall)
print(rich_trend_kendall)


# =======================================
# = Prevalence and Proximity to Absence =
# =======================================
r2.corr.mer <- function(m) {
   lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
   summary(lmfit)$r.squared
}
omega2 <- function(m){
	1-var(residuals(m))/(var(model.response(model.frame(m))))
}

# ---- testing on ebs and column that i think might be relevant ----
# b <- spp_master[reg=='ebs' & !is.na(stretch_type) & propStrata!=0, list(spp, year, propStrata, ext_dist_sign, stretch_type, event_year, ce_categ)]

# too complicated to think about stretch_type while making sure i have less-obvious components right
# (tb <- lmer(propStrata~stretch_type*ext_dist_sign+(ext_dist_sign|spp:stretch_type), data=b))
# summary(tb)
#
# (tb_pc <- lmer(propStrata~ext_dist_sign+(ext_dist_sign|spp), data=b[stretch_type=="post_col"]))
# summary(tb_pc)

big <- spp_master[!is.na(stretch_type) & propStrata!=0]
big<-big[,list(prevalence=propStrata, time=ext_dist, type=as.character(stretch_type), reg=reg, event=as.character(event_year), spp=spp)]

# LMER Model Components:
# Begin Model

# mm_form <- formula(
# 	# Repsonse Variable
# 	prevalence ~
#
# 	# Fixed Effects
# 	time + time:type +
# 	# ext_dist_sign*stretch_type +  # stretch type has no discernable effect on the intercept, which makes sense, because 0 years for either colonization or extinction means 0% strata occupied!
#
# 	 # Random Effects
# 	(1 + time + type:time | reg) + (1 + time + type:time | reg:spp) # + (1 | reg:event))
# )

# mm_form <- formula(
# 	# Repsonse Variable
# 	prevalence ~
#
# 	# Fixed Effects
# 	time*reg +
#
# 	 # Random Effects
# 	(1 + time | spp)
# )

mm_form <- formula(
	# Repsonse Variable
	prevalence ~

	# Fixed Effects
	time + 

	 # Random Effects
	(1 + time | reg) + (1 + time | reg:spp) # + (1 | reg:event))
)

# am <- afex::mixed(mm_form,data=big)
(prev_prox_mod <- lmer(mm_form,data=big))

# LMER Model Metrics
summary(prev_prox_mod) # summary :p
r2.corr.mer(prev_prox_mod) # R^2, v1
omega2(prev_prox_mod) # Omega_naught_squared (R^2 v2)
Anova(prev_prox_mod)
anova(prev_prox_mod)
coef(prev_prox_mod)$reg


# ---- Meaning of intercept in prevalence - proximity model ----
spp_master[ce_categ!='neither' & propStrata!=0,min(propStrata),by=c("spp","reg")][,summary(V1)]
spp_master[ce_categ!='neither' & propStrata!=0,min(propStrata),by=c("spp","reg")][,ecdf(V1)(0.07)]

