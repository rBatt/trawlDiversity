
library('trawlDiversity')
library('rbLib')
library("EnvStats")
library('lme4')

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

rich_trend_kendall[,BH:=adjP(p.value)[,"BH"]]
setcolorder(rich_trend_kendall, c("reg","estimate","BH","p.value"))
rich_naive_trend_kendall[,BH:=adjP(p.value)[,"BH"]]
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
b <- spp_master[reg=='ebs' & !is.na(stretch_type) & propStrata!=0, list(spp, year, propStrata, ext_dist_sign, stretch_type, event_year, ce_categ)]

# too complicated to think about stretch_type while making sure i have less-obvious components right
# (tb <- lmer(propStrata~stretch_type*ext_dist_sign+(ext_dist_sign|spp:stretch_type), data=b))
# summary(tb)

(tb_pc <- lmer(propStrata~ext_dist_sign+(ext_dist_sign|spp), data=b[stretch_type=="post_col"]))
summary(tb_pc)

big <- spp_master[!is.na(stretch_type) & propStrata!=0]

# LMER Model Components:
# Begin Model
(prev_prox_mod <- lmer(
	
	# Repsonse Variable
	propStrata ~ 
	
	# Fixed Effects
	ext_dist_sign + ext_dist_sign:stretch_type + 
	# ext_dist_sign*stretch_type +  # stretch type has no discernable effect on the intercept, which makes sense, because 0 years for either colonization or extinction means 0% strata occupied!
	
	 # Random Effects
	(ext_dist_sign | reg/spp) + (1 | reg:event_year)
	
	# End Model
,data=big))

# LMER Model Metrics
summary(prev_prox_mod) # summary :p
r2.corr.mer(prev_prox_mod) # R^2, v1
omega2(prev_prox_mod) # Omega_naught_squared (R^2 v2)



