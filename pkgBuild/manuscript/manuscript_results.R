#' ---
#' title: "Marine species' range shifts drive long-term changes in richness"
#' author: "Ryan Batt"
#' date: "2016-08-15"
#' abstract: |
#'       Several regions show long-term changes in richness, generally increases. The rest of this study focuses on understanding these long-term changes through the lens of the geographic ranges of individual species and how geographic range relates to the colonization and extinction of those species, and thus to changes in richness.   
#' 
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 4
#'     fig_caption: true
#'     theme: "readable"
#'     template: default
#'   pdf_document:
#'     toc: true
#'     toc_depth: 4
#'     template: latex-ryan.template
#'     fig_caption: true
#' geometry: margin=1.0in
#' lineno: true
#' lineSpacing: false
#' titlesec: true
#' documentclass: article
#' placeins: true
#' ---

#+ deleted-pandoc-headers, include=FALSE, echo=FALSE
# #'      pandoc_args: [
# #'      "--chapters"
# #'      ]


#+ setup, include=FALSE, echo=FALSE
# =================
# = Load Packages =
# =================
# #' date: "2016-08-10"
# #' title: "manuscript_results.R"
library('trawlDiversity')
library('rbLib')
library('lme4')
library('car')
library('multcomp')
library("data.table")
library('maps')
library('raster')
library('spatstat')
library("spdep")
library("piecewiseSEM")
library("viridis")

# Report
library(knitr)
library(rmarkdown)
library(xtable)
library(kfigr) # devtools::install_github("github mkoohafkan/kfigr") # ?
library(stargazer)
library(texreg)

# Other
library(rbLib) # library(devtools); install_github("rBatt/rbLib")


# ================
# = Report Setup =
# ================
doc_type <- c("html", "pdf")[1]
table_type <- c("html"="html", "pdf"="latex")[doc_type]
options("digits"=5) # rounding output to 4 in kable() (non-regression tables)
o_f <- paste(doc_type, "document", sep="_")

# problem with pdflatex in El Capitan? It might be directories. Check http://pages.uoregon.edu/koch/FixLink.pkg

# render!
# rmarkdown::render(
# 	"~/Documents/School&Work/pinskyPost/trawlDiversity/pkgBuild/manuscript/manuscript_results.R",
# 	output_format=o_f,
# 	output_dir='~/Documents/School&Work/pinskyPost/trawlDiversity/pkgBuild/manuscript',
# 	clean = TRUE
# )

Sys.setenv(PATH=paste0("/Library/TeX/texbin:",Sys.getenv("PATH")))

opts_chunk$set(
	fig.path = 'manuscript_report/', 
	cache.path='manuscript_report/',
	echo=TRUE, 
	include=TRUE, 
	cache=F,
	results='asis',
	warning=FALSE,
	fig.show="hold",
	fig.lp = if(o_f=="html_document"){"**Figure.**"}else{NULL}
)


# setwd("~/Documents/School&Work/pinskyPost/trawlDiversity/pkgBuild/manuscript")
source("../manuscript/manuscript_figures_functions.R")
source("../manuscript/manuscript_stats_functions.R")
source("../manuscript/manuscript_data_functions.R")
# source("../manuscript/fig_tbl_number.R")
# eval(fig_tbl_number())

#' \FloatBarrier  
#' ***  
#+ analysisOptions, echo=TRUE, results="markup"
rsType <- c("observed", "rarefied", "msom")[2]
rSMet_base <- c("observed"="propStrata", "rarefied"="range_size_samp", "msom"="range_size_mu")[rsType]
rSMet_annComm <- c("observed"="propStrata_avg", "rarefied"="range_size_samp_avg", "msom"="range_size_mu_avg")[rsType]
rSMet_ltComm <- c("observed"="propStrata_avg_ltAvg", "rarefied"="range_size_samp_avg_ltAvg", "msom"="range_size_mu_avg_ltAvg")[rsType]

rR <-  c("observed"="naive_rich", "rarefied"=NA, "msom"="reg_rich")[3]
lR <-  c("observed"="local_rich_obs", "rarefied"="local_rich_samp", "msom"="local_rich")[2]
betaD <-  c("observed"="beta_div_obs", "msom"="beta_div_mu")[1]


#' There are four primary metrics:  
#'  - range size  
#'  - regional richness (gamma)  
#'  - local richness (alpha)  
#'  - spatial beta diversity  
#'   
#' Each of these metrics can be based on raw observations or on MSOM estimates, and some can also be based on a rarefied (resampled, or 'samp') set of observations.  
#'   
#' Furthermore, the range size metric can be based on a few combinations of cross-year and/or cross-species averages. The base metric for range size (no averaging) is used to to examine how range size changes with time until extinction and time after colonization. A cross-species average of a given year's observed range sizes (_avg) is used for the rare species to see how this sub-group changes in its geographic distribution, or when making specialized comparisons. A cross-year average (ltAvg) is used to relate a species' 'typical' range size to the total number of colonizations and extinctions it experienced during the time series. A cross-year (all years a species was observed) followed by a cross-species (species present in a given year) average (_avg_ltAvg) is used to compare species richness to the typical range size of species in the community.


# ============
# = Richness =
# ============
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Basic Info
#' ##Number of Speices
#+ basicInfo, echo=TRUE
#' Number of species across all regions: `r spp_master[, length(unique(spp))]`  
#' 
#' 



#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Species Richness
#' ###Table. Summary of richness mean and variability
#+ richnessSummary-basic-table, echo=FALSE
kable(
	data.table(
		comm_master[,list(mu=mean(reg_rich), sd=stats::sd(reg_rich)),by='reg'], 
		comm_master[, list(mu_sd=mean(.SD[,stats::sd(reg_rich),by='reg'][,V1]))]
	), 
	caption="Long-term mean and standard deviation in MSOM richness, and average of those sd's. Gmex and Shelf have the highest and lowest long-term averages in species richness (MSOM)."
)

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Figure 1. MSOM richness time series
#+ Richness-ts-fig, fig.height=5, fig.width=3.5, ecal=TRUE, echo=TRUE, fig.cap="**Figure 1.** Time series of MSOM estimates of region richness. Each point is the posterior mean of regional richness in a year. Lines indicate long-term trends from fitted values of linear regression models predicting richness from time."
# load regression results
load("../../pkgBuild/results/rich_naive_trend_kendall.RData")
rich_naive_trend_kendall[reg!="wcann",BH:=p.adjust(taup, method='BH')]
rich_naive_trend_kendall <- rich_naive_trend_kendall[
	reg!="wcann",list(reg=reg, estimate=tau, BH=BH, pvalue=taup)
]
load("../../pkgBuild/results/rich_trend_kendall.RData")
rich_trend_kendall[reg!="wcann",BH:=p.adjust(pvalue, method="BH")]
rich_trend_kendall <- rich_trend_kendall[
	reg!="wcann",list(reg=reg, estimate=tau, BH=BH, pvalue=pvalue)
]

# make plot
load("../results/rich_trend_kendall.RData") # load rich trend stats
richness_ts()
pdf("../../submitted_pdfs/final/Fig1.pdf", width=3.228, height=4.6114)
richness_ts()
dev.off()

#' Looking for trends in species richness in both the naive estimates and the MSOM estimates. Using Kendall's Tau_b, calculated for 1E4 resamplings of the posterior of richness. Tau is also calculated using a method that removes serial correlation to achieve independence of observations so that the p-values are correct.  
#'   
#' In the Naive estimates, `r rich_naive_trend_kendall[reg!='wcann', sum(pvalue<=0.05)]` regions had significant $\tau_b$.  
#'   
#' ##Table 1. MSOM tau
#+ msomTau-table, echo=FALSE
kable(rich_trend_kendall, 
	caption="**Table 2.** MSOM richness trends in each region. Estimate is Kendall's Tau_b, BH is the Benjamini-Hochberg corrected p-value, and p.value is the original p-value. Trends are calculated for 1E4 resampled combinations of the richness posterior. For each resampling, independence of observations is achieved before estimating Tau by removing serial correlation from resampled time series. This procedure retains the integrity of p-values."
)
#' In the MSOM estimates, `r rich_trend_kendall[reg!='wcann', sum(pvalue<=0.05)]` regions had significant $\tau_{b}$.  
#'   
#' Overall, ~half the regions show positive trends. No regions show significant negative trends (although SEUS is close, depending on the analysis). It's a bit surprising how much removing the autocorrelation seemed to impact some of the trends (AI in particular, I think).  
#'   
#' **Manuscript paragraph:**  
#' Estimated slopes for long-term trends (Kendall's $\tau_b$) in richness were positive for most regions. For Naive estimates, all the seven positive $\tau$ were also significantly different from 0, whereas the two negative $\tau_b$, Aleutian Islands and Southeast US, were not (**Table S2**). Long-term trends in MSOM estimates of richness had the same sign as Naive trends, except for Aleutian Islands, but now the only significant $\tau_b$ were the following four regions: Eastern Bering Sea (**$\tau_b$ = `r rich_trend_kendall[reg=='ebs',tau]`**), West Coast US (**$\tau_b$ = `r rich_trend_kendall[reg=='wctri',tau]`**), Scotian Shelf (**$\tau_b$ = `r rich_trend_kendall[reg=='shelf',tau]` **), and Newfoundland (**$\tau_b$ = `r rich_trend_kendall[reg=='newf',tau]` **) (**Table 1**). Richness trends were not significant in the Gulf of Mexico, Gulf of Alaska, and Northeast US, Aleutian Islands, Southeast US.   

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Colonization and Extinction
#' ##Colonization and Extinction Summary
#+ col-ext-smry-categs, echo=TRUE
#' The number of species-regions: `r spp_master[,length(unique(paste(reg,spp)))]`  
#' The number of *core* species: `r spp_master[ce_categ=='neither',length(unique(paste(reg,spp)))]`  
#' The number of *both* species: `r spp_master[ce_categ=='both',length(unique(paste(reg,spp)))]`  
#' The number of *colonizing* species: `r spp_master[ce_categ=='colonizer',length(unique(paste(reg,spp)))]`  
#' The number of *leaving* species: `r spp_master[ce_categ=='leaver',length(unique(paste(reg,spp)))]`  


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' #Trends in Range Size Near Colonizations and Extinctions
#' ##Figure 2. Changes in Geographic Range before Extinction and after Colonization
#+ rangeSize_ColExt, fig.width=6, fig.height=3, fig.cap="**Figure 4.** Geographic range size vs years until extinction (A) and years after colonization (B). For visualization purposes, range size is averaged across species for each unique value on each axis, and a linear model fit through this average. Statistics in main text use unaggregated data. The horizontal axes were formulated as time until (since) the nearest upcoming (previous) absence. Because range size must be zero when either horizontal axis has a value of zero, points at (0,0) were excluded from figures and analyses."
rangeSize_absenceTime(rSMet_base)
pdf("../../submitted_pdfs/final/Fig2.pdf", width=6.81102, height=3*(6.81102/6))
rangeSize_absenceTime(rSMet_base)
dev.off()

print(rSMet_base)
#' Range size declines near an absence much more consistently than does range density; both are (relatively) low just before extinction and just after colonization. However, range density has much more variable intercepts among regions, whereas range size does not.   
#'   
#' This makes sense, at least somewhat, because colonization and extinction events are defined at the site level; though the outcome isn't necessitated by this formulation, because size could drop suddenly. In fact, when a species is absent, both its range size and its range density must be 0 (though, range density is technically calculated for only those sites that are occupied, so I supposed it's technically undefined according to the equations I'm using).  
#'   
#' I think the regressions for range size should omit an intercept, while the regressions for range density should have it. This might be hard to justify fully *a priori* (though see my thinking in previous paragraph), so I'll probably just do a model selection and maybe discuess the difference if one model has an intercept and the other does not.  
#'   
#'   
#' ###Range Size Time to Absence Regressions
#+ rangeSizeDensity-ColExtTime-reg-simple
rangeTimeDT <- make_rangeTime(sizeName=rSMet_base)

# Set up -- region names and empty named lists
ur <- rangeTimeDT[,unique(reg)]
sTime_reg_mods0 <- structure(vector("list",length(ur)), .Names=ur)
sTime_reg_mods1 <- structure(vector("list",length(ur)), .Names=ur)
sTime_reg_mods2 <- structure(vector("list",length(ur)), .Names=ur)
sTime_reg_mods3 <- structure(vector("list",length(ur)), .Names=ur)

# Fit 3 models for each region
for(r in 1:length(ur)){
	t_reg <- ur[r]
	t_dat <- rangeTimeDT[reg==t_reg]
	sTime_reg_mods0[[t_reg]] <- lme4::lmer(size ~ time + (time-1|spp), data=t_dat)
	sTime_reg_mods1[[t_reg]] <- lme4::lmer(size ~ time + (time|spp), data=t_dat)
	sTime_reg_mods2[[t_reg]] <- lme4::lmer(size ~ time + type + (time|spp), data=t_dat)
	sTime_reg_mods3[[t_reg]] <- lme4::lmer(size ~ time * type + (time|spp), data=t_dat)
}

# Summarize each model
sTime_reg_smry0 <- smry_modList2(sTime_reg_mods0)
sTime_reg_smry1 <- smry_modList2(sTime_reg_mods1)
sTime_reg_smry2 <- smry_modList2(sTime_reg_mods2)
sTime_reg_smry3 <- smry_modList2(sTime_reg_mods3)

# Captions for each model
# Also a caption in case I want to report all models for all regions
capS0 <- sTime_reg_smry0[,unique(mod_call)]
capS1 <- sTime_reg_smry1[,unique(mod_call)]
capS2 <- sTime_reg_smry2[,unique(mod_call)]
capS3 <- sTime_reg_smry3[,unique(mod_call)]
timeRegMods_cap <- "Summary statistics for fits of predicting range size from years before extinction and years after colonization. These are mixed effect models of the form "

# Compare AIC values from each model fit
compAIC <- rbind(
	cbind(sTime_reg_smry0[,list(reg,MargR2,CondR2,AIC)], mod=capS0),
	cbind(sTime_reg_smry1[,list(reg,MargR2,CondR2,AIC)], mod=capS1),
	cbind(sTime_reg_smry2[,list(reg,MargR2,CondR2,AIC)], mod=capS2),
	cbind(sTime_reg_smry3[,list(reg,MargR2,CondR2,AIC)], mod=capS3)
)
setkey(compAIC, reg) # sort

# Get best model from each region, and from each region list best overall model
bestEach <- compAIC[,list("mod_call"=unique(mod)[which.min(AIC)]),by=c('reg')] # combinations of each region and its best model
bestOverall_name <- bestEach[,j={bm <- table(mod_call); names(bm)[which.max(bm)]}] # name of best overall model
bestOverall <- bestEach[,CJ(reg=reg,mod_call=bestOverall_name)] # get combos for each region w/ best overall model
bestEachOverall <- unique(rbind(bestEach,bestOverall)) # combinations for best overall and each
setkey(bestEachOverall, reg, mod_call) # sort

# get the worthy models
allSmry <- rbind(sTime_reg_smry0, sTime_reg_smry1, sTime_reg_smry2, sTime_reg_smry3, fill=TRUE) # combine results from all models and regions
bestModels <- allSmry[bestEachOverall, on=c('reg','mod_call')] # for each region, only select models that were best in that region or best overall

# remove columns that are NA for all rows b/c models w/ those parameters were never winners
loserNames <- sapply(bestModels, function(x)all(is.na(x))) # names of columns only pertaining to non-best models
loserNames <- names(loserNames)[loserNames] # okay, get the names for real, not just logic
bestModels[,c(loserNames):=list(NULL)] # drop the loser columns

#' ###Range Size Time to Absence Table
#+ rangeSizeDensity-ColExtTime-reg-compareModels, echo=FALSE
# kable(bestEach, caption="Which models for predicting range size from time to event were best in each region?")
# kable(compAIC, caption="Comparing fit of models of varying complexity; more complex models test for differences in the slopes or intercepts of pre-extinction and post-colonization trends.")
kable(bestModels, caption="Shows each region's best model, and shows the model that was most often the best for all regions (including those for which it wasn't the best). If all regions had the same best model, only 1 model is shown.")

#' ###Text summarizing pre-extinction contraction and post-colonization expansion slopes and fits
#+ rangeSizeDensity-ColExtTime-reg-compareModels-text
chosenModels <- bestModels[!is.na(MargR2) & mod_call!=("size~time+type+(time|spp)")]
mR2Range <-  chosenModels[,range(MargR2)]
cR2Range <- chosenModels[,range(CondR2)]
pRange <- chosenModels[,range(pval_time)]
slopeRange <- chosenModels[,range(time)]
slopeAvg <- chosenModels[,mean(time)]

sd_ranef <- sapply(sTime_reg_mods1, function(x)sd(coef(x)[[1]][,'time']))
mean_ranef <- sapply(sTime_reg_mods1, function(x)mean(coef(x)[[1]][,'time']))
dRange_ranef <- sapply(sTime_reg_mods1, function(x)diff(range(coef(x)[[1]][,'time'])))
ranef_smry <- data.table(reg=names(sd_ranef), sd_ranef=sd_ranef, mean_ranef=mean_ranef, dRange_ranef=dRange_ranef)
ranef_smry[,c("sd_mean","dR_mean"):=list(sd_ranef/mean_ranef, dRange_ranef/mean_ranef)]
#' ####Summary of Range Size ~ Time Until Extinction Regression Fits
#' after correcting for multiple tests, all **$p \leq `r pRange[2]`; `r mR2Range[1]` \leq mR^2 \leq `r mR2Range[2]`; `r cR2Range[1]` \leq cR^2 \leq `r cR2Range[2]`$**. For slope, **$`r slopeRange[1]*10*100` \leq \beta \leq `r slopeRange[2]*10*100`, \bar{\beta} = `r slopeAvg*10*100`$** percent per decade.   
#' "Indeed, variation among species' slopes was similar to the average slope(**$\hat{\sigma}_{\beta} = `r mean(sd_ranef*1000)`$** percent per decade)"



#' ###Conclusion for Range Change before Extinction/ after Colonization
#' As stated with the larger (pooled regional data) regressions, range size responds more consistently/ strongly to an approaching extinction/ departure from colonization than does range density. Range size shrinks as extinction approaches, increases as time since colonization increases. In both cases, there was only 1 region for which the direction (into extinction, out of colonization) mattered: for range density, it was the Southeast US (effect of pre-extinction direction = 0.037, p=0.043 [BH correction], Table 12), and for range size it was Scotial Shelf (effect = -0.021, p = 0.002). Time was a significant predictor of range size in all regions; time was *not* a significant predictor of range density in Newfoundland (p = 0.846 [BH]), and the Scotial Shelf (p = 0.916).  
#'  
#' These results support the hypothesis that species rarity is closely associated with proximity to extinction/ colonization, and therefore the probability that the species will contribute to a change in species richness. Furthermore, these relationships have not been directly quantified for entire assemblages. Competing hypotheses exist regarding how range should change approaching extinction and after colonization. We find that the change in range is similar regardless of direction. However, spatial scale was important. Although slopes were similar, the intercept for density was far larger than for range size --- the fraction of tows containing a species in occupied sites may remain high even when extinction is imminent, and may similar be large even if colonization was recent. Range size, on the other hand, was much closer to 0 near an absence.  
#'   
#' Overall, the results suggest that the spatial footprint of individual species is important for understanding changes in species richness. Furthermore, because species contributing most to the dynamics of richness are those that repeatedly colonize and go extinct, it is meaningful to look at a species' long-term rarity in order to gauge whether it is likely to contribute to those long-term richness changes. Determining what drives the geographic range of individual species is probably a powerful way to anticipate richness changes.  
#'   

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Richness and Range Size
#' ##Do Transient Species Have Smaller Ranges than Core Species?
#+ comp-transCore-range
transRange <- spp_master[,list(reg=reg, spp=spp, transient=(ce_categ!='neither'),rangeSize=range_size_samp)]
ur <- names(pretty_reg)
nr <- length(ur)
tR_mod <- structure(vector('list', length(ur)), .Names=ur)
for(r in 1:nr){
	tR_mod[[ur[r]]] <- lmer(rangeSize~transient + (1|spp), data=transRange[reg==ur[r]])
}
tR_mod_smry <- smry_modList2(tR_mod)
kable(tR_mod_smry)
transDiffMean <- tR_mod_smry[,mean(transientTRUE)]*100
sign_tDM <- sign(transDiffMean)
#' In mixed effects models with intecepts varying among species and transient *versus* core as a categorical predictor, transient species had range sizes that were **`r abs(transDiffMean)`** (% occupancy) **`r c("smaller","larger")[(sign_tDM>0)+1]`** than the ranges of core species (average intercept adjustment; all **$p \leq$ `r tR_mod_smry[,max(BH_transient)]`** after correcting for multiple tests).  
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Richness and Geographic Range
#' ##Figure 3. Richness versus Community Range Index
#+ rich-geo-rangeSize, fig.width=3.5, fig.height=3.5, fig.cap="**Figure 3.** Species richness vs geographic range size. Range size is presented as each species' long-term average of the proportion of sites it occupied. Solid lines are linear regressions with MSOM richness as the response and the horizontal axis and an intercept as the predictors."
rich_geoRange(rSMet_ltComm, leg=TRUE, legPan=1, panLab=FALSE)
pdf("../../submitted_pdfs/final/Fig3.pdf", width=3.228, height=3.228)
rich_geoRange(rSMet_ltComm, leg=TRUE, legPan=1, panLab=FALSE)
dev.off()

#' Range size is a pretty good predictor of species richness. I think I had originally missed the range size relationship b/c I hadn't done the same aggregating procedure. The interpretation I have is that richness is highest when you have a bunch of rare species.  
#'   
#'   
#' ###Regressions relating Richness to Community Range Index
#+ rich-rangeSize
range_reg <- make_range_reg(dens="propTow_occ_avg", size=rSMet_ltComm)

# Fit different models to the whole data set
# Models vary in which/ how parameters vary among regions
rSize_mods <- list()
rSize_mods[[1]] <- lm(rich ~ size, data=range_reg) # simple
rSize_mods[[2]] <- lm(rich ~ size*reg, data=range_reg) # slope and intercept ind. among regs
rSize_mods[[3]] <- lme4::lmer(rich ~ size + (1|reg), data=range_reg) # intercept varies randomly among regs
rSize_mods[[4]] <- lme4::lmer(rich ~ size + (size|reg), data=range_reg) # slope & int. vary randomly among regs
rich_size_smry <- smry_modList(rSize_mods) # summarize to get R^2, AIC, p-vals, etc

# Fit simple model to each region separately 
ur <- range_reg[,unique(reg)]
rSize_reg_mods <- structure(vector("list",length(ur)), .Names=ur)
for(r in 1:length(ur)){
	rSize_reg_mods[[r]] <- lm(rich ~ size, data=range_reg[reg==ur[r]])
}
rich_s_reg_smry <- data.table(reg=ur, smry_modList(rSize_reg_mods)) # summarize
setnames(rich_s_reg_smry, old=c("Marginal","Conditional","p.value"), new=c("MargR2","CondR2","pval")) # rename
setkey(rich_s_reg_smry, reg, Class, mod_call, predictor) # sort
rich_s_reg_smry[, BH:=round(p.adjust(pval, "BH"), 3), by=c("mod_call", "predictor")] # correct p-vals for multiple tests
smry_formula <- formula(reg+Class+mod_call+MargR2+CondR2+AIC~predictor) # formula to organize
rich_s_reg_smry <- dcast(rich_s_reg_smry, smry_formula, value.var=c("pval", "BH")) # organize; helpful if multiple predictors

#+ rich-rangeSize-tables, echo=FALSE
# kable(
# 	rbind(rich_size_smry),
# 	caption="Summary of regression models involving richness and community metric of range size."
# ) # different kinds of models
kable(
	rich_s_reg_smry, 
	caption="Regression models of richness predicted by community metric of range size, but each region has a separate model."
) # same model applied to each reg sep
kable(
	rich_s_reg_smry[,lapply(.SD, base::mean), by="mod_call"], 
	caption="Average of above region-specific rich ~ size models."
) 

rich_cri_maxP <- rich_s_reg_smry[,max(BH_size)]
rich_cri_rangeR2 <- rich_s_reg_smry[,range(MargR2)]
rich_cri_meanR2 <- rich_s_reg_smry[,mean(MargR2)]
#' **Manuscript Text:**  
#'  Furthermore, species richness was negatively correlated with the community range index (CRI) in each region (Fig. 3; separate linear regression for each region, for slope all corrected **$p \leq$ `r rich_cri_maxP`**, **$`r rich_cri_rangeR2[1]` \leq R^2 \leq `r rich_cri_rangeR2[2]`$**, average **$R^2 = `r rich_cri_meanR2`$ **).  
#'   
#' **Other Notes:**  
#' All models are pretty good predictors. Well, the most basic model kinda sucks I guess. It needs to account for some of the between-region variation.  
#'   
#' *Predicting Richness: Range Size or Density?*
#' As far as picking one or the other, it doesn't end up mattering much. Range size is a lot better than density in NEUS, and density outperforms size in AI. Otherwise, size as a slight edge over density on average, although both predictors are significant in all regions.  
#'   
#' *Conclusion for Richness and Geographic Range*
#' The two metrics of geographic range are well correlated. Furthermore, richness can be predicted pretty well using regressions with either as a predictor. There are large differences among regions, though. This is probably because richness is not readily comparable among most regions. Regions vary mostly in their intercept values, and they have fairly similar slopes (though they are not identical, and model fits improve when allowing slopes to vary among regions; it's just that the improvement is small compared to allowing intercepts to vary among regions).  
#'  
#' The interpretation of the result that geographic distribution predicts species richness is likely associated with species rarity. When the average range density or range size of a community is low, it means it has a lot of species that are rare (at either spatial scale). It's these rare species that come and go, and form the dynamics of richness that we observe. When that dynamical value is high, it implies that an above-average number of the dynamic species are present. Because those species are transient (dynamic), they are also probably rare.  
#'   

#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' #Range Size Over Time For Core & Transients
#' ##Figure 4. Range size over time for transient & core
#+ rareExpansion, fig.width=6.5, fig.height=6.5, results='asis'
eval(neitherPres_subColor())
plot_rangeSize_FullTrans(range_type="range_size_samp")
pdf("../../submitted_pdfs/final/Fig4.pdf", width=6.81102, height=6.81102)
plot_rangeSize_FullTrans(range_type="range_size_samp")
dev.off()

#' ###Range Size over time for transient & core -- Stats
#+ rareExpansion-stats
transMod <- structure(vector("list",length(ur)), .Names=ur)
coreMod <- structure(vector("list",length(ur)), .Names=ur)
ctMod <- structure(vector("list",length(ur)), .Names=ur)
ctMod2 <- ctMod
for(r in 1:length(ur)){
	rr <- ur[r]
	td <- spp_master[eval(rE_logic1)]
	td[,year:=year-min(year)]
	transMod[[r]] <- lmer(range_size_samp~year+(1|spp), data=td)
	
	td2 <- spp_master[eval(rE_logic3)]
	td2[,year:=year-min(year)]
	coreMod[[r]] <- lmer(range_size_samp~year+(1|spp), data=td2)
	
	td3 <- spp_master[reg==rr & present==1]
	td3[,year:=year-min(year)]
	td3[,categ2:=c("core","transient")[(ce_categ!="neither")+1L]]
	ctMod[[r]] <- lmer(range_size_samp~year*categ2+(1|spp), data=td3)
	
	ctMod2[[r]] <- lm(range_size_samp~year*categ2, data=td3)
	# ctMod2[[r]] <- lmer(range_size_samp~year+(year|spp), data=td)
	# ctMod2[[r]] <- lmer(range_size_samp~year+(year|spp), data=td3)
	# ctMod2[[r]] <- lmer(range_size_samp~year*categ2+(year|spp), data=td3)
	# td3[,auto.arima(range_size_samp, xreg=as.matrix(cbind(year=year, categ2=(categ2=="transient"), year_categ2=(year*(categ2=="transient")))), allowdrift=FALSE, d=0)]

}
transExpansionStats <- smry_modList2(transMod)
coreExpansionStats <- smry_modList2(coreMod)
ctExpansionStats <- smry_modList2(ctMod)
# ctExpansionStats2 <- smry_modList2(ctMod2)

meanNum <- function(x, ...){
	if(storage.mode(x) %in% c("double", "integer", "logical")){
		mean(x, ...)
	}else{
		NA
	}
}

kable(rbind(transExpansionStats,lapply(transExpansionStats, meanNum), transExpansionStats[BH_year<0.05, lapply(.SD, meanNum)][,reg:="signifYear"]))
kable(rbind(coreExpansionStats,lapply(coreExpansionStats, meanNum), coreExpansionStats[BH_year<0.05, lapply(.SD, meanNum)][,reg:="signifYear"]))
kable(rbind(ctExpansionStats,lapply(ctExpansionStats, meanNum), ctExpansionStats[BH_year<0.05, lapply(.SD, meanNum)][,reg:="signifYear"], ctExpansionStats[`BH_year:categ2`<0.05, lapply(.SD, meanNum)][,reg:="signifYear:categ2"]))

#' ###Range Size over time for transient vs core -- Text 
#+ rareExpansion-description, fig.width=6.5, fig.height=6.5, results='markup'
#' Examining how range size changed over time for "core" and "transient" groups of species. Core species are those that are present in all years. Transient species are those that are not present in all years. **`r transExpansionStats[,sum(BH_year<0.05)]`** regions had significant trends in range size for transient species (average slope for regions with significant trend = **`r transExpansionStats[,mean(year[BH_year<0.05])]`**), and **`r coreExpansionStats[,sum(BH_year<0.05)]`** regions had significant trends in range size for core species (**`r coreExpansionStats[,paste(reg[(BH_year>=0.05)], collapse = ", ")]`** did not have significant trends for core species; averge slope for regions with significant trend = **`r coreExpansionStats[,mean(year[BH_year<0.05])]`**).  
#'   
#' Given that many regions showed positive trends in range size for both species groups, it is prudent to compare trends between these groups, to see if they differ. To this end, I analyzed the groups together, and tested the significance of an interaction term that adjusted the a baseline slope if the species was in the transient group.  
#'   
#' Note that transient species are very much expected to have a smaller range size than the core group, particularly at the beginning of the time series; I'm not dynamically reporting the significance of the difference in intercepts here, but at the time of writing all regions showed significant differences in intercept between groups, with the transient group having a smaller intercept than the core group.
#'   
#' Interestingly, **`r ctExpansionStats[,sum(BH_year<0.05)]`** regions had significant trends (main effect), and **`r ctExpansionStats[,paste(reg[(BH_year<0.05 & year <0)], collapse = ", ")]`** had significant negative trends, while **`r ctExpansionStats[,paste(reg[(BH_year<0.05 & year >0)], collapse = ", ")]`** had significant positive trends. **`r ctExpansionStats[,paste(reg[(get("BH_year:categ2")>=0.05)], collapse = ", ")]`** were the regions without significant interactions -- (at the time of this writing) these are not the same regions that lacked significant trends in the core group. If we sum the main effect and transient-interaction together (ignoring significance), we see that **`r ctExpansionStats[,paste(reg[(year+get("year:categ2transient"))<0], collapse = ", ")]`** was/were the only region/s with negative estimated trends for transient species; all others had positive trends in range for transients. Some regions had significant trends for both core and transients, but the direction of these trends differed between the two groups: **`r ctExpansionStats[,paste(reg[(BH_year<0.05 & year <0 & get("BH_year:categ2")< 0.05 & get("year:categ2transient")>0)], collapse = ", ")]`** had core species with trends in range size that were significantly less than 0 but transient species with trends in range size that were significantly greater than 0.  
#'   
#' Southeast U.S. core species long-term slope **$\beta_Y=$`r ctExpansionStats[reg=='sa', year]*1E3`** (percent per decade)  
#' Newfoundland core species slope **$\beta_Y=$`r ctExpansionStats[reg=='newf', year]*1E3`** (percent per decade)  
#' Scotian Shelf core species slope **$\beta_Y=$`r ctExpansionStats[reg=='shelf', year]*1E3`** (percent per decade)  
#'   
#' For core species with positive long-term slopes, the average was **$\bar{\beta}_{Y+}=$`r ctExpansionStats[year>0, mean(year)]*1E3`** (percent per decade). After correcting for multipel testing, all **$p \leq$ `r ctExpansionStats[,max(BH_year)]`**
#'   
#' However, the interaction term indicated that the slopes of core and transient species were different in **`r ctExpansionStats[,sum(get("BH_year:categ2")<0.05)]` of the regions** (**$p \leq$ `r ctExpansionStats[get("BH_year:categ2")<0.05,max(get("BH_year:categ2"))]`**), including positive interactions in Newfoundland (**$\beta_{Y \times T} =$ `r ctExpansionStats[reg=='newf', get("year:categ2transient")]*1E3`**) and Scotian Shelf (**$\beta_{Y \times T} =$ `r ctExpansionStats[reg=='shelf', get("year:categ2transient")]*1E3`**). Excluding the Southeast U.S., the average long-term trend in range size for transient species was **`r ctExpansionStats[reg!='sa', mean(get("year:categ2transient")+year)]*1E3`**  



#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Supplement
#' ##Figure S1. Map of Site Locations
#+ siteMap, fig.width=7.42, fig.height=3
plotSiteMap()


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##List of Species sampled
#+ spp_sampled_list
spp_master[,j={cat('\n\n\n\n\n');cat('**',reg[1],'**',sep=''); cat('\n\n'); cat(paste(sort(unique(spp)),collapse=', '));cat('\n');NULL},by='reg']

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Figure S2. MSOM - naive scatter
#+ Richness-msom-naive-scatter-fig, fig.height=3.5, fig.width=3.5, fig.cap="**Figure S1.** MSOM richness vs naive richness"
naive_msom_scatter()
#' MSOM richness and Naive richness are pretty similar. MSOM richness is probably more accurate, or is at least more conservative b/c it has fewer significant trends. But their similarity should help justify using observed presences/ absences in other analyses.  
#'   
#' **Manuscript paragraph:**  
#' MSOM estimates of richness were greater than Naive estimates, but the two methods produced similar temporal dynamics in richness (Figure S1). Henceforth, we report MSOM richness estimates. The greatest long-term average richness was the Gulf of Mexico (142.3), lowest in the Scotian Shelf (46.3). The inter-region average of long-term standard deviations of richness was 5.3; the Aleutian Islands showed the lowest variability (sd = 2.4), and the Gulf of Alaska was the most variable (sd = 8.5).  
#'   
#+ Richness-msom-naive-scatter-stats, results='markup'
ur <- comm_master[,unique(reg)]
nr <- length(ur)
naive_msom_mod <- list()
for(r in 1:nr){
	td <- comm_master[reg==ur[r]]
	naive_msom_mod[[ur[r]]] <- lm(reg_rich~naive_rich, data=td)
}

diff1 <- function(x){
	co <- summary(x)$coeff
	est <- co["naive_rich","Estimate"]
	se <- co["naive_rich","Std. Error"]
	# est + sign(est)*1.96*se
	
	pnorm(1-est, sd=se, lower.tail=FALSE)
}
rbindlist(lapply(naive_msom_mod, mod_smry, pred_name="naive_rich"), idcol='reg')
sapply(naive_msom_mod, diff1)
#' All regions except Northeast U.S. and Newfoundland have estimates that are less than 1, although West Coast U.S. is not signifiantly less than 1. Therefore, for two-thirds of the regions, the MSOM tended to estimate that more undetected species were present when naive estimates were low.  

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Table S1. Summary of Sampling in Each Region
#+ regionalSamplingInfo
sumry_counts <- data_all[reg!="wcann", 
	j = {
		yi <- table(diff(sort(unique(year))))
		yi_form <- paste0(names(yi), paste("(",yi,")",sep=""))
		data.table(
			"N" = trawlData::lu(year),
			"Years" = paste(min(year),max(year),sep=" - "),
			"Frequency" = paste(yi_form,collapse=", "),
			"Sites" = trawlData::lu(stratum), 
			"Tows" = paste0(
				.SD[,trawlData::lu(haulid),by=c("stratum","year")][,max(V1)],
				paste0("(",.SD[,trawlData::lu(haulid),by=c("stratum","year")][,round(mean(V1),1)],")")
				),
			"Species" = trawlData::lu(spp)
		)
	}, 
	by=c('reg')
]
setnames(sumry_counts, 'reg', "Region")
sumry_counts[,Region:=factor(Region, levels=c("ebs","ai","goa","wctri","gmex","sa","neus","shelf","newf"), labels=c("E. Bering Sea","Aleutian Islands","Gulf of Alaska","West Coast US","Gulf of Mexico","Southeast US", "Northeast US", "Scotian Shelf","Newfoundland"))]
setkey(sumry_counts, Region)

sumry_counts[,Org:=c("AFSC","AFSC","AFSC","NMFS","GSMFC","SCDNR","NEFSC","DFO","DFO")]
stargazer(sumry_counts, summary=FALSE, rownames=FALSE, column.sep.width="2pt", digits.extra=0, digits=NA)

# \begin{table}[!htbp] \centering
# %  \caption{}
# %  \label{}
# \begin{tabular}{@{\extracolsep{2pt}} cccccccc}
# \\[-1.8ex]\hline
# \hline \\[-1.8ex]
# Region & N & Years & Frequency & Sites & Tows & Species & Org \\
# \hline \\[-1.8ex]
# E. Bering Sea & $31$ & 1984 - 2014 & 1(30) & $138$ & 4(1.6) & $110$ & $\text{AFSC}^{*}$ \\
# Aleutian Islands & $12$ & 1983 - 2014 & 2(5), 3(4), 4(1), 5(1) & $82$ & 12(2.8) & $54$ & $\text{AFSC}^{*}$ \\
# Gulf of Alaska & $13$ & 1984 - 2013 & 2(7), 3(5) & $89$ & 15(4.4) & $98$ & $\text{AFSC}^{*}$ \\
# West Coast US & $10$ & 1977 - 2004 & 3(9) & $84$ & 91(4.7) & $92$ & $\text{NMFS}^{\dagger}$ \\
# Gulf of Mexico & $17$ & 1984 - 2000 & 1(16) & $39$ & 39(5.7) & $144$ & $\text{GSMFC}^{\ddagger}$ \\
# Southeast US & $25$ & 1990 - 2014 & 1(24) & $24$ & 13(3.9) & $104$ & $\text{SCDNR}^{\S}$ \\
# Northeast US & $32$ & 1982 - 2013 & 1(31) & $100$ & 10(3) & $141$ & $\text{NEFSC}^{\P}$ \\
# Scotian Shelf & $41$ & 1970 - 2010 & 1(40) & $48$ & 11(2.5) & $48$ &$\text{DFO}^{\nabla}$ \\
# Newfoundland & $16$ & 1996 - 2011 & 1(15) & $191$ & 9(2.2) & $72$ & $\text{DFO}^{\nabla}$ \\
# \hline \\[-1.8ex]
# \end{tabular}
# \end{table}

# knitr::kable(sumry_counts, caption="Years = Total number of years sampled, Year Range = the minimum and maximum of years sampled, Year Interval = the number of years elapsed between samples (and the frequency of this interval in parentheses), Sites = the total number of sites in the region, Max. Tows = maximum number of tows per site per year, Average Tows = the average number of tows per site per year, Total Species = the total number of species observed in the region across all sites and years.")


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Table S2. Naive tau  
#+ naiveTau-table, echo=FALSE
kable(rich_naive_trend_kendall, 
	caption="**Table S1.** Naive richness trends in each region. Estimate is Kendall's Tau_b, BH is the Benjamini-Hochberg corrected p-value, and p.value is the original p-value."
)


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Figure S2. Barplot of categories  
#+ Richness-col-ext-barplot, include=TRUE, echo=TRUE, fig.height=3.5, fig.width=3.5, fig.cap="**Figure S2.** Number of species beloning to the categories of both, neither, colonizer, leaver in each region"
# par(cex=1, mar=c(3,1,1,0), ps=8, mgp=c(1,0.25,0), tcl=-0.15, lheight=0.7)
categ_barplot()
# png("~/Desktop/bp.png", width=3.5, height=3.5, res=200, units='in'); categ_barplot(); dev.off()
#+ Richness-col-ext-barplot-table, echo=FALSE
categ_tbl <- t(spp_master[!duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ)])[c(4,1,2,3),]
kable(categ_tbl, caption = "Number of species in each category in each region.")
#' It's the same pattern, whichever way you split it. However, AI is the only region that had more *colonizers* than *both* species. An interesting way to think about some of this is that the average sd in richness was `r comm_master[,stats::sd(reg_rich),by='reg'][,mean(V1)]`, so when the number of *colonizer* or *leaver* species exceed's that region's sd, the impact of those categories, which I consider to be dubious, might start being relevant (though it's not necessarily problematic, nor is this even close to an actual test for the significance of those categories to the trend). EBS and Shelf had significant positive trends in richness and very low numbers in the *colonizer* category. WCTRI and NEWF had similar numbers in the *both* and *colonizer* category.  


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Table S3. Attributing Richness Change to Colonizers (only)  
#+ attribute_richness_categ, echo=FALSE
eval(figure_setup()) # contains pretty_reg names
categ_tbl2 <- spp_master[!duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ, dnn=NULL)]
class(categ_tbl2) <- "matrix" # matrix converts to a data.table better
categ_dt <- data.table(reg=rownames(categ_tbl2),categ_tbl2)

predRich <- bquote(predict(lm(reg_rich~year))) # predicted richness
dels <- comm_master[, # del = delta richness
	list(
		rangeRich=diff(range(reg_rich)), # not using
		delPred=rev(eval(predRich))[1] - eval(predRich)[1] # net change in predicted richness
	),
	by=c('reg')
]

att_categ <- merge(dels, categ_dt)[order(names(pretty_reg))] # merge and sort
att_categ[,c("attr_col","attr_ext"):=list(delPred-colonizer, delPred-leaver)] # new columns

att_categ[,reg:=pretty_reg] # change to pretty region names
setnames(att_categ, old=c("delPred","reg"), c("Richness Change","Region")) # rename columns
att_categ_print <- att_categ[,c("Region", "Richness Change", "colonizer"), with=FALSE] # change order for printing
kable(att_categ_print)
#+ attribute_richness_categ_print, results='asis'
stargazer(att_categ_print, summary=FALSE, rownames=FALSE, column.sep.width="2pt", digits=1, digits.extra=0) # print table
# \begin{table}[!htbp] \centering
#   %\caption{}
#   %\label{}
# \begin{tabular}{@{\extracolsep{2pt}} ccc}
# \\[-1.8ex]\hline
# \hline \\[-1.8ex]
# Region & Richness Change & Colonizing Species \\
# \hline \\[-1.8ex]
# E. Bering Sea & $12.5$ & $2$ \\
# Aleutian Islands & $4.3$ & $8$ \\
# Gulf of Alaska & $17.5$ & $18$ \\
# West Coast US & $18.8$ & $12$ \\
# Gulf of Mexico & $7.4$ & $8$ \\
# Southeast US & $$-$5.0$ & $0$ \\
# Northeast US & $9.5$ & $1$ \\
# Scotian Shelf & $7.6$ & $1$ \\
# Newfoundland & $12.7$ & $10$ \\
# \hline \\[-1.8ex]
# \end{tabular}
# \end{table}
#'   
#'   

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Figure S4. Total Colonizations/ Extinctions vs Georaphic Range  
#+ ceEvents-vs-rangeSize, fig.width=3.5, fig.height=3.5, fig.cap="**Figure S7.** Number of colonizations and extinctions as a function of range size and range density."
ceEventRange("mean_size")
#' Yup, this is definitely a thing. Long-term average range size and range density predict how many colonizations and extinctions a species is likely to have. This will lead nicely into examing how range changes prior to an extinction or after a colonization.  


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Figure S5. Range Size Absence Time Detailed  
#+ rangeSize-timeUntil-detailedFigure, fig.width=6.8, fig.height=6.8
par(mfrow=c(3,3), mgp=c(0.85,0.2,0), ps=8, cex=1, mar=c(1.75,1.5,0.75,0.5), oma=c(0.5,0.5,0.1,0.1), tcl=-0.15)
ur <- names(pretty_reg)
for(r in 1:length(ur)){
	t_reg <- ur[r]
	t_dat <- rangeTimeDT[reg==t_reg]
	t_dat[,spp:=paste0(spp,event)]
	t_dat[type=="pre_ext", c("spp","time"):=list(spp=paste0(spp,type),time=-time)]
	setkey(t_dat, spp, type, time)
	scatterLine(t_dat, x="time", y="size", lineBy="spp", colBy=adjustcolor(pretty_col[t_reg], 1), lwdBy=0.5, type='p', xlab="", ylab="")
	if(t_reg=="ai"){
		t_dat[time <=-10,lines(time, size, col='black')]
	}
	t_dat[,points(mean(time), mean(size), pch=19, col='black', cex=1.2),by=c('time')]
	t_dat[,points(mean(time), mean(size), bg=pretty_col[t_reg], pch=21, col='white'),by=c('time')]

	mtext(pretty_reg[t_reg], line=0, side=3, adj=0.1, font=2)
}
mtext("Time to Event", side=1, line=-0.5, outer=TRUE)
mtext("Range Size", side=2, line=-0.5, outer=TRUE)



#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Extra Stuff (not main text, not supplement)  

#' ##Figure. Time series of colonizations and extinctions
#+ Richness-col-ext-ts, include=TRUE, echo=TRUE, fig.height=6.8, fig.width=6.8, fig.cap="**Figure NotIncluded.** Number of colonizations (blue) and extinctions (red) over time in each region."
col_ext_ts()

#' In most regions the differences in colonization and extinction numbers are similar over time. The most obvious exceptions are for the 3 regions that showed large initial spikes in richness; the GOA, GMEX, and AI regions initially have much larger numbers of colonizers than leavers, but this number shrinks rapidly until the two rates are ~equal.  
#'   
#' For the regions with significant positive slopes, there is no visually obvious increase in colonizations relative to extinctions over time. Because the colonization and extinction numbers tend to track each other over the long-term, it it would be difficult to attribute the long-term changes in richness to a change in just colonization or extinction rates.  
#'   
#' **Manuscript paragraph:**  
#' A time series of richness can be decomposed into the colonizations and extinctions of individual species over time. We categorized species according to the following colonization extinction patterns: present in all years = neither (536 species), colonized and went extinct = both (263 species),  initially absent but present every year after its colonization = colonizer (61 species), initially present but absent every year after its extinction = leaver (4 species). Most regions had the same overall ranking (neither > both > colonizer > leaver), except in the Northeast US where both was the most common and neither second, and in the Aleutian Islands where colonizer was the second most common and both third (**Figure S2**). In general, changes in richness were not due to permanent departures or introductions of species to the region. Furthermore, colonization and extinction rates did not become more dissimilar over time for any region (**Figure NotIncluded**). Colonizations were initially greater than extinctions in Aleutian Islands, Gulf of Alaska, and Gulf of Mexico, but this difference disappeared in the latter portion of these time series, as evidenced by these regions’ initially rapid increase in richness that later plateaued. The other regions did not show strong trends in the difference between colonizations and extinctions over time, making it difficult to attribute the long-term trends in richness to changes in just colonization or just extinction rates.  


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Figure. Changes in Geographic Density before Extinction and after Colonization  
#+ rangeDensity_ColExt, fig.width=6, fig.height=3, fig.cap="**Figure 4b.** Geographic range density vs years until extinction (A) and years after colonization (B). For visualization purposes, range density is averaged across species for each unique value on each axis, and a linear model fit through this average. Statistics in main text use unaggregated data. The horizontal axes were formulated as time until (since) the nearest upcoming (previous) absence. Because range density must be zero when either horizontal axis has a value of zero, points at (0,0) were excluded from figures and analyses."
rangeSize_absenceTime("rangeDensity")

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Figure. Novel Colonizations Per Year
#+ novelCol-ts, include=TRUE, echo=TRUE, fig.height=6.8, fig.width=6.8
fp <- spp_master[present==1,.SD[which.min((year))] ,by=c('reg','spp')]
fp2 <- fp[,list(new_col=sum(col)),keyby=c('reg','year')]
ry_skele <- unique(spp_master[,list(reg,year)])
fp3 <- merge(ry_skele, fp2, all=TRUE, by=c('reg','year'))
fp3[is.na(new_col), new_col:=0]
par(mfrow=c(3,3))
fp3[,plot(year, new_col, type='o', main=reg[1]),by=c('reg')]
#' Reveals similar information as the number of colonizations per year, except it isolates the number of colonizations by species that have never been seen in the region before. This is the rate at which new species are introduced to the region. An important caveat here is that the trends here are probably biased downward because we onlly included species that were observed in at leat 10 year-sites in order to avoid the inclusion of ultra-rare species. These exclusions reduce the odds that species that colonized a limited number of sites towards the end of the time series met the threshold of 10, especially for regions with relatively few strata.  
#'   
#' E.g., say that a species colonized Southeast U.S., which has 24 sites, in the last year of sampling --- for this species to be included, it would have had to occupy 42% of the region (10 sites) in its first year present. Across all regions, the average colonizing species occupies only 5% of strata in its first year, and for the 3 species that colonized the Southeast U.S., they colonized only 4.2% of sites, a full order of magnitude lower than what was needed. Depending on the expansion rate, a species would need to first colonize several years before the end of the time series to be included. It's also worth noting that no region has a new colonizing species towards the end of the time series --- this makes a similar point.  
#'   
#' Because of the way we eliminated species from our data set, our analysis was extremely conservative with regard to detecting changes in richness, particularly if those changes were going to be the result of changes in the rate of novel introductions (right now I'm a little less clear on what sources of species declines might be affected, other than to say that generally fewer species would have been included due to only being present at the beginning [affects declines more] or end [affects increases more] of the time series). Regions with short time series and few sites were especially prone to having muted changes in richness.  


#'   
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'  
#' ##Figure. Relationship between range size and range density
#+ geo-range-densVsize, fig.width=3.5, fig.height=5.5, fig.cap="**Figure S6**. Range density versus range size. In panel A each point is a species-region-year combination. In panel B, each point is a region-year. Range size is the proportion of sites occupied, range density the tows in occupied sites. The community metrics in B is calculated by take each species' long-term average from A, then taking the average across all species present in the community in a given year. Fitted lines in A are from a loess fit."
rangeSizeDens()
#' First posterity, I'll show a figure relating range size and range density. Then I'll present the figure of how size can 'predict' richness. Finally, I'll explore models more formally expressing the relationship between richness and size.   
#'   
#' These plots show that there is a strong relationship between range size and range density. Interestingly, in **Figure S6B**, the cross-region relationship is negative (if each color had 1 pt), whereas the within-region relationship is positive.  
#'   
#' The positive relationship between size and density is not surprising. My interpretation of density is ~population size. Often population size is correlated with range size. I think this is a standard result, but I need to double-check.  

#' ##Figure. Range size vs absence time as individual regressions
#+ rangeSize_time_sepRegs, fig.width=5, fig.height=6, fig.cap="**Exploration Figure** histograms of separate regressions of size ~ time; this is for each run-up to an extinction and reach follow-up to a colonization. Trying to understand how regularly the pattern might be observed. Hard to answer because adjusting the restriction on number of events in the run-up or follow-up (nTime) greatly affects the proportion that are significant."
rangeTimeDT[,nTime:=length(time),by=c("reg","type","event","spp")]
o <- rangeTimeDT[nTime>=3,j={
	getEPR(lm(size~time))
	} ,by=c("reg","type","event","nTime","spp")
]
rangeTime_signif <- o[,list(propSignificant=(sum(Pr<0.05)/sum(!is.na(Pr))),n=sum(!is.na(Pr))),by=c("reg","type")]

# histograms of estimates, p-value, and r-squared split by pre-extinction and post-colonization
par(mfrow=c(4,2), mar=c(2,2,0.5,0.5), cex=1, ps=10, mgp=c(1,0.1,0), tcl=-0.1)
o[,j={hist(Estimate, main=type[1]);NULL},by='type']
o[,j={hist(Pr, main=type[1]);NULL},by='type']
o[,j={hist(Rsquared, main=type[1]);NULL},by='type']

# p value vs number of years in series
setorder(o, nTime)
# o[,j={plot(nTime,Pr); lines(spline(Pr~nTime),col='red',lwd=2)}]
# abline(h=0.05, lwd=1.5)
# abline(h=0.05, col='white', lwd=0.5)
o[,mean(Pr),by='nTime'][,plot(nTime, V1)]; 
abline(h=0.05, col='blue', lty=2)

# slope vs number of years in series
# o[,j={plot(nTime,Estimate); lines(spline(Estimate~nTime),col='red',lwd=2)}]
# abline(h=0, lwd=1.5)
# abline(h=0, col='white', lwd=0.5)
o[,mean(Rsquared),by='nTime'][,plot(nTime, V1)]; 

#+ rangeSize_time_sepRegs-table1, echo=FALSE
kable(rangeTime_signif, caption="Significant trends in range size before ext/ after col. Sample size is number of 'stretches'.")
kable(rangeTime_signif[,list(proportionSignificant=mean(propSignificant)),by=c("type")])
#' 
#' The mixed effect models show a strong tendency for range size to be smaller near an absence. I figured that this relationship wouldn't be apparent for all cases, which ended up being the result. However, the significance of these relationships depends on the number of years for each event stretch. So sample size matters. On one hand, this could hold implications for the rate of decline, so excluding short stretches might also exclude more cases with abrupt declines/ increases (although, these may well be preceded by long durations of stability, so not all abrupt declines/ increases would be excluded necessarily). However, the slope estimate doesnt change a whole lot across sample size: from 0 to ~10 years, increasing sample size (duration) also increases the slope (going from slightly negative to pretty consistently positive). Then it levels out, and doesn't continue increasing. So it is *not* the case that longer stretches are long because they have more gradual slopes, necessarily. Actually, in the range of data for which there does appear to be a trend (0-10 samples), the slope becomes larger which is the opposite of "longer stretches result from more gradual processes" notion (thinking in terms of extinction).  
#' 
#' Ultimately, I think the mixed effects model is far more appropriate so long as conclusions are cauched in general patterns, not in being able to detect the pattern for any single species. While it could have been interesting to present the many-regression results too, I think it is a separate analysis to consider how often this rule would apply to individual species. Such an analysis would benefit from understanding when the rule does and does not apply, not just the apparent frequency of relevance.  
#'  


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Figure. Long-term Range Trends for Core, Leaving, Colonizing, Both
#+ rare-more-common, fig.width=7, fig.height=7
blah <- spp_master[present==1, j={
	list(propSlope=summary(lm(range_size_mu~year))$coeff[2,1], mu_propStrat=mean(range_size_mu), ce_categ=ce_categ[1])
} ,by=c("reg","spp")]
par(mfrow=c(3,3))
ureg <- blah[,unique(reg)]
nreg <- length(ureg)
for(r in 1:nreg){
	blah[reg==ureg[r],j={boxplot(propSlope~ce_categ, main=reg[1], outline=FALSE);abline(h=0);NULL}]
}
#'   
#' I wanted to determine if species in the region are becoming more widespread over time. In other words, when the species is present, what fraction of the region does it occupy? Does that fraction change over time? This slope is the y-axis for the boxplots. It indicates the long-term slope of the proportion of sites occupied (excluding any year when that proportion was 0). The x-axis represents the categories of colonization-extinction patterns. Species in the "both" category experienced both colonization and extinction events at some point in the time series.  
#' 
#' Other figures and analyses suggest that species in the "both" category play an important role in *increases* in species richness. In other words, the coming-and-going of these organisms is not neutral in the long term. This conjecture is supported by 1) the magnitude of increase for regions with increasing richness cannot be explained by "colonizing" species alone; 2) geographically constrained species have more colonizations and extinctions that widespread species.  
#' 
#' Looking at the different box categories, if historically rare species (i.e., those that were in the both category and had most of the colonizations/ extinctions) were appearing more consistently in the time series because they were becoming more geogrpahically widespread (and therefore less likely to go extinct), then we'd expect to see that the "both" category to have more-positive slopes in their annual range sizes.  
#' 
#' That's kinda? what I see. It's definitely not a very strong relationship. One issue with this analysis is that you wouldn't expect any given species in the "both" category to become more prevalent. In fact, it should only be a small fraction of those species that are becomign more widespread and thus contributing to long-term increases in species richness.  


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Figure. Time Series of Local and Regional Richness
#+ alpha-gamma-TimeSeries, fig.width=7, fig.height=7
localR <- data_all[,list(lR=length(unique(spp))),by=c("reg","stratum","year")]
bothR <- merge(localR, comm_master[,list(reg,year,naive_rich,reg_rich)],by=c("reg","year"))
par(mfrow=c(3,3), cex=1, ps=8, mar=c(2,2,2,2), mgp=c(1,0.2,0), tcl=-0.2)
ureg <- bothR[,unique(reg)]
nreg <- length(ureg)
for(r in 1:nreg){
	comm_master[reg==ureg[r],j={
		plot(year, get(lR), type='l', ylab=lR)
		mtext(reg[1], side=3, line=0, adj=0.1, font=2)
		par(new=TRUE)
		plot(year, get(rR), type='l', col='red', xaxt='n', yaxt='n', xlab='',ylab='')
		axis(side=4, col='red')
		mtext(rR, side=4, line=1, col='red')
		NULL
	}]
}
mtext("Observed Regional (red) and Mean Local (black) Richness", side=3, line=-0.75, outer=TRUE, font=2)
#' In general, local richness and regional richness are similar. There are differences, possibly in some cases the regional slope might be significant whereas local would not be (I haven't checked, though). However, the trends aren't in opposite directions at the two scales.  


#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Comparisons of Alpha, Beta, Gamma Diversity
#' ###Figure. Alpha, beta, and gamma diversity (Naive)
#+ abgDiversity-Naive, fig.width=3.5, fig.height=7, caption="Alpha, beta, gamma diversity."
eval(figure_setup())
par(mfrow=c(3,1), mar=c(1.75,1.5,0.25,0.25),mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
ureg <- comm_master[,unique(reg)]
nreg <- length(ureg)
ptCol <- bquote(adjustcolor(pretty_col[reg],0.5))
scatterLine(Data=comm_master, x="naive_rich", y="beta_div_obs", lineBy="reg", col=ptCol, pch=20)
scatterLine(Data=comm_master, x="local_rich_obs", y="beta_div_obs", lineBy="reg", col=ptCol, pch=20)
comm_master[,legend("topright",ncol=2,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.01, -0.03), bty='n', x.intersp=0.15, y.intersp=0.65)]
scatterLine(Data=comm_master, x="naive_rich", y="local_rich_obs", lineBy="reg", col=ptCol, pch=20)

#' ###Figure. Alpha, beta, and gamma diversity (MSOM)
#+ abgDiversity-MSOM, fig.width=3.5, fig.height=7, caption="Alpha, beta, gamma diversity."
eval(figure_setup())
par(mfrow=c(3,1), mar=c(1.75,1.5,0.25,0.25),mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
ureg <- comm_master[,unique(reg)]
nreg <- length(ureg)
ptCol <- bquote(adjustcolor(pretty_col[reg],0.5))
scatterLine(Data=comm_master, x="reg_rich", y="beta_div_mu", lineBy="reg", col=ptCol, pch=20)
scatterLine(Data=comm_master, x="local_rich", y="beta_div_mu", lineBy="reg", col=ptCol, pch=20)
comm_master[,legend("topright",ncol=2,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.01, -0.03), bty='n', x.intersp=0.15, y.intersp=0.65)]
scatterLine(Data=comm_master, x="reg_rich", y="local_rich", lineBy="reg", col=ptCol, pch=20)


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Richness-range size, varying long-term and annual, all spp and transients only
#+ annualLongTerm-allSppNoNeitherSpp, fig.width=6.5, fig.height=6.5
noNeither <- spp_master[,j={
	td <- .SD[ce_categ!="neither" & present==1]
	td1.0 <- td[,list(propStrata_noNeither_ltAvg=mean(range_size_mu)),keyby=c("reg","spp")]
	setkey(td, reg, spp)
	td1.1 <- td[td1.0]
	td1.2 <- td1.1[,list(propStrata_noNeither_avg_ltAvg=mean(propStrata_noNeither_ltAvg)),by=c("reg","year")]
	td2 <- td[,list(propStrata_noNeither_avg=mean(range_size_mu)),by=c("reg","year")]
	td_out <- merge(td1.2, td2, by=c('reg','year'))
}]
cmNN <- merge(comm_master, noNeither)

par(mfrow=c(2,2))
ureg <- comm_master[,unique(reg)]
nreg <- length(ureg)
ptCol <- bquote(adjustcolor(pretty_col[reg],0.5))
scatterLine(Data=cmNN, x="range_size_mu_avg_ltAvg", y="reg_rich", lineBy="reg", col=ptCol, pch=20)
scatterLine(Data=cmNN, x="range_size_mu_avg", y="reg_rich", lineBy="reg", col=ptCol, pch=20)
scatterLine(Data=cmNN, x="propStrata_noNeither_avg_ltAvg", y="reg_rich", lineBy="reg", col=ptCol, pch=20)
scatterLine(Data=cmNN, x="propStrata_noNeither_avg", y="reg_rich", lineBy="reg", col=ptCol, pch=20)


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Range Size and Extinction Probability
#' ###Statistics for Range Size Year before Extinction
#+ rangeProbExt-modelFit, results='markup'
std_range_t1 <- spp_master[ce_categ!='neither' & ce_categ!="colonizer", j={
	
	std_range <- .SD[,list(year, ext, std_range=c((get(rSMet_base))), ext_dist=ext_dist),by=c("spp")]
	setkey(std_range, spp, year)
	std_range_t1 <- std_range[present==1,list(year=year[-1], ext=ext[-1], std_range=std_range[-1], std_range_t1=head(std_range, -1), ext_dist=ext_dist[-1]),by=c("spp")]
	
	std_range_t1[,nspp:=length(unique(spp)),by=c('year')]
	std_range_t1
	
},by="reg"]
# print(std_range_t1[,{print(summary(glm(ext~std_range, family='binomial')));NULL},by='reg'])
# print(std_range_t1[,{print(summary(glmer(ext~std_range+(std_range|reg), family='binomial')));NULL}])

# Set up -- region names and empty named lists
ur <- std_range_t1[,unique(reg)]
rangeExtProb1 <- structure(vector("list",length(ur)), .Names=ur)
rangeExtProb2 <- structure(vector("list",length(ur)), .Names=ur)
rangeExtProb3 <- structure(vector("list",length(ur)), .Names=ur)

# Fit 3 models for each region
for(r in 1:length(ur)){
	t_reg <- ur[r]
	t_dat <- std_range_t1[reg==t_reg]
	rangeExtProb1[[t_reg]] <- (glm(ext~std_range, family='binomial', data=t_dat))
	rangeExtProb2[[t_reg]] <- summary(glm(ext~std_range+nspp, family='binomial', data=t_dat))
	rangeExtProb3[[t_reg]] <- summary(glm(ext~std_range*nspp, family='binomial', data=t_dat))
}
rangeExtProb1_smry <- (smry_modList(rangeExtProb1, pred_name="std_range"))[,reg:=names(rangeExtProb1)]
rangeExtProb1_smry[,BH:=p.adjust(p.value, method='BH')]

#' ###Figure for Range Size Year before Extinction
#+ rangeProbExt-modelFit-fig, fig.width=6, fig.height=6
par(mfrow=c(3,3), mar=c(2.5,2.5,1.5,0.5), ps=10, mgp=c(0.75,0.15,0), tcl=-0.15)
for(r in 1:length(ur)){
	std_range_t1[reg==ur[r],plot(std_range, ext)]
	new_r <- std_range_t1[reg==ur[r],sort(std_range)]
	lines(new_r, predict(rangeExtProb1[[ur[r]]], newdata=data.frame(std_range=new_r), type='response'))
	mtext(pretty_reg[ur[r]], side=3, line=0, adj=0.1, font=2)
}


#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ###Statistics for Range Size Predicting Years until Extinction
#+ rangeProbExt-modelFit2, results='markup'
make_timeRange <- function(densName=c("propTow_occ"), sizeName=c("range_size_samp","range_size_mu","propStrata")){
	rangeTimeDT <-  spp_master[stretch_type=="pre_ext" & !is.na(stretch_type) & propStrata!=0]
	rangeTimeDT <- rangeTimeDT[,list(
		reg=reg, 
		event=as.character(event_year), 
		spp=spp, 
		type=as.character(stretch_type),
		time=ext_dist, 
		size=get(sizeName),
		density=get(densName)
	)]
	return(rangeTimeDT)
}

timeRangeDT <- make_timeRange()
ur <- timeRangeDT[,unique(reg)]
sTime_reg_mods_Inv <- structure(vector("list",length(ur)), .Names=ur)
for(r in 1:length(ur)){
	t_reg <- ur[r]
	t_dat <- timeRangeDT[reg==t_reg]
	sTime_reg_mods_Inv[[t_reg]] <- lme4::lmer(time ~ size + (size|spp), data=t_dat)
}
sTime_reg_smry_Inv <- smry_modList2(sTime_reg_mods_Inv)
kable(sTime_reg_smry_Inv)

#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Beta Div vs CRI
#+ betaDiv-CRI-plot, fig.width=3.5, fig.height=3.5
par(mfrow=c(3,3), mar=c(1.85,1.85, 0.5,0.5), mgp=c(0.65, 0.15, 0), oma=c(0.5, 0.5, 0.25, 0.1), tcl=-0.15, ps=8, cex=1)
ur <- names(pretty_reg)[names(pretty_reg)%in%comm_master[,unique(reg)]]
nr <- length(ur)
for(r in 1:nr){
	comm_master[reg==ur[r],j={
		plot(beta_div_obs, propStrata_avg_ltAvg, xlab='', ylab='')
		mod <- lm(propStrata_avg_ltAvg~beta_div_obs)
		abline(mod)
	}]
	mtext(pretty_reg[ur[r]], side=3, adj=0.1, font=2)
}





#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Boxplot Time Series of Range Size, Color is Richness of Transient Species (Resampled)
#+ rangeTS-colorAlpha-samp, fig.width=6.5, fig.height=6.5, caption="Range size and richness of transient species. Boxplots represent the cross-species distribution of range sizes. The color of the boxplot indicates the number of  transient species present in that year."
eval(figure_setup())
eval(neitherPres_subColor())
boxRange_colRich(range_type="range_size_samp")



#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Richness Trend vs Transient Range Trend
#+ richTrend-vs-transientRangeTrend, fig.width=5, fig.height=5
dat <- merge(transExpansionStats, comm_master[,coef(lm(reg_rich~year))[2], by='reg'], by='reg')
dat[,summary(lm(V1~year))]
dat[,plot(year, V1, col=pretty_col[reg], pch=19, xlab="Rate of Range Expansion by Transient Species", ylab="Rate of Change in Species Richness")]
#' The significance of this trend depends on Southeast US, though. I'm not sure I want to get into defending a plot based on an outlier, and the fact that I didn't account for uncertainty in the x-axis.  


#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#+ systemSettings, results='markup'
Sys.time()
sessionInfo()