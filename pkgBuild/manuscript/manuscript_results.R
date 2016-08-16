#' ---
#' title: "Shifting distributions and long-term changes in marine assemblage richness"
#' author: "Ryan Batt"
#' date: "2016-08-15"
#' abstract: |
#'   These are the results and figures.
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 4
#'     fig_caption: true
#'     theme: "readable"
#'     template: default
#'   pdf_document:
#'     toc: true
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
o_f <- c("html_document", "pdf_document")[1]

# render!
# rmarkdown::render(
# 	"~/Documents/School&Work/pinskyPost/trawlDiversity/pkgBuild/manuscript/manuscript_results.R",
# 	output_format=o_f,
# 	output_dir='~/Documents/School&Work/pinskyPost/trawlDiversity/pkgBuild/manuscript',
# )

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
# source("../manuscript/fig_tbl_number.R")
# eval(fig_tbl_number())


# ============
# = Richness =
# ============
#+ Richness-basic, include=TRUE, echo=TRUE, eval=TRUE
#' \FloatBarrier   
#' 
#' ***  
#' 
#' #Results
#' ##Species Richness
#' ###Richness Summary
#' Give a basic summary of species richness stats. For example, the regions with the lowest and highest long-term averages in species richness are:  
kable(comm_master[,mean(reg_rich), by='reg'][reg%in%c("gmex","shelf")])
#' Can also measure long-term variability of the different regions. Here is the standard deviation for each region:
kable(comm_master[,stats::sd(reg_rich),by='reg'])
#' And here is the cross-region average of those long-term standard deviations:
kable(comm_master[,stats::sd(reg_rich),by='reg'][,mean(V1)])

#+ Richness-ts-fig.cap, echo=FALSE
richness.ts.fig.cap <- "**Figure 1.** Time series of MSOM estimates of region richness. Each point is the posterior mean of regional richness in a year. Lines indicate long-term trends from fitted values of linear regression models predicting richness from time."
#' ####Figure 1. MSOM richness time series
#+ Richness-ts-fig, fig.height=5, fig.width=3.5, fig.cap=richness.ts.fig.cap, ecal=TRUE, echo=TRUE
richness_ts()

#+ Richness-msom-naive-scatter-cap, echo=FALSE
msom_naive_scatter_cap <- "**Figure S1.** MSOM richness vs naive richness"
#' ####Figure S1. MSOM - naive scatter
#+ Richness-msom-naive-scatter-fig, fig.height=3.5, fig.width=3.5, fig.cap=msom_naive_scatter_cap
naive_msom_scatter()
#' MSOM richness and Naive richness are pretty similar. MSOM richness is probably more accurate, or is at least more conservative b/c it has fewer significant trends. But their similarity should help justify using observed presences/ absences in other analyses.  
#'   
#' **Manuscript paragraph:**  
#' MSOM estimates of richness were greater than Naive estimates, but the two methods produced similar temporal dynamics in richness (Figure S1). Henceforth, we report MSOM richness estimates. The greatest long-term average richness was the Gulf of Mexico (142.3), lowest in the Scotian Shelf (46.3). The inter-region average of long-term standard deviations of richness was 5.3; the Aleutian Islands showed the lowest variability (sd = 2.4), and the Gulf of Alaska was the most variable (sd = 8.5).  
#' 
#' \FloatBarrier   
#' 
#' ***  
#' 
#+ Richness-trend, include=TRUE, echo=TRUE, results='asis'
#' ###Richness Trend
#' Examine trends for naive richness:  
load("../../pkgBuild/results/rich_naive_trend_kendall.RData")
rich_naive_trend_kendall[reg!="wcann",BH:=p.adjust(taup, method='BH')]
rich_naive_trend_kendall <- rich_naive_trend_kendall[reg!="wcann",list(reg=reg, estimate=tau, BH=BH, p.value=taup)]
#' ####Table S1. Naive tau
kable(rich_naive_trend_kendall)
#' In the Naive estimates, `r rich_naive_trend_kendall[reg!='wcann', sum(p.value<=0.05)]` regions had significant $\tau_b$.  

#' Examine trends for MSOM richness:  
load("../../pkgBuild/results/rich_trend_kendall.RData")
rich_trend_kendall[reg!="wcann",BH:=p.adjust(pvalue, method="BH")]
rich_trend_kendall <- rich_trend_kendall[reg!="wcann",list(reg=reg, estimate=tau, BH=BH, p.value=pvalue)]
#' ####Table 2. MSOM tau
kable(rich_trend_kendall)
#' In the MSOM estimates, `r rich_trend_kendall[reg!='wcann', sum(p.value<=0.05)]` regions had significant $\tau_{b}$.  
#'   
#' Overall, ~half the regions show positive trends. No regions show significant negative trends (although SEUS is close, depending on the analysis). It's a bit surprising how much removing the autocorrelation seemed to impact some of the trends (AI in particular, I think).  
#'   
#' **Manuscript paragraph:**  
#' Estimated slopes for long-term trends (Kendall’s τb) in richness were positive for most regions. For Naïve estimates, all the seven positive τ were also significantly different from 0, whereas the two negative τb, Aleutian Islands and Southeast US, were not (Table S2). Long-term trends in MSOM estimates of richness had the same sign as Naïve trends, except for Aleutian Islands, but now the only significant τb were the following four regions: Eastern Bering Sea (τb = 0.41), West Coast US (τb = 0.61), Scotian Shelf (τb = 0.45), and Newfoundland (τb = 0.72) (Table 1). Richness trends were not significant in the Gulf of Mexico, Gulf of Alaska, and Northeast US, Aleutian Islands, Southeast US.   
#'   
#' \FloatBarrier   
#' 
#' ***  
#' ##Richness and Geographic Range
#' Looking for spatial footprint that will predict richness. Here I characterize the geographic distribution of a species in two ways: its density and its size. This is tricky for richness because richness is a community level trait, but I've characterized these two aspects of geographic distribution at the species level. In fact, these "species level" attributes are also potentially dyanmic -- a species need not have fixed range density or size. So there are two levels of averaging -- for each species take its long-term average value, then for each community take the average across species.  
#' 
#' First I'll show a figure relating range size and range density. Then I'll present the figure of how size and density can predict richness. Finally, I'll explore models more formally expressing the relationship between richness and size/ density.
#' 
#'   
#' ####Figure S2
#+ geo-range-densVsize, fig.width=3.5, fig.height=5.5, fig.cap="**Figure S2**. Range density versus range size. In panel A each point is a species-region-year combination. In panel B, each point is a region-year. Range size is the proportion of sites occupied, range density the tows in occupied sites. The community metrics in B is calculated by take each species' long-term average from A, then taking the average across all species present in the community in a given year. Fitted lines in A are from a loess fit."
rangeSizeDens()
#' These plots show that there is a strong relationship between range size and range density. Interestingly, in Figure S2B, the cross-region relationship is negative (if each color had 1 pt), whereas the within-region relationship is positive.  
#'   
#' The positive relationship between size and density is not surprising. My interpretation of density is ~population size. Often population size is correlated with range size. I think this is a standard result, but I need to double-check.
#'   
#' ####Figure 2
#+ rich-geo-range, fig.width=3.5, fig.height=5.5, fig.cap="**Figure 2.** Species richness vs A) geographic range size, and B) geographic range density. Both metrics are based on each species' long-term average of a statistic; range size is the proportion of sites occupied, range density is the proportion of tows in occupied sites. Solid lines are linear regressions with MSOM richness as the response and the horizontal axis and an intercept as the predictors."
rich_geoRange()
#' Both range size and range density are pretty good predictors of species richness. I think I had originally missed the range size relationship b/c I hadn't done the same aggregating procedure. The interpretation I have is that richness is highest when you have a bunch of rare species.  
#'   
#' ###Richness and Range Density
#' The goal here is to see if species richness is predicted by the typical range density of community's constituent species. First I'll run different types of models just to explore whether this is true, in general (across regions). Then I'll drill in to each region individually to answer the same question.  
#'   
#+ rich-rangeDensity
# This is a function that'll help summarize model fit and coefficients and parameter significance:  
mod_smry <- function(m, pred_name=c("density","size")){
	pred_name <- match.arg(pred_name)
	sc <- sem.coefs(m)
	mod_call <- switch(class(m), lmerMod=m@call, lm=m$call)
	mod_call <- as.character(mod_call)[2]
	out <- cbind(
		mod_call = mod_call,
		sc[sc[,"predictor"]==pred_name,],
		sem.model.fits(m)
	)
	out[,c("Class","N","mod_call","predictor","estimate","std.error","p.value","Marginal","Conditional")]
}

# Make a data set that is useful for these regressions (short variable names, etc):  
range_reg <- comm_master[,list(reg, year, rich=reg_rich, density=propTow_occ_avg, size=propStrata_avg_ltAvg)]

# Fit different models to the whole data set
rich_dens_mods <- list()
rich_dens_mods[[1]] <- lm(rich ~ 1 + density, data=range_reg)
rich_dens_mods[[2]] <- lm(rich ~ 1 + density*reg, data=range_reg)
rich_dens_mods[[3]] <- lmer(rich ~ 1 + density + (1|reg), data=range_reg)
rich_dens_mods[[4]] <- lmer(rich ~ 1 + density + (1+density|reg), data=range_reg)
rich_dens_smry <- rbindlist(lapply(rich_dens_mods, mod_smry, pred_name="density"))

# Fit same model to each region separately 
rich_dens_reg_mods <- list()
ur <- range_reg[,unique(reg)]
for(r in 1:length(ur)){
	rich_dens_reg_mods[[r]] <- lm(rich ~ 1 + density, data=range_reg[reg==ur[r]])
}
rich_dens_reg_smry <- data.table(reg=ur, rbindlist(lapply(rich_dens_reg_mods, mod_smry, pred_name="density")))

# unlist(lapply(rich_dens_mods, stargazer, type='html'))
stargazer(rich_dens_mods[[1]], rich_dens_mods[[2]], rich_dens_mods[[3]], rich_dens_mods[[4]], type='html')
kable(rich_dens_smry) # different kinds of models
kable(rich_dens_reg_smry) # same model applied to each reg sep
kable(as.data.frame(as.list( # avg of above
	sapply(rich_dens_reg_smry, function(x)suppressWarnings(mean(x)))
))) 
#' All models are pretty good predictors. Well, the most basic model kinda sucks I guess. It needs to account for some of the between-region variation.  
#'   
#' ###Richness and Range Size
#' Same as above, but now let's look at range size as a predictor of species richness.  
#'   
#' "Is range size a good predictor of species richness?"  
#+ rich-rangeSize
rich_size_mods <- list()
rich_size_mods[[1]] <- lm(rich ~ 1 + size, data=range_reg)
rich_size_mods[[2]] <- lm(rich ~ 1 + size*reg, data=range_reg)
rich_size_mods[[3]] <- lmer(rich ~ 1 + size + (1|reg), data=range_reg)
rich_size_mods[[4]] <- lmer(rich ~ 1 + size + (1+size|reg), data=range_reg)
rich_size_smry <- rbindlist(lapply(rich_size_mods, mod_smry, pred_name="size"))

rich_size_reg_mods <- list()
ur <- range_reg[,unique(reg)]
for(r in 1:length(ur)){
	rich_size_reg_mods[[r]] <- lm(rich ~ 1 + size, data=range_reg[reg==ur[r]])
}
rich_size_reg_smry <- data.table(reg=ur, rbindlist(lapply(rich_size_reg_mods, mod_smry, pred_name="size")))

# unlist(lapply(rich_size_mods, stargazer, type='html'))
stargazer(rich_size_mods[[1]], rich_size_mods[[2]], rich_size_mods[[3]], rich_size_mods[[4]], type='html')
kable(rich_size_smry)
kable(rich_size_reg_smry)
kable(as.data.frame(as.list(sapply(rich_size_reg_smry, function(x)suppressWarnings(mean(x))))))

#' ###Predicting Richness: Range Size or Density?
#' As far as picking one or the other, it doesn't end up mattering much. Range size is a lot better than density in NEUS, and density outperforms size in AI. Otherwise, size as a slight edge over density on average, although both predictors are significant in all regions.  
#'   
#+ rich-range-compareSizeDens
rich_ds_smry <- rbind(rich_dens_reg_smry, rich_size_reg_smry)
setkey(rich_ds_smry, reg, Class, predictor)
kable(rich_ds_smry)
fight <- rich_ds_smry[,list(marginal_density_minus_size=.SD[predictor=="density", Marginal] - .SD[predictor=="size", Marginal]),by='reg']
kable(data.table(fight, mu=fight[,mean(marginal_density_minus_size)]))
#' ###Conclusion for Richness and Geographic Range
#' The two metrics of geographic range are well correlated. Furthermore, richness can be predicted pretty well using regressions with either as a predictor. There are large differences among regions, though. This is probably because richness is not readily comparable among most regions. Regions vary mostly in their intercept values, and they have fairly similar slopes (though they are not identical, and model fits improve when allowing slopes to vary among regions; it's just that the improvement is small compared to allowing intercepts to vary among regions).  
#'  
#' The interpretation of the result that geographic distribution predicts species richness is likely associated with species rarity. When the average range density or range size of a community is low, it means it has a lot of species that are rare (at either spatial scale). It's these rare species that come and go, and form the dynamics of richness that we observe. When that dynamical value is high, it implies that an above-average number of the dynamic species are present. Because those species are transient (dynamic), they are also probably rare.  
#'   
#'   
#' \FloatBarrier   
#' 
#' ***  
#' 
#' ##Colonization and Extinction  
#' ###Colonization and Extinction Summary
#' ####Figure S3.  Barplot of categories
#+ Richness-col-ext-barplot, include=TRUE, echo=TRUE, fig.height=5, fig.width=3.5, fig.cap="**Figure S3.** Number of species beloning to the categories of both, neither, colonizer, leaver in each region"
categ_barplot()
# And here's the table behind the barplot:  
kable(t(
	spp_master[!duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ)]
)[c(4,1,2,3),])
# Or, to see the most common categories overall:  
kable(data.frame(as.list(
	rowSums(spp_master[!duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ)][,c(4,1,2,3)])
)))
# Does that pattern change if we just look at regions with positive slopes? 
pos_reg <- c("ebs","newf","shelf","wctri")
kable(data.frame(as.list(
	rowSums(t(
		spp_master[reg%in%pos_reg & !duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ)]
	)[c(4,1,2,3),])
)))
kable(data.frame(as.list(
	rowSums(t(
		spp_master[!reg%in%pos_reg & !duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ)]
	)[c(4,1,2,3),])
)))
#' It's the same pattern, whichever way you split it. However, AI is the only region that had more *colonizers* than *both* species. An interesting way to think about some of this is that the average sd in richness was `r comm_master[,stats::sd(reg_rich),by='reg'][,mean(V1)]`, so when the number of *colonizer* or *leaver* species exceed's that region's sd, the impact of those categories, which I consider to be dubious, might start being relevant (though it's not necessarily problematic, nor is this even close to an actual test for the significance of those categories to the trend). EBS and Shelf had significant positive trends in richness and very low numbers in the *colonizer* category. WCTRI and NEWF had similar numbers in the *both* and *colonizer* category.  
#'   
#' ####Figure S4. Time series of colonizations and extinctions
#+ Richness-col-ext-ts, include=TRUE, echo=TRUE, fig.height=3.5, fig.width=3.5, fig.cap="**Figure S4.** Number of colonizations (blue) and extinctions (red) over time in each region."
col_ext_ts()
#' In most regions the differences in colonization and extinction numbers are similar over time. The most obvious exceptions are for the 3 regions that showed large initial spikes in richness; the GOA, GMEX, and AI regions initially have much larger numbers of colonizers than leavers, but this number shrinks rapidly until the two rates are ~equal.  
#'   
#' For the regions with significant positive slopes, there is no visually obvious increase in colonizations relative to extinctions over time. Because the colonization and extinction numbers tend to track each other over the long-term, it it would be difficult to attribute the long-term changes in richness to a change in just colonization or extinction rates.
#'   
#' **Manuscript paragraph:**  
#' A time series of richness can be decomposed into the colonizations and extinctions of individual species over time. We categorized species according to the following colonization extinction patterns: present in all years = neither (536 species), colonized and went extinct = both (263 species),  initially absent but present every year after its colonization = colonizer (61 species), initially present but absent every year after its extinction = leaver (4 species). Most regions had the same overall ranking (neither > both > colonizer > leaver), except in the Northeast US where both was the most common and neither second, and in the Aleutian Islands where colonizer was the second most common and both third (Figure S4). In general, changes in richness were not due to permanent departures or introductions of species to the region. Furthermore, colonization and extinction rates did not become more dissimilar over time for any region (Figure S5). Colonizations were initially greater than extinctions in Aleutian Islands, Gulf of Alaska, and Gulf of Mexico, but this difference disappeared in the latter portion of these time series, as evidenced by these regions’ initially rapid increase in richness that later plateaued. The other regions did not show strong trends in the difference between colonizations and extinctions over time, making it difficult to attribute the long-term trends in richness to changes in just colonization or just extinction rates.
#'   
#' ###Total Colonizations + Extinctions and Geographic Range
#' The result that richness is predicted by geographic range implied an underlying association between range, colonization/ extinction, and richness itself. This paper begins by explaining richness with range. It will end by explaining how range changes near a colonization/ extinction event. Here, between the two, I'll show how the number of colonizations is related to range.  
#' ####Figure S5
#+ ceEvents-vs-rangeSizeDensity, fig.width=3.5, fig.height=6, fig.cap="**Figure S5.** Number of colonizations and extinctions as a function of range size and range density."
ceEventRange()
#' Yup, this is definitely a thing. Long-term average range size and range density predict how many colonizations and extinctions a species is likely to have. This will lead nicely into examing how range changes prior to an extinction or after a colonization.
#'   
#' \FloatBarrier   
#' 
#' ***  
#' 
#' ##Spatial Clustering of Colonization and Extinction
#+ col-ext-maps-cap, echo=FALSE
col_map_cap <- "**Figure 2.** Maps of long-term averages of colonizations per site per decade for each region: A) E. Bering Sea, B) Gulf of Alaska, C) Aleutian Islands, D) Scotian Shelf, E) West Coast US, F) Newfoundland, G) Gulf of Mexico, H) Northeast US, I) Southeast US. Values of colonization rate were smoothed using a Gaussian kernel smoother. The smoothed colonization rate is indicated by the color bars in each panel; colors are scaled independently for each region."
ext_map_cap <- "**Figure S7.** Maps of long-term averages of extinctions per site per decade for each region: A) E. Bering Sea, B) Gulf of Alaska, C) Aleutian Islands, D) Scotian Shelf, E) West Coast US, F) Newfoundland, G) Gulf of Mexico, H) Northeast US, I) Southeast US. Values of extinction rate were smoothed using a Gaussian kernel smoother. The smoothed extinction rate is indicated by the color bars in each panel; colors are scaled independently for each region."
col_nb_cap <- "**Figure S8.** Connectivity and local spatial autocorrelation of colonization events in each region. Each site is represented by a point. Points connected by a line are neighbors. For each region, neighbors were determined by first calculating the minimum distance required to allow each site to have at least 1 neighbor. Neighbors of a focal point were then defined as the points within this minimum distance from the focal point. Local spatial autocorrelation is local Moran’s I, significant LMI is indicated by a solid point, the color of which indicates the value of the LMI statistic. The outline is the region boundary used for smoothing in Figure 3 (main text), but does not affect calculations of LMI."
ext_nb_cap <- "**Figure S9.** Connectivity and local spatial autocorrelation of extinction events in each region. Each site is represented by a point. Points connected by a line are neighbors. For each region, neighbors were determined by first calculating the minimum distance required to allow each site to have at least 1 neighbor. Neighbors of a focal point were then defined as the points within this minimum distance from the focal point. Local spatial autocorrelation is local Moran’s I, significant LMI is indicated by a solid point, the color of which indicates the value of the LMI statistic. The outline is the region boundary used for smoothing in Figure 3 (main text), but does not affect calculations of LMI."
#'   
#' ####Figure 3. Colonization map
#+ col-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap=col_map_cap
ceRate_map(ce="colonization")
#'   
#' ####Figure S7. Extinction map
#+ ext-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap=ext_map_cap
ceRate_map(ce="extinction")
#' Hotspots can be seen in most regions. Newfoundland also has high values around its edge (as opposed to interior), it seems. NEUS and Gmex show very strong hotspots, and other locations tend to be much much lower. Other regions show more of a continuum.  
#'     
#+ col-ext-intensities, echo=TRUE, fig.width=7, fig.height=3, fig.cap=ext_map_cap, cache=TRUE
rel_col_ext_rate <- mapDat[,j={
	map_smooth_col <- spatstat::Smooth(spatstat::ppp(x=.SD[,lon], y=.SD[,lat], marks=.SD[,n_spp_col_weighted], window=mapOwin[[reg]]), hmax=1)
	mark_range_col <- range(map_smooth_col, na.rm=TRUE)*10
	
	map_smooth_ext <- spatstat::Smooth(spatstat::ppp(x=.SD[,lon], y=.SD[,lat], marks=.SD[,n_spp_ext_weighted], window=mapOwin[[reg]]), hmax=1)
	mark_range_ext <- range(map_smooth_ext, na.rm=TRUE)*10
	ol <- list(
		minval_col=mark_range_col[1], maxval_col=mark_range_col[2], max_o_min_col=do.call("/",as.list(rev(mark_range_col))),
		minval_ext=mark_range_ext[1], maxval_ext=mark_range_ext[2], max_o_min_ext=do.call("/",as.list(rev(mark_range_ext)))
	)
	lapply(ol, function(x)if(is.numeric(x)){signif(x,3)}else{x})
},by=c("reg")]
kable(rel_col_ext_rate, caption="Table. The colonization and extinction intensity range and max/min ratio")
kable(rel_col_ext_rate[,lapply(.SD, median)], caption="Table. Median of above table")


#' ####Figure S8. Colonization neighborhood
#+ col-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap=col_nb_cap
nb_moranI(ce="colonization")
#'   
#' ####Figure S9. Extinction neighborhood
#+ ext-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap=ext_nb_cap
nb_moranI(ce="extinction")
#'   
#' \FloatBarrier   
#' 
#' ***  
#' 
#' ##Geographic Range Near Colonization and Extinction
#'   
#' ####Figure 4
#+ rangeColExt, fig.width=6, fig.height=6, fig.cap="Figure 4. Geographic range size (A,B) and geographic range density (C,D) vs years until extinction (A,C) and years after colonization (B,D). For each unique value on the horizontal axis, the cross-species average for the range metric is displayed, and a linear model fit through this average. Statistics in main text do not use this aggregation. Extinction events are identified as occuring the year before the species is absent (?right?), colonization the first year it is present after an absence."
rangeSize_absenceTime()
#' Range size declines near an absence much more consistently than does range density; both are (relatively) low just before extinction and just after colonization. However, range density has much more variable intercepts among regions, whereas range size does not.   
#'   
#' This makes sense, at least somewhat, because colonization and extinction events are defined at the site level; though the outcome isn't necessitated by this formulation, because size could drop suddenly. In fact, when a species is absent, both its range size and its range density must be 0 (though, range density is technically calculated for only those sites that are occupied, so I supposed it's technically undefined according to the equations I'm using).  
#'   
#' I think the regressions for range size should omit an intercept, while the regressions for range density should have it. This might be hard to justify fully *a priori* (though see my thinking in previous paragraph), so I'll probably just do a model selection and maybe discuess the difference if one model has an intercept and the other does not.  
#'   
#' ###Changes in Range Size Near Colonization or Extinction
#+ rangeSize-ColExt
rangeTimeDT <-  spp_master[!is.na(stretch_type) & propStrata!=0]
rangeTimeDT <- rangeTimeDT[,list(
	reg=reg, 
	event=as.character(event_year), 
	spp=spp, 
	type=as.character(stretch_type),
	time=ext_dist, 
	size=propStrata, 
	density=propTow_occ
)]

sizeColExt_mods <- list()
sizeColExt_mods[[1]] <- lmer(size ~ time*type + (time*type | spp/reg), data=rangeTimeDT)
sizeColExt_mods[[2]] <- lmer(size ~ time*type + (time | spp/reg) + (type|reg), data=rangeTimeDT)
sizeColExt_mods[[3]] <- lmer(size ~ time + type + (time | spp/reg) + (type|reg), data=rangeTimeDT)
sizeColExt_mods[[4]] <- lmer(size ~ time + (time | spp/reg), data=rangeTimeDT) # this is what I settled on previously i think
sizeColExt_mods[[5]] <- lmer(size ~ time + (1 | spp/reg), data=rangeTimeDT)
sizeColExt_mods[[6]] <- lmer(size ~ time + (time - 1 | spp/reg), data=rangeTimeDT)

# do.call(stargazer, c(sizeColExt_mods, list(type="html")))



#' ###Changes in Range Density Near Colonization or Extinction
#+ rangeDens-ColExt
densityColExt_mods <- list()
densityColExt_mods[[1]] <- lmer(density ~ time*type + (time*type | spp/reg), data=rangeTimeDT)
densityColExt_mods[[2]] <- lmer(density ~ time*type + (time | spp/reg) + (type|reg), data=rangeTimeDT)
densityColExt_mods[[3]] <- lmer(density ~ time + type + (time | spp/reg) + (type|reg), data=rangeTimeDT)
densityColExt_mods[[4]] <- lmer(density ~ time + (time | spp/reg), data=rangeTimeDT) 
densityColExt_mods[[5]] <- lmer(density ~ time + (1 | spp/reg), data=rangeTimeDT)
densityColExt_mods[[6]] <- lmer(density ~ time + (time - 1 | spp/reg), data=rangeTimeDT)

do.call(stargazer, c(sizeColExt_mods, densityColExt_mods, list(type="html")))


