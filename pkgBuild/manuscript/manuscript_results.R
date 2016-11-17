#' ---
#' title: "Marine species' range shifts drive long-term changes in richness"
#' author: "Ryan Batt"
#' date: "2016-08-15"
#' abstract: |
#'       Several regions show long-term changes in richness, generally increases (**Figure 1**, **Table 2**). The rest of this study focuses on understanding these long-term changes through the lens of the geographic ranges of individual species and how geographic range relates to the colonization and extinction of those species, and thus to changes in richness.   
#'         
#'       Any change in richness must be tied to an assemblage-level imbalance between the rates of colonization and extinction of individual species. Most species, if not present every year, both colonized and went extinct, possibly repeatedly (**Figure S2**). However, the long-term geography of these two processes shows hotspots (**Figure 2**, **Figure S3**) of changing richness in the form of spatial clustering of sites with high rates of colonization and extinction (**Figure S4**, **Figure S5**). In other words, richness changes take place at a limited number of nearby sites that are characteristically occupied by species that have just gone extinct or have just colonized the region. This result suggests that the geographic distributions of these species might be closely linked to observed changes in richness.  
#'         
#'       We analyzed geographic distribution of a species in two ways: range size (proportion of sites) and range density (proportion of tows in occupied sites). These metrics are year- and species-specific; their long-term average yields a value characteristic of the species, and averaging characteristic values across species present in a year provides a community-scale metric of geographic range. The two metrics are closely associated (**Figure S6**), and both predict predict species richness (**Figure 3**).   
#'         
#'       Why is richness predicted by geographic range? Richness is highest when the typical range of present species is small and/or sparse (**Figure 3**), which means that richness is higher when more rare species are present. Colonization and extinction are essential for richness to change, and basic macroecological theory states that species that are about to go extinct or that have just colonized are relatively rare. Therefore, when more rare species are present, richness will be higher because a large fraction of the assemblage has either just colonized or is about to go extinct, which means there will be more species in addition to the baseline common species (which, because they are common, probably don't go extinct or colonize often in the time series, and so are typically present, ergo "baseline").      
#'         
#'       The explanation for why richness is related to geographic range hinges on the notion that as a species becomes more rare it is closer to extinction, and that after a species colonizes it becomes more common. If this pattern holds, and if the intercept and slope doesn't change much among species, then it should be true that species that are more rare on the long-term should colonize and go extinct more often. Indeed, we find that species with higher total colonization and extinction numbers do tend to be more rare over the long-term (**Figure S7**). Furthermore, we find that species range size and density reliably decline as extinction approaches, and increase after colonization (**Figure 4**).   
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
options("digits"=3) # rounding output to 4 in kable() (non-regression tables)
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
# source("../manuscript/fig_tbl_number.R")
# eval(fig_tbl_number())


# ============
# = Richness =
# ============
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Results
#' ##Species Richness
#' ###Richness Time Series
#' ####Table. Summary of richness mean and variability
#+ richnessSummary-basic-table, echo=FALSE
kable(
	data.table(
		comm_master[,list(mu=mean(reg_rich), sd=stats::sd(reg_rich)),by='reg'], 
		comm_master[, list(mu_sd=mean(.SD[,stats::sd(reg_rich),by='reg'][,V1]))]
	), 
	caption="Long-term mean and standard deviation in MSOM richness, and average of those sd's. Gmex and Shelf have the highest and lowest long-term averages in species richness (MSOM)."
)

#' ####Figure 1. MSOM richness time series
#+ Richness-ts-fig, fig.height=5, fig.width=3.5, ecal=TRUE, echo=TRUE, fig.cap="**Figure 1.** Time series of MSOM estimates of region richness. Each point is the posterior mean of regional richness in a year. Lines indicate long-term trends from fitted values of linear regression models predicting richness from time."
richness_ts()
#' ####Figure S1. MSOM - naive scatter
#+ Richness-msom-naive-scatter-fig, fig.height=3.5, fig.width=3.5, fig.cap="**Figure S1.** MSOM richness vs naive richness"
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
#+ Richness-trend
#' ###Richness Trend
#' Looking for trends in species richness in both the naive estimates and the MSOM estimates. Using Kendall's Tau_b, calculated for 1E4 resamplings of the posterior of richness. Tau is also calculated using a method that removes serial correlation to achieve independence of observations so that the p-values are correct.
load("../../pkgBuild/results/rich_naive_trend_kendall.RData")
rich_naive_trend_kendall[reg!="wcann",BH:=p.adjust(taup, method='BH')]
rich_naive_trend_kendall <- rich_naive_trend_kendall[
	reg!="wcann",list(reg=reg, estimate=tau, BH=BH, p.value=taup)
]
load("../../pkgBuild/results/rich_trend_kendall.RData")
rich_trend_kendall[reg!="wcann",BH:=p.adjust(pvalue, method="BH")]
rich_trend_kendall <- rich_trend_kendall[
	reg!="wcann",list(reg=reg, estimate=tau, BH=BH, p.value=pvalue)
]
#' ####Table S1. Naive tau
#+ naiveTau-table, echo=FALSE
kable(rich_naive_trend_kendall, 
	caption="**Table S1.** Naive richness trends in each region. Estimate is Kendall's Tau_b, BH is the Benjamini-Hochberg corrected p-value, and p.value is the original p-value."
)
#' In the Naive estimates, `r rich_naive_trend_kendall[reg!='wcann', sum(p.value<=0.05)]` regions had significant $\tau_b$.  
#'   
#' ####Table 2. MSOM tau
#+ msomTau-table, echo=FALSE
kable(rich_trend_kendall, 
	caption="**Table 2.** MSOM richness trends in each region. Estimate is Kendall's Tau_b, BH is the Benjamini-Hochberg corrected p-value, and p.value is the original p-value. Trends are calculated for 1E4 resampled combinations of the richness posterior. For each resampling, independence of observations is achieved before estimating Tau by removing serial correlation from resampled time series. This procedure retains the integrity of p-values."
)
#' In the MSOM estimates, `r rich_trend_kendall[reg!='wcann', sum(p.value<=0.05)]` regions had significant $\tau_{b}$.  
#'   
#' Overall, ~half the regions show positive trends. No regions show significant negative trends (although SEUS is close, depending on the analysis). It's a bit surprising how much removing the autocorrelation seemed to impact some of the trends (AI in particular, I think).  
#'   
#' **Manuscript paragraph:**  
#' Estimated slopes for long-term trends (Kendall's $\tau_b$) in richness were positive for most regions. For Naive estimates, all the seven positive $\tau$ were also significantly different from 0, whereas the two negative $\tau_b$, Aleutian Islands and Southeast US, were not (**Table S2**). Long-term trends in MSOM estimates of richness had the same sign as Naive trends, except for Aleutian Islands, but now the only significant $\tau_b$ were the following four regions: Eastern Bering Sea ($\tau_b$ = 0.41), West Coast US ($\tau_b$ = 0.61), Scotian Shelf ($\tau_b$ = 0.45), and Newfoundland ($\tau_b$ = 0.72) (**Table 1**). Richness trends were not significant in the Gulf of Mexico, Gulf of Alaska, and Northeast US, Aleutian Islands, Southeast US.   
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Colonization and Extinction  
#' ###Colonization and Extinction Summary
#' ####Figure S2.  Barplot of categories
#+ Richness-col-ext-barplot, include=TRUE, echo=TRUE, fig.height=3.5, fig.width=3.5, fig.cap="**Figure S2.** Number of species beloning to the categories of both, neither, colonizer, leaver in each region"
categ_barplot()
#+ Richness-col-ext-barplot-table, echo=FALSE
categ_tbl <- t(spp_master[!duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ)])[c(4,1,2,3),]
kable(categ_tbl, caption = "Number of species in each category in each region.")
#' It's the same pattern, whichever way you split it. However, AI is the only region that had more *colonizers* than *both* species. An interesting way to think about some of this is that the average sd in richness was `r comm_master[,stats::sd(reg_rich),by='reg'][,mean(V1)]`, so when the number of *colonizer* or *leaver* species exceed's that region's sd, the impact of those categories, which I consider to be dubious, might start being relevant (though it's not necessarily problematic, nor is this even close to an actual test for the significance of those categories to the trend). EBS and Shelf had significant positive trends in richness and very low numbers in the *colonizer* category. WCTRI and NEWF had similar numbers in the *both* and *colonizer* category.  
#'   
#' ####Figure NotIncluded. Time series of colonizations and extinctions
#+ Richness-col-ext-ts, include=TRUE, echo=TRUE, fig.height=3.5, fig.width=3.5, fig.cap="**Figure NotIncluded.** Number of colonizations (blue) and extinctions (red) over time in each region."
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
#' ##Spatial Clustering of Colonization and Extinction
#' ###Heat Maps
#' ####Figure 2. Colonization map
#+ col-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure 2.** Maps of long-term averages of colonizations per site per decade for each region: A) E. Bering Sea, B) Gulf of Alaska, C) Aleutian Islands, D) Scotian Shelf, E) West Coast US, F) Newfoundland, G) Gulf of Mexico, H) Northeast US, I) Southeast US. Values of colonization rate were smoothed using a Gaussian kernel smoother. The smoothed colonization rate is indicated by the color bars in each panel; colors are scaled independently for each region."
ceRate_map(ce="colonization")
#'   
#' ####Figure S3. Extinction map
#+ ext-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure S3.** Maps of long-term averages of extinctions per site per decade for each region: A) E. Bering Sea, B) Gulf of Alaska, C) Aleutian Islands, D) Scotian Shelf, E) West Coast US, F) Newfoundland, G) Gulf of Mexico, H) Northeast US, I) Southeast US. Values of extinction rate were smoothed using a Gaussian kernel smoother. The smoothed extinction rate is indicated by the color bars in each panel; colors are scaled independently for each region."
ceRate_map(ce="extinction")
#' Hotspots can be seen in most regions. Newfoundland also has high values around its edge (as opposed to interior), it seems. NEUS and Gmex show very strong hotspots, and other locations tend to be much much lower. Other regions show more of a continuum.  
#'     
#+ col-ext-intensities, echo=TRUE,  cache=FALSE
rel_col_ext_rate <- mapDat[,j={
	map_smooth_col <- spatstat::Smooth(spatstat::ppp(
		x=.SD[,lon], y=.SD[,lat], 
		marks=.SD[,n_spp_col_weighted], window=mapOwin[[reg]]
	), hmax=1)
	mark_range_col <- range(map_smooth_col, na.rm=TRUE)*10
	
	map_smooth_ext <- spatstat::Smooth(spatstat::ppp(
		x=.SD[,lon], y=.SD[,lat], 
		marks=.SD[,n_spp_ext_weighted], window=mapOwin[[reg]]
	), hmax=1)
	mark_range_ext <- range(map_smooth_ext, na.rm=TRUE)*10
	ol <- list(
		minval_col=mark_range_col[1], maxval_col=mark_range_col[2], 
		max_o_min_col=do.call("/",as.list(rev(mark_range_col))),
		minval_ext=mark_range_ext[1], maxval_ext=mark_range_ext[2], 
		max_o_min_ext=do.call("/",as.list(rev(mark_range_ext)))
	)
	lapply(ol, function(x)if(is.numeric(x)){signif(x,3)}else{x})
},by=c("reg")]
#+ col-ext-intensities-table, echo=FALSE
kable(
	rbind(rel_col_ext_rate, rel_col_ext_rate[,lapply(.SD, median)][,reg:="MEDIAN"]), 
	caption="The colonization and extinction intensity range and max/min ratio, and median among regions. Useful for assessing how big of a difference there is between red and blue for each region."
)
#'   
#' ###Neighborhoods and Local Moran's I
#' ####Figure S4. Colonization neighborhood
#+ col-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure S4.** Connectivity and local spatial autocorrelation of colonization events in each region. Each site is represented by a point. Points connected by a line are neighbors. For each region, neighbors were determined by first calculating the minimum distance required to allow each site to have at least 1 neighbor. Neighbors of a focal point were then defined as the points within this minimum distance from the focal point. Local spatial autocorrelation is local Moran’s I, significant LMI is indicated by a solid point, the color of which indicates the value of the LMI statistic. The outline is the region boundary used for smoothing in Figure 3 (main text), but does not affect calculations of LMI."
nb_moranI(ce="colonization")
#'   
#' ####Figure S5. Extinction neighborhood
#+ ext-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap="**Figure S5.** Connectivity and local spatial autocorrelation of extinction events in each region. Each site is represented by a point. Points connected by a line are neighbors. For each region, neighbors were determined by first calculating the minimum distance required to allow each site to have at least 1 neighbor. Neighbors of a focal point were then defined as the points within this minimum distance from the focal point. Local spatial autocorrelation is local Moran’s I, significant LMI is indicated by a solid point, the color of which indicates the value of the LMI statistic. The outline is the region boundary used for smoothing in Figure 3 (main text), but does not affect calculations of LMI."
nb_moranI(ce="extinction")
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Richness and Geographic Range
#' ###Richness and Range Density, Range Size
#' Looking for spatial footprint that will predict richness. Here I characterize the geographic distribution of a species as its range size. This is tricky for richness because richness is a community level trait, but I've characterized range size at the species level. In fact, this "species level" attribute is also potentially dyanmic -- a species need not have fixed range size. So there are two levels of averaging -- for each species take its long-term average value, then for each community take the average across species.  
#' 
#' Note that I'd previously used range density in addition to range size. I've dropped the density metric in most places in order to simplify the results, because the range size results were closely related to the range density results.  
#'  
#' First posterity, I'll show a figure relating range size and range density. Then I'll present the figure of how size can 'predict' richness. Finally, I'll explore models more formally expressing the relationship between richness and size. 
#'  
#' ####Figure S6. Relationship between range size and range density
#+ geo-range-densVsize, fig.width=3.5, fig.height=5.5, fig.cap="**Figure S6**. Range density versus range size. In panel A each point is a species-region-year combination. In panel B, each point is a region-year. Range size is the proportion of sites occupied, range density the tows in occupied sites. The community metrics in B is calculated by take each species' long-term average from A, then taking the average across all species present in the community in a given year. Fitted lines in A are from a loess fit."
rangeSizeDens()
#' These plots show that there is a strong relationship between range size and range density. Interestingly, in **Figure S6B**, the cross-region relationship is negative (if each color had 1 pt), whereas the within-region relationship is positive.  
#'   
#' The positive relationship between size and density is not surprising. My interpretation of density is ~population size. Often population size is correlated with range size. I think this is a standard result, but I need to double-check.
#'   
#' ####Figure 3. Species richness versus geographic range size
#+ rich-geo-rangeSize, fig.width=3.5, fig.height=3.5, fig.cap="**Figure 3.** Species richness vs geographic range size. Range size is presented as each species' long-term average of the proportion of sites it occupied. Solid lines are linear regressions with MSOM richness as the response and the horizontal axis and an intercept as the predictors."
rich_geoRange("size", leg=TRUE, legPan=1, panLab=FALSE)

#' Range size is a pretty good predictor of species richness. I think I had originally missed the range size relationship b/c I hadn't done the same aggregating procedure. The interpretation I have is that richness is highest when you have a bunch of rare species.  
#'   
#' The goal here is to see if species richness is predicted by the typical range size of community's constituent species. First I'll run different types of models just to explore whether this is true, in general (across regions). Then I'll drill in to each region individually to answer the same question.  
#'   
#' ####Table. Regressions relating richness to range size
#+ rich-rangeSize
# This is a function that'll help summarize model fit and coeffs and parameter significance:  
mod_smry <- function(m, pred_name=c("density","size","time","type","time:type")){
	pred_name <- match.arg(pred_name)
	sc <- sem.coefs(m)
	sc[,c("estimate","p.value")] <- lapply(sc[,c("estimate","p.value")], signif, 4)
	mod_call <- switch(class(m), lmerMod=m@call, lm=m$call)
	mod_call <- as.character(mod_call)[2]
	fits <- sem.model.fits(m)
	fits[,c("estimate","p.value","Marginal","Conditional")] <- lapply(fits[,c("Marginal","Conditional")], signif, 4)
	out <- cbind(
		mod_call = mod_call,
		sc[sc[,"predictor"]==pred_name,],
		fits, 
		AIC=round(AIC(m), getOption("digits"))
	)
	out[,c(
		# "std.error","N",
		"Class","mod_call","predictor","estimate","p.value","Marginal","Conditional","AIC"
	)]
}

# Make a data set that is useful for these regressions (short variable names, etc):  
range_reg <- comm_master[,list(
	reg, year, rich=reg_rich, density=propTow_occ_avg, size=propStrata_avg_ltAvg
)]

# Fit different models to the whole data set
rSize_mods <- list()
rSize_mods[[1]] <- lm(rich ~ size, data=range_reg)
rSize_mods[[2]] <- lm(rich ~ size*reg, data=range_reg)
rSize_mods[[3]] <- lme4::lmer(rich ~ size + (1|reg), data=range_reg)
rSize_mods[[4]] <- lme4::lmer(rich ~ size + (size|reg), data=range_reg)
rich_size_smry <- rbindlist(lapply(rSize_mods, mod_smry, pred_name="size"))

# Fit same model to each region separately 
rSize_reg_mods <- list()
ur <- range_reg[,unique(reg)]
for(r in 1:length(ur)){
	rSize_reg_mods[[r]] <- lm(rich ~ size, data=range_reg[reg==ur[r]])
}
rich_s_reg_smry <- data.table(
	reg=ur, 
	rbind(
		rbindlist(lapply(rSize_reg_mods, mod_smry, pred_name="size"))
	)
)
setkey(rich_s_reg_smry, reg, Class, predictor)

setnames(rich_s_reg_smry, old=c("Marginal","Conditional","p.value"), new=c("MargR2","CondR2","pval"))
setkey(rich_s_reg_smry, reg, Class, mod_call, predictor)
rich_s_reg_smry[, BH:=round(p.adjust(pval, "BH"), 3), by=c("mod_call", "predictor")]
rich_s_reg_smry <- dcast(rich_s_reg_smry, reg+Class+mod_call+MargR2+CondR2+AIC~predictor, value.var=c("pval", "BH"))

#+ rich-rangeSize-tables, echo=FALSE
kable(
	rbind(rich_size_smry), 
	caption="Summary of regression models involving richness and community metric of range size."
) # different kinds of models
kable(
	rich_s_reg_smry, 
	caption="Regression models of richness predicted by community metric of range size, but each region has a separate model."
) # same model applied to each reg sep
kable(
	rich_s_reg_smry[,lapply(.SD, base::mean), by="mod_call"], 
	caption="Average of above region-specific rich ~ size models."
) 
#' All models are pretty good predictors. Well, the most basic model kinda sucks I guess. It needs to account for some of the between-region variation.  
#'   
#' ###Predicting Richness: Range Size or Density?
#' As far as picking one or the other, it doesn't end up mattering much. Range size is a lot better than density in NEUS, and density outperforms size in AI. Otherwise, size as a slight edge over density on average, although both predictors are significant in all regions.  
#'   
#' ###Conclusion for Richness and Geographic Range
#' The two metrics of geographic range are well correlated. Furthermore, richness can be predicted pretty well using regressions with either as a predictor. There are large differences among regions, though. This is probably because richness is not readily comparable among most regions. Regions vary mostly in their intercept values, and they have fairly similar slopes (though they are not identical, and model fits improve when allowing slopes to vary among regions; it's just that the improvement is small compared to allowing intercepts to vary among regions).  
#'  
#' The interpretation of the result that geographic distribution predicts species richness is likely associated with species rarity. When the average range density or range size of a community is low, it means it has a lot of species that are rare (at either spatial scale). It's these rare species that come and go, and form the dynamics of richness that we observe. When that dynamical value is high, it implies that an above-average number of the dynamic species are present. Because those species are transient (dynamic), they are also probably rare.  
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Geographic Range and Colonizations and Extinctions
#'   
#' ###Total Colonizations vs Extinctions and Geographic Range
#' The result that richness is predicted by geographic range implied an underlying association between range, colonization/ extinction, and richness itself. Above, richness was explained with range. Later, the results will explain how range changes near a colonization/ extinction event. Here, between the two, the results show how the number of colonizations is related to range.  
#'   
#' ####Figure S7. Total Colonizations/ Extinctions vs Georaphic Range
#+ ceEvents-vs-rangeSizeDensity, fig.width=3.5, fig.height=6, fig.cap="**Figure S7.** Number of colonizations and extinctions as a function of range size and range density."
ceEventRange()
#' Yup, this is definitely a thing. Long-term average range size and range density predict how many colonizations and extinctions a species is likely to have. This will lead nicely into examing how range changes prior to an extinction or after a colonization.  
#'   
#' ###Range Change after Colonization/ before Extinction
#' ####Figure 4. Changes in Geographic Range before Extinction and after Colonization
#+ rangeSize_ColExt, fig.width=6, fig.height=3, fig.cap="**Figure 4.** Geographic range size vs years until extinction (A) and years after colonization (B). For visualization purposes, range size is averaged across species for each unique value on each axis, and a linear model fit through this average. Statistics in main text use unaggregated data. The horizontal axes were formulated as time until (since) the nearest upcoming (previous) absence. Because range size must be zero when either horizontal axis has a value of zero, points at (0,0) were excluded from figures and analyses."
rangeSize_absenceTime("rangeSize")
#' ####Figure 4b. Changes in Geographic Range before Extinction and after Colonization
#+ rangeDensity_ColExt, fig.width=6, fig.height=3, fig.cap="**Figure 4b.** Geographic range density vs years until extinction (A) and years after colonization (B). For visualization purposes, range density is averaged across species for each unique value on each axis, and a linear model fit through this average. Statistics in main text use unaggregated data. The horizontal axes were formulated as time until (since) the nearest upcoming (previous) absence. Because range density must be zero when either horizontal axis has a value of zero, points at (0,0) were excluded from figures and analyses."
rangeSize_absenceTime("rangeDensity")
#' Range size declines near an absence much more consistently than does range density; both are (relatively) low just before extinction and just after colonization. However, range density has much more variable intercepts among regions, whereas range size does not.   
#'   
#' This makes sense, at least somewhat, because colonization and extinction events are defined at the site level; though the outcome isn't necessitated by this formulation, because size could drop suddenly. In fact, when a species is absent, both its range size and its range density must be 0 (though, range density is technically calculated for only those sites that are occupied, so I supposed it's technically undefined according to the equations I'm using).  
#'   
#' I think the regressions for range size should omit an intercept, while the regressions for range density should have it. This might be hard to justify fully *a priori* (though see my thinking in previous paragraph), so I'll probably just do a model selection and maybe discuess the difference if one model has an intercept and the other does not.  
#'   
#' ####Table. Regressions w/ regions pooled relating time until extinction or after colonization to geographic range
#+ rangeSize-ColExtTime-data
# handy data set for regressions
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

# summary function used in tables
mod_smry2 <- function(m){
	mod_call <- switch(class(m), lmerMod=m@call, lm=m$call)
	mod_call <- as.character(mod_call)[2]
	fits <- sem.model.fits(m)
	fits[,c("estimate","p.value","Marginal","Conditional")] <- lapply(fits[,c("Marginal","Conditional")], signif, 4)
	anova_sum <- car::Anova(m)
	
	out <- data.frame(
		mod_call = mod_call,
		data.frame(predictor=rownames(anova_sum), anova_sum)[,c("predictor","Pr..Chisq.")],
		fits, 
		AIC=AIC(m)#signif(AIC(m), getOption("digits"))
	)
	out[,c("Class","mod_call", "predictor", "Pr..Chisq.","Marginal","Conditional","AIC")]
}

#+ rangeSize-ColExtTime-models, results='markup'
# models for range size
sizeCE_mods <- list()
sizeCE_mods[[1]] <- lme4::lmer(size ~ time + (time|spp/reg), data=rangeTimeDT)
# sizeCE_mods[[2]] <- lme4::lmer(size ~ time + (time|spp/reg) + (1|type/reg), data=rangeTimeDT) # error
sizeCE_mods[[2]] <- lme4::lmer(size ~ time*type + (time|spp/reg), data=rangeTimeDT)
# lmerTest::step(sizeCE_mods[[2]]) # doing most complicated throws an error

#+ rangeSize-ColExtTime-table, echo=FALSE
do.call(stargazer, c(
		sizeCE_mods, 
		report = c("vc*"), header=FALSE,
		title = "Mixed effect models predicting range size (size) from the temporal proximity (time) after colonization and/or until extinction (factor name = type). Each model uses data from all regions (reg) and all species (spp), and includes both time to extinction and time from colonization (type).",
		list(type=table_type)
))

tbl_ColExtTime <- rbindlist(lapply(c(sizeCE_mods), mod_smry2))
setnames(tbl_ColExtTime, old=c("Pr..Chisq.","Marginal","Conditional"), new=c("pval","MargR2","CondR2"))
kable(tbl_ColExtTime, caption="Same as above, but shows slightly different metrics")

#' Range size seems to be important to include type as a fixed effect. It does not need an interaction between time*type. I did additional testing beyond what's presenting here, and I can confirm that having spp as a random factor is useful, too.  
#'   
#' ####Table. Regressions separate regions -- geographic range vs time until extinction or after colonization
#+ rangeSizeDensity-ColExtTime-reg
# Fit same model to each region separately 
sTime_reg_mods <- list()
ur <- range_reg[,unique(reg)]
for(r in 1:length(ur)){
	sTime_reg_mods[[r]] <- lme4::lmer(size ~ time + type + (time|spp), data=rangeTimeDT[reg==ur[r]])
}
sTime_reg_smry <- data.table(rbind(
	rbindlist(structure(lapply(sTime_reg_mods, mod_smry2), .Names=ur), idcol=TRUE)
))
setnames(sTime_reg_smry, old=c(".id","Pr..Chisq.","Marginal","Conditional"), new=c('reg',"pval","MargR2","CondR2"))
setkey(sTime_reg_smry, reg, Class, mod_call, predictor)

# adjust p-values for multiple tests
sTime_reg_smry[, BH:=round(p.adjust(pval, "BH"), 3), by=c("mod_call", "predictor")]

# rearrange so each model on 1 line
sTime_reg_smry <- dcast(sTime_reg_smry, reg+Class+mod_call+MargR2+CondR2+AIC~predictor, value.var=c("pval", "BH"))

# add coefficients to the sTime_reg_smry data.table
getCoefs <- function(mList){
	outList <- list()
	for(r in 1:length(ur)){
		outList[[ur[r]]] <- rbindlist(lapply(
			coef(mList[[r]]), function(x)as.list(colMeans(x))
		), idcol=TRUE)
		setnames(outList[[ur[r]]], c(".id"), c("randomGroup"))
		mod_call <- switch(class(mList[[r]]), lmerMod=mList[[r]]@call, lm=mList[[r]]$call)
		mod_call <- as.character(mod_call)[2]
		outList[[ur[r]]][,mod_call:=mod_call]
	}
	outList <- rbindlist(outList, idcol=TRUE)
	setnames(outList, c(".id"), c("reg"))
	# setcolorder(outList, c("reg", "mod_call", "randomGroup", "(Intercept)", "time", "typepre_ext"))
	outList
}
modCoef <- getCoefs(sTime_reg_mods)
sTime_reg_smry <- merge(sTime_reg_smry, modCoef, by=c("reg","mod_call"))
setnames(sTime_reg_smry, old=c("(Intercept)"), new=c("Int"))

#+ rangeSizeDensity-ColExtTime-reg-table, echo=FALSE
capS <- sTime_reg_smry[grepl("size",mod_call),mod_call[1]]
sTime_reg_smry <- sTime_reg_smry[grepl("size",mod_call)]
sTime_reg_smry[,c("mod_call","randomGroup","Class"):=NULL]
kable(sTime_reg_smry, caption=paste0("Summary statistics for fits of predicting range SIZE from years before extinction and years after colonization. These are mixed effect models of the form ", capS))

#' ####Table. Regressions separate regions -- range size vs time until, WITH INTERACTION
#+ rangeSizeDensity-ColExtTime-reg-interaction
sTime_reg_mods2 <- list()
ur <- range_reg[,unique(reg)]
for(r in 1:length(ur)){
	sTime_reg_mods2[[r]] <- lme4::lmer(size ~ time * type + (time|spp), data=rangeTimeDT[reg==ur[r]])
}
sTime_reg_smry2 <- data.table(rbind(
	rbindlist(structure(lapply(sTime_reg_mods2, mod_smry2), .Names=ur), idcol=TRUE)
))
setnames(sTime_reg_smry2, old=c(".id","Pr..Chisq.","Marginal","Conditional"), new=c('reg',"pval","MargR2","CondR2"))
setkey(sTime_reg_smry2, reg, Class, mod_call, predictor)

# adjust p-values for multiple tests
sTime_reg_smry2[, BH:=round(p.adjust(pval, "BH"), 3), by=c("mod_call", "predictor")]

# rearrange so each model on 1 line
sTime_reg_smry2 <- dcast(sTime_reg_smry2, reg+Class+mod_call+MargR2+CondR2+AIC~predictor, value.var=c("pval", "BH"))

modCoef2 <- getCoefs(sTime_reg_mods2)
sTime_reg_smry2 <- merge(sTime_reg_smry2, modCoef2, by=c("reg","mod_call"))
setnames(sTime_reg_smry2, old=c("(Intercept)"), new=c("Int"))

#+ rangeSizeDensity-ColExtTime-reg-interaction-table, echo=FALSE
capS <- sTime_reg_smry2[grepl("size",mod_call),mod_call[1]]
sTime_reg_smry2 <- sTime_reg_smry2[grepl("size",mod_call)]
sTime_reg_smry2[,c("mod_call","randomGroup","Class"):=NULL]
kable(sTime_reg_smry2, caption=paste0("Summary statistics for fits of predicting range SIZE from years before extinction and years after colonization. These are mixed effect models of the form ", capS))


#' ###Conclusion for Range Change before Extinction/ after Colonization
#' As stated with the larger (pooled regional data) regressions, range size respondes more consistently/ strongly to an approaching extinction/ departure from colonization than does range density. Range size shrinks as extinction approaches, increases as time since colonization increases. In both cases, there was only 1 region for which the direction (into extinction, out of colonization mattered): for range density, it was the Southeast US (effect of pre-extinction direction = 0.037, p=0.043 [BH correction], Table 12), and for range size it was Scotial Shelf (effect = -0.021, p = 0.002). Time was a significant predictor of range size in all regions; time was *not* a significant predictor of range density in Newfoundland (p = 0.846 [BH]), and the Scotial Shelf (p = 0.916).  
#'  
#' These results support the hypothesis that species rarity is closely associated with proximity to extinction/ colonization, and therefore the probability that the species will contribute to a change in species richness. Furthermore, these relationships have not been directly quantified for entire assemblages. Competing hypotheses exist regarding how range should change approaching extinction and after colonization. We find that the change in range is similar regardless of direction. However, spatial scale was important. Although slopes were similar, the intercept for density was far larger than for range size --- the fraction of tows containing a species in occupied sites may remain high even when extinction is imminent, and may similar be large even if colonization was recent. Range size, on the other hand, was much closer to 0 near an absence.  
#'   
#' Overall, the results suggest that the spatial footprint of individual species is important for understanding changes in species richness. Furthermore, because species contributing most to the dynamics of richness are those that repeatedly colonize and go extinct, it is meaningful to look at a species' long-term rarity in order to gauge whether it is likely to contribute to those long-term richness changes. Determining what drives the geographic range of individual species is probably a powerful way to anticipate richness changes.  
#'   
#' ####Exploring Range size vs absence time as individual regressions
#+ rangeSize_time_sepRegs, fig.width=5, fig.height=6, fig.cap="**Exploration Figure** histograms of separate regressions of size ~ time; this is for each run-up to an extinction and reach follow-up to a colonization. Trying to understand how regularly the pattern might be observed. Hard to answer because adjusting the restriction on number of events in the run-up or follow-up (nTime) greatly affects the proportion that are significant."
rangeTimeDT[,nTime:=length(time),by=c("reg","type","event","spp")]
getEPR <- function(x){
	col_select <- c("Estimate","Pr(>|t|)")
	sx <- summary(x)
	EP <- sx$coefficients[2,col_select]
	names(EP) <- c("Estimate","Pr")
	R <- sx$r.squared
	data.table(as.data.table(as.list(EP)), Rsquared=R)
}
o <- rangeTimeDT[nTime>=3,j={
	getEPR(lm(size~time))
	} ,by=c("reg","type","event","nTime","spp")
]

par(mfcol=c(3,2), mar=c(2,2,0.5,0.5), cex=1, ps=10, mgp=c(1,0.1,0), tcl=-0.1)
hist(o[,Estimate])
hist(o[,(Pr)])
hist(o[,(Rsquared)])
setorder(o, nTime)
o[,j={plot(nTime,Pr); lines(spline(Pr~nTime),col='red',lwd=2)}]
o[,j={plot(nTime,Estimate); ; lines(spline(Estimate~nTime),col='red',lwd=2)}]


o[,list(propSignificant=(sum(Pr<0.05)/sum(!is.na(Pr))),n=sum(!is.na(Pr))),by=c("reg","type")]
o[,list(propSignificant=(sum(Pr<0.05)/sum(!is.na(Pr))),n=sum(!is.na(Pr))),by=c("reg","type")][,mean(propSignificant),by=c("type")]

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
#' #Supplement
#' ##Bin Size
#+ bin_size, fig.width=7, fig.height=7, fig.cap="**Bin Size** How Bin size (lon, lat, depth) affects the number of sites in a region. The vertical axis is the size of the depth bin in meters, horizontal axis is the size of the lon-lat bin in degrees, and the colors indicate the number of sites in that region that are sampled for at least 85% of years."
load("../manuscript/manuscript_binSize.RData")
par(mfrow=c(3,3), mar=c(3,3,0.5,3), oma=c(1,1,1,0.5))
for(r in 1:length(regs)){
	image(bin_sites[[r]], axes=FALSE, col=fields::tim.colors())
	mtext(regs[r], side=3, line=0, font=2)
	xl <- as.numeric(rownames(bin_sites[[r]]))
	yl <- as.numeric(colnames(bin_sites[[r]]))
	axis(1, at=seq(0,1, length.out=length(xl)), labels=xl)
	axis(2, at=seq(0,1, length.out=length(yl)), labels=yl)
	if(r==length(regs)){
		mtext("Longitude-latitude bin size (degrees)", side=1, line=0, outer=TRUE)
		mtext("Depth bin size (meters)", side=2, line=-0.5, outer=TRUE)
	}
	fields::image.plot(bin_sites[[r]], axes=FALSE, legend.only=TRUE, graphics.reset=TRUE)
}

#+ bin_size_scaled, fig.width=4, fig.height=7, fig.cap="**Bin Size** How Bin size (lon, lat, depth) affects the number of sites in a region. To compare across regions, we can first z-score each region, then take the cross-region average for each bin size combination."
scaled_bin_sites <- lapply(bin_sites, function(x){x2 <- x; x2[,] <- scale(c(x)); x2})
mean_sbs <- apply(simplify2array(scaled_bin_sites), c(1,2), mean, na.rm=TRUE)
min_sbs <- apply(simplify2array(scaled_bin_sites), c(1,2), min, na.rm=TRUE)
par(mfrow=c(2,1), mar=c(4,4,1,5))

image(mean_sbs, axes=FALSE, col=fields::tim.colors())
mtext("cross-region average of z-scores", side=3, line=0, font=2)
xl <- as.numeric(rownames(bin_sites[[1]]))
yl <- as.numeric(colnames(bin_sites[[1]]))
axis(1, at=seq(0,1, length.out=length(xl)), labels=xl)
axis(2, at=seq(0,1, length.out=length(yl)), labels=yl)
mtext("Longitude-latitude bin size (degrees)", side=1, line=2)
mtext("Depth bin size (meters)", side=2, line=2)
fields::image.plot(mean_sbs, axes=FALSE, legend.only=TRUE, graphics.reset=TRUE)

image(min_sbs, axes=FALSE, col=fields::tim.colors())
mtext("cross-region minimum of z-scores", side=3, line=0, font=2)
axis(1, at=seq(0,1, length.out=length(xl)), labels=xl)
axis(2, at=seq(0,1, length.out=length(yl)), labels=yl)
mtext("Longitude-latitude bin size (degrees)", side=1, line=2)
mtext("Depth bin size (meters)", side=2, line=2)
fields::image.plot(min_sbs, axes=FALSE, legend.only=TRUE, graphics.reset=TRUE)


#' ##Tows per site, sites per region
#+ counts_years_sites_tows_spp
sumry_counts <- data_all[reg!="wcann", 
	j = {
		yi <- table(diff(sort(unique(year))))
		yi_form <- paste0(names(yi), paste("(",yi,")",sep=""))
		data.table(
			"Years" = trawlData::lu(year),
			"Year Range" = paste(min(year),max(year),sep=" - "),
			"Year Interval" = paste(yi_form,collapse=", "),
			"Sites" = trawlData::lu(stratum), 
			# .SD[,list("max. Tows"=max(trawlData::lu(haulid)), "Average Tows" = mean(trawlData::lu(haulid))),by=c("stratum","year")]
			"Max. Tows" = .SD[,trawlData::lu(haulid),by=c("stratum","year")][,max(V1)],
			"Average Tows" = .SD[,trawlData::lu(haulid),by=c("stratum","year")][,mean(V1)],
			"Total Species" = trawlData::lu(spp)
		)
	}, 
	by=c('reg')
]
setnames(sumry_counts, 'reg', "Region")
knitr::kable(sumry_counts, caption="Years = Total number of years sampled, Year Range = the minimum and maximum of years sampled, Year Interval = the number of years elapsed between samples (and the frequency of this interval in parentheses), Sites = the total number of sites in the region, Max. Tows = maximum number of tows per site per year, Average Tows = the average number of tows per site per year, Total Species = the total number of species observed in the region across all sites and years.")




#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   