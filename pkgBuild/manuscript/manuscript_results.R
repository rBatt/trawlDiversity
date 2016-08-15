#' ---
#' title: "Shifting distributions and long-term changes in marine assemblage richness"
#' author: "Ryan Batt"
#' date: "2015-08-23"
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

# Report
library(knitr)
library(rmarkdown)
library(xtable)
library(kfigr)
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
#' 
#' ##Colonization and Extinction  
#+ Richness-col-ext, include=TRUE, echo=FALSE
col_ext_bp_cap <- "**Figure S2.** Number of species beloning to the categories of both, neither, colonizer, leaver in each region"
col_ext_ts_cap <- "**Figure S3.** Number of colonizations (blue) and extinctions (red) over time in each region."
#'   
#' ####Figure S2.  Barplot of categories
#+ Richness-col-ext-barplot, include=TRUE, echo=TRUE, fig.height=5, fig.width=3.5, fig.cap=col_ext_bp_cap
categ_barplot()
#' And here's the table behind the barplot:  
kable(t(
	spp_master[!duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ)]
)[c(4,1,2,3),])
#'   
#' Or, to see the most common categories overall:  
kable(data.frame(as.list(
	rowSums(spp_master[!duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ)][,c(4,1,2,3)])
)))
#'   
#' Does that pattern change if we just look at regions with positive slopes? 
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
#' No, it's the same pattern, whichever way you split it. However, AI is the only region that had more *colonizers* than *both* species. An interesting way to think about some of this is that the average sd in richness was `r comm_master[,stats::sd(reg_rich),by='reg'][,mean(V1)]`, so when the number of *colonizer* or *leaver* species exceed's that region's sd, the impact of those categories, which I consider to be dubious, might start being relevant (though it's not necessarily problematic, nor is this even close to an actual test for the significance of those categories to the trend). EBS and Shelf had significant positive trends in richness and very low numbers in the *colonizer* category. WCTRI and NEWF had similar numbers in the *both* and *colonizer* category.  
#'   
#' ####Figure S3. Time series of colonizations and extinctions
#+ Richness-col-ext-ts, include=TRUE, echo=TRUE, fig.height=3.5, fig.width=3.5, fig.cap=col_ext_ts_cap
col_ext_ts()
#' In most regions the differences in colonization and extinction numbers are similar over time. The most obvious exceptions are for the 3 regions that showed large initial spikes in richness; the GOA, GMEX, and AI regions initially have much larger numbers of colonizers than leavers, but this number shrinks rapidly until the two rates are ~equal.  
#'   
#' For the regions with significant positive slopes, there is no visually obvious increase in colonizations relative to extinctions over time. Because the colonization and extinction numbers tend to track each other over the long-term, it it would be difficult to attribute the long-term changes in richness to a change in just colonization or extinction rates.
#'   
#' **Manuscript paragraph:**  
#' A time series of richness can be decomposed into the colonizations and extinctions of individual species over time. We categorized species according to the following colonization extinction patterns: present in all years = neither (536 species), colonized and went extinct = both (263 species),  initially absent but present every year after its colonization = colonizer (61 species), initially present but absent every year after its extinction = leaver (4 species). Most regions had the same overall ranking (neither > both > colonizer > leaver), except in the Northeast US where both was the most common and neither second, and in the Aleutian Islands where colonizer was the second most common and both third (Figure S2). In general, changes in richness were not due to permanent departures or introductions of species to the region. Furthermore, colonization and extinction rates did not become more dissimilar over time for any region (Figure S3). Colonizations were initially greater than extinctions in Aleutian Islands, Gulf of Alaska, and Gulf of Mexico, but this difference disappeared in the latter portion of these time series, as evidenced by these regions’ initially rapid increase in richness that later plateaued. The other regions did not show strong trends in the difference between colonizations and extinctions over time, making it difficult to attribute the long-term trends in richness to changes in just colonization or just extinction rates.
#'   
#' \FloatBarrier   
#' 
#' ***  
#' 
#' ##Spatial Clustering of Colonization and Extinction
#+ col-ext-maps-cap, echo=FALSE
col_map_cap <- "**Figure 2.** Maps of long-term averages of colonizations per site per decade for each region: A) E. Bering Sea, B) Gulf of Alaska, C) Aleutian Islands, D) Scotian Shelf, E) West Coast US, F) Newfoundland, G) Gulf of Mexico, H) Northeast US, I) Southeast US. Values of colonization rate were smoothed using a Gaussian kernel smoother. The smoothed colonization rate is indicated by the color bars in each panel; colors are scaled independently for each region."
ext_map_cap <- "**Figure S4.** Maps of long-term averages of extinctions per site per decade for each region: A) E. Bering Sea, B) Gulf of Alaska, C) Aleutian Islands, D) Scotian Shelf, E) West Coast US, F) Newfoundland, G) Gulf of Mexico, H) Northeast US, I) Southeast US. Values of extinction rate were smoothed using a Gaussian kernel smoother. The smoothed extinction rate is indicated by the color bars in each panel; colors are scaled independently for each region."
col_nb_cap <- "**Figure S5.** Connectivity and local spatial autocorrelation of colonization events in each region. Each site is represented by a point. Points connected by a line are neighbors. For each region, neighbors were determined by first calculating the minimum distance required to allow each site to have at least 1 neighbor. Neighbors of a focal point were then defined as the points within this minimum distance from the focal point. Local spatial autocorrelation is local Moran’s I, significant LMI is indicated by a solid point, the color of which indicates the value of the LMI statistic. The outline is the region boundary used for smoothing in Figure 3 (main text), but does not affect calculations of LMI."
ext_nb_cap <- "**Figure S6.** Connectivity and local spatial autocorrelation of extinction events in each region. Each site is represented by a point. Points connected by a line are neighbors. For each region, neighbors were determined by first calculating the minimum distance required to allow each site to have at least 1 neighbor. Neighbors of a focal point were then defined as the points within this minimum distance from the focal point. Local spatial autocorrelation is local Moran’s I, significant LMI is indicated by a solid point, the color of which indicates the value of the LMI statistic. The outline is the region boundary used for smoothing in Figure 3 (main text), but does not affect calculations of LMI."
#'   
#' ####Figure 2. Colonization map
#+ col-map, echo=TRUE, fig.width=7, fig.height=3, fig.cap=col_map_cap
ceRate_map(ce="colonization")
#'   
#' ####Figure S4. Extinction map
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


#' ####Figure S5. Colonization neighborhood
#+ col-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap=col_nb_cap
nb_moranI(ce="colonization")
#'   
#' ####Figure S6. Extinction neighborhood
#+ ext-nb, echo=TRUE, fig.width=7, fig.height=3, fig.cap=ext_nb_cap
nb_moranI(ce="extinction")
#'   
#' \FloatBarrier   
#' 
#' ***  
#' 
##+ predict-rangeSize

##+ predict-rangeSize




