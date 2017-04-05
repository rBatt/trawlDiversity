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
#+ attribute_richness_categ, echo=FALSE
categ_tbl <- t(spp_master[!duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ)])[c(4,1,2,3),]
predRich <- bquote(predict(lm(reg_rich~year)))

categ_tbl2 <- spp_master[!duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ, dnn=NULL)]
class(categ_tbl2) <- "matrix"
dels <- comm_master[,
	list(
		# minRich=min(reg_rich),
		# maxRich=max(reg_rich), 
		rangeRich=diff(range(reg_rich)),
		# minPred=min(eval(predRich)),
		# maxPred=max(eval(predRich)),
		delPred=rev(eval(predRich))[1] - eval(predRich)[1]
	),
	by=c('reg')
]
att_categ <- merge(dels, data.table(reg=rownames(categ_tbl2),categ_tbl2))
att_categ[,attr_col:=(delPred-colonizer)]
att_categ[,attr_ext:=(delPred-leaver)]

setnames(att_categ, 'reg', "Region")
att_categ[,Region:=factor(Region, levels=c("ebs","ai","goa","wctri","gmex","sa","neus","shelf","newf"), labels=c("E. Bering Sea","Aleutian Islands","Gulf of Alaska","West Coast US","Gulf of Mexico","Southeast US", "Northeast US", "Scotian Shelf","Newfoundland"))]
setkey(att_categ, Region)
setnames(att_categ, c("delPred"), "Richness Change")
att_categ_print <- att_categ[,c("Region", "Richness Change", "colonizer"), with=FALSE]
stargazer(att_categ_print, summary=FALSE, rownames=FALSE, column.sep.width="2pt", digits=1, digits.extra=0)
# \begin{table}[!htbp] \centering
#   %\caption{}
#   %\label{}
# \begin{tabular}{@{\extracolsep{2pt}} ccc}
# \\[-1.8ex]\hline
# \hline \\[-1.8ex]
# Region & Richness Change & Colonizing Species \\
# \hline \\[-1.8ex]
# E. Bering Sea & $12.5$ & $2$ \\
# Aleutian Islands & $5.6$ & $9$ \\
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
# range_reg <- comm_master[,list(
# 	reg, year, rich=reg_rich, density=propTow_occ_avg, size=propStrata_avg_ltAvg
# )]
range_reg <- comm_master[,list(
	reg, year, rich=reg_rich, density=propTow_occ_avg, size=range_size_mu_avg_ltAvg
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
#+ ceEvents-vs-rangeSize, fig.width=3.5, fig.height=3.5, fig.cap="**Figure S7.** Number of colonizations and extinctions as a function of range size and range density."
ceEventRange("mean_size")
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
	# size=propStrata,
	size=range_size_mu,
	density=propTow_occ
)]

# summary function used in tables
mod_smry2 <- function(m){
	mod_call <- switch(class(m), lmerMod=m@call, lm=m$call)
	mod_call <- as.character(mod_call)[2]
	fits <- tryCatch(sem.model.fits(m), error=function(cond)NA)
	if(all(is.na(fits))){
		warning("Error in sem.model.fits")
		fits <- data.frame(Class=class(m), Family="gaussian", Link="identity", N=as.numeric(nobs(m)),"Marginal"=NA_real_,"Conditional"=NA_real_) # family always gaussian for lmer and lm
	}else{
		fits[,c("Marginal","Conditional")] <- lapply(fits[,c("Marginal","Conditional")], signif, 4)
	}
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
#' ####Table. Regressions separate regions -- range size vs time until, SIMPLE
#+ rangeSizeDensity-ColExtTime-reg-simple
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

sTime_reg_mods3 <- list()
ur <- range_reg[,unique(reg)]
for(r in 1:length(ur)){
	sTime_reg_mods3[[r]] <- lme4::lmer(size ~ time + (time|spp), data=rangeTimeDT[reg==ur[r]])
}
sTime_reg_smry3 <- data.table(rbind(
	rbindlist(structure(lapply(sTime_reg_mods3, mod_smry2), .Names=ur), idcol=TRUE)
))
setnames(sTime_reg_smry3, old=c(".id","Pr..Chisq.","Marginal","Conditional"), new=c('reg',"pval","MargR2","CondR2"))
setkey(sTime_reg_smry3, reg, Class, mod_call, predictor)
sTime_reg_smry3[, BH:=round(p.adjust(pval, "BH"), 3), by=c("mod_call", "predictor")]# adjust p-values for multiple tests
sTime_reg_smry3 <- dcast(sTime_reg_smry3, reg+Class+mod_call+MargR2+CondR2+AIC~predictor, value.var=c("pval", "BH"))# rearrange so each model on 1 line

modCoef3 <- getCoefs(sTime_reg_mods3)
sTime_reg_smry3 <- merge(sTime_reg_smry3, modCoef3, by=c("reg","mod_call"))
setnames(sTime_reg_smry3, old=c("(Intercept)"), new=c("Int"))

#+ rangeSizeDensity-ColExtTime-reg-simple-table, echo=FALSE
capS3 <- sTime_reg_smry3[grepl("size",mod_call),mod_call[1]]
sTime_reg_smry3 <- sTime_reg_smry3[grepl("size",mod_call)]
sTime_reg_smry3[,c("mod_call","randomGroup","Class"):=NULL]
kable(sTime_reg_smry3, caption=paste0("Summary statistics for fits of predicting range SIZE from years before extinction and years after colonization. These are mixed effect models of the form ", capS3))
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
sTime_reg_smry[, BH:=round(p.adjust(pval, "BH"), 3), by=c("mod_call", "predictor")]# adjust p-values for multiple tests
sTime_reg_smry <- dcast(sTime_reg_smry, reg+Class+mod_call+MargR2+CondR2+AIC~predictor, value.var=c("pval", "BH"))# rearrange so each model on 1 line

# add coefficients to the sTime_reg_smry data.table
modCoef <- getCoefs(sTime_reg_mods)
sTime_reg_smry <- merge(sTime_reg_smry, modCoef, by=c("reg","mod_call"))
setnames(sTime_reg_smry, old=c("(Intercept)"), new=c("Int"))

#+ rangeSizeDensity-ColExtTime-reg-table, echo=FALSE
capS <- sTime_reg_smry[grepl("size",mod_call),mod_call[1]]
sTime_reg_smry <- sTime_reg_smry[grepl("size",mod_call)]
sTime_reg_smry[,c("mod_call","randomGroup","Class"):=NULL]
kable(sTime_reg_smry, caption=paste0("Summary statistics for fits of predicting range SIZE from years before extinction and years after colonization. These are mixed effect models of the form ", capS))
#'   
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
sTime_reg_smry2[, BH:=round(p.adjust(pval, "BH"), 3), by=c("mod_call", "predictor")]# adjust p-values for multiple tests
sTime_reg_smry2 <- dcast(sTime_reg_smry2, reg+Class+mod_call+MargR2+CondR2+AIC~predictor, value.var=c("pval", "BH"))# rearrange so each model on 1 line
modCoef2 <- getCoefs(sTime_reg_mods2)
sTime_reg_smry2 <- merge(sTime_reg_smry2, modCoef2, by=c("reg","mod_call"))
setnames(sTime_reg_smry2, old=c("(Intercept)"), new=c("Int"))

#+ rangeSizeDensity-ColExtTime-reg-interaction-table, echo=FALSE
capS2 <- sTime_reg_smry2[grepl("size",mod_call),mod_call[1]]
sTime_reg_smry2 <- sTime_reg_smry2[grepl("size",mod_call)]
sTime_reg_smry2[,c("mod_call","randomGroup","Class"):=NULL]
kable(sTime_reg_smry2, caption=paste0("Summary statistics for fits of predicting range SIZE from years before extinction and years after colonization. These are mixed effect models of the form ", capS2))

#+ rangeSizeDensity-ColExtTime-reg-compareModels, echo=FALSE
compAIC <- rbind(
	cbind(sTime_reg_smry3[,list(reg,MargR2,CondR2,AIC)], mod=capS3),
	cbind(sTime_reg_smry[,list(reg,MargR2,CondR2,AIC)], mod=capS),
	cbind(sTime_reg_smry2[,list(reg,MargR2,CondR2,AIC)], mod=capS2)
)
setkey(compAIC, reg)
kable(compAIC)



#' ###Conclusion for Range Change before Extinction/ after Colonization
#' As stated with the larger (pooled regional data) regressions, range size responds more consistently/ strongly to an approaching extinction/ departure from colonization than does range density. Range size shrinks as extinction approaches, increases as time since colonization increases. In both cases, there was only 1 region for which the direction (into extinction, out of colonization) mattered: for range density, it was the Southeast US (effect of pre-extinction direction = 0.037, p=0.043 [BH correction], Table 12), and for range size it was Scotial Shelf (effect = -0.021, p = 0.002). Time was a significant predictor of range size in all regions; time was *not* a significant predictor of range density in Newfoundland (p = 0.846 [BH]), and the Scotial Shelf (p = 0.916).  
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

kable(o[,list(propSignificant=(sum(Pr<0.05)/sum(!is.na(Pr))),n=sum(!is.na(Pr))),by=c("reg","type")])
kable(o[,list(propSignificant=(sum(Pr<0.05)/sum(!is.na(Pr))),n=sum(!is.na(Pr))),by=c("reg","type")][,mean(propSignificant),by=c("type")])

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
			"N" = trawlData::lu(year),
			"Years" = paste(min(year),max(year),sep=" - "),
			"Frequency" = paste(yi_form,collapse=", "),
			"Sites" = trawlData::lu(stratum), 
			# .SD[,list("max. Tows"=max(trawlData::lu(haulid)), "Average Tows" = mean(trawlData::lu(haulid))),by=c("stratum","year")]
			"Tows" = paste0(
				.SD[,trawlData::lu(haulid),by=c("stratum","year")][,max(V1)],
				paste0("(",.SD[,trawlData::lu(haulid),by=c("stratum","year")][,round(mean(V1),1)],")")
				),
			# "Max. Tows" = .SD[,trawlData::lu(haulid),by=c("stratum","year")][,max(V1)],
# 			"Avg Tows" = .SD[,trawlData::lu(haulid),by=c("stratum","year")][,mean(V1)],
			"Species" = trawlData::lu(spp)
		)
	}, 
	by=c('reg')
]
setnames(sumry_counts, 'reg', "Region")
sumry_counts[,Region:=factor(Region, levels=c("ebs","ai","goa","wctri","gmex","sa","neus","shelf","newf"), labels=c("E. Bering Sea","Aleutian Islands","Gulf of Alaska","West Coast US","Gulf of Mexico","Southeast US", "Northeast US", "Scotian Shelf","Newfoundland"))]
setkey(sumry_counts, Region)

sumry_counts[,Org:=c("AFSC","AFSC","AFSC","NMFS","GSMFC","SCDNR","NEFSC","DFO","DFO")]
# sumry_counts[,"Avg Tows":=round(eval(s2c("Avg Tows"))[[1]],1)]
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
# Aleutian Islands & $12$ & 1983 - 2014 & 2(5), 3(4), 4(1), 5(1) & $82$ & 12(2.8) & $55$ & $\text{AFSC}^{*}$ \\
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


#' ##Years Excluded & Strata Sampled
#+ excluding_years_strata, fig.width=7, fig.height=7, fig.cap="**Trimming Years b/c of Strata, Counts of Analyzed and Excluded Strata** Along the vertical axis is a stratum ID. Along the horizontal axis is a year of sampling. The colors are binary: the light color indicates that the stratum was sampled in that year, and the red color indicates that the stratum was not sampled. The vertical line indicates that years at the line or to the left (if the line is on the left half of graph) or to the right (if line is on right half of graph) were excluded from the analysis, inclusively. E.g., for E. Bering Sea (ebs), 1982 and 1983 were excluded. Three regions had no years excluded (ai = Aleutian Islands, goa = Gulf of Alaska, wctri West Coast US).  The last three regions had years excluded because of changes in the number or location of strata sampled: gmex = Gulf of Mexico 1983 < years < 2001; newf = Newfoundland 1995 < years; sa = Southeast US 1989 < years. Three regions had years excluded for reasons other than number of strata sampled: ebs = E. Bering Sea years < 1984, years excluded due to massive early increase in number of species sampled each year, change in identification suspected; neus = Northeast US 1981 < years < 2014; shelf = Scotian Shelf years < 2011, bottom temperature not available in final year (used as covariate in occupancy model). Finally, note that in GoA strata had to be present in 100% of years or else they were not included; this was done b/c in 2001 the western-most stratum that was sampled was much further to the east than in other years. Similarly, in E. Bering Sea strata had to be present in all years, otherwise the Northeastern extent of the sampled range was reduced substantially in several years."

reg_depthStratum <- c(
	"ebs" = 500,
	"ai" = 100,
	"goa" = 500,
	"wctri" = 100, 
	"wcann" = 500, 
	"gmex" = 500, 
	"sa" = 500, 
	"neus" = 500, 
	"shelf" = 500, 
	"newf" = 500
)

reg_tolFraction <- c(
	"ebs" = 0,
	"ai" = 0.15,
	"goa" = 0,
	"wctri" = 0.15, 
	"wcann" = 0.15, 
	"gmex" = 0.15, 
	"sa" = 0.15, 
	"neus" = 0.15, 
	"shelf" = 0.15, 
	"newf" = 0.15
)
regs <- names(reg_tolFraction)

yr_subs <- list(
	ebs = bquote((year)>1983),
	ai = NA,
	goa = NA,
	gmex = bquote((year)>1983 & (year)<2001),
	neus = bquote((year)>1981 & (year)<2014),
	newf = bquote((year)>1995),
	sa = bquote((year)>=1990),
	shelf = bquote((year)!=2011),
	wctri = NA
)
yr_ablin <- list(
	ebs = 1983,
	ai = NA,
	goa = NA,
	gmex = c(1983,2001),
	neus = c(1981,2014),
	newf = 1995,
	sa = 1989,
	shelf = 2011,
	wctri = NA
)
yregs <- names(yr_subs)
par(mfrow=c(3, 3), mar=c(2,2,1,0.25), cex=1, mgp=c(1,0.15,0), tcl=-0.15, ps=8)
for(r in 1:length(yregs)){
	
	b <- trim_msom(yregs[r], gridSize=0.5, grid_stratum=TRUE, depthStratum=reg_depthStratum[yregs[r]], tolFraction=0.5, plot=FALSE, cull_show_up=FALSE, trimYears=FALSE)
	
	b_tbl <- b[,table(year,stratum)>0]
	b_tbl <- b_tbl[,order(colSums(b_tbl))]
	image(b_tbl, axes=FALSE)
	at_vals <- seq(0,1,length.out=nrow(b_tbl))
	axis(1, at=at_vals, labels=rownames(b_tbl))
	abline(v=at_vals[rownames(b_tbl)%in%yr_ablin[[r]]])
	mtext(yregs[r], side=3, line=0, adj=0, font=2)
}

#' ##Years Excluded & Strata Excluded, Part 2
#+ excluding_years_strata_part2_allStrata, fig.width=5.5, fig.height=10, fig.cap="**Trimming Years b/c of Stratum Locations, Stats for Analyzed and Excluded Strata**  The number and coordinate extreme of strata in each region if all years had been included in the analysis, and if strata were only removed if they were absent in more than 15% of years. Note that the 15%  tolerance rule used here can differ from the implemented rule in two ways: 1) for regions with years excluded, the increased number of years may change whether or not a stratum meets the 15% cutoff; 2) two of the regions (E Bering Sea and Gulf of Alaska) had a 0% tolerance, therefore the same strata were sampled in each year for these regions (though the changes in coordinate extrema and stratum count illustrate the need for the unique tolerance rule in these regions)."
par(mfrow=c(9, 5), mar=c(2,2,1,0.25), cex=1, mgp=c(1,0.15,0), tcl=-0.15, ps=8)
for(r in 1:length(yregs)){
	
	b <- trim_msom(yregs[r], gridSize=0.5, grid_stratum=TRUE, depthStratum=reg_depthStratum[yregs[r]], tolFraction=0.15, plot=FALSE, cull_show_up=FALSE, trimYears=FALSE)
	plot(b[,list("# strata sampled"=trawlData::lu(stratum)),by=c("year")])
	mtext(yregs[r], side=3, line=0, adj=0, font=2)
	abline(v=yr_ablin[[r]])
	plot(b[,list("min latitude"=min(lat)),by=c("year")])
	mtext(yregs[r], side=3, line=0, adj=0, font=2)
	abline(v=yr_ablin[[r]])
	plot(b[,list("max latitude"=max(lat)),by=c("year")])
	mtext(yregs[r], side=3, line=0, adj=0, font=2)
	abline(v=yr_ablin[[r]])
	plot(b[,list("min longitude"=min(lon)),by=c("year")])
	mtext(yregs[r], side=3, line=0, adj=0, font=2)
	abline(v=yr_ablin[[r]])
	plot(b[,list("max longitude"=max(lon)),by=c("year")])
	mtext(yregs[r], side=3, line=0, adj=0, font=2)
	abline(v=yr_ablin[[r]])
}

#+ excluding_years_strata_part2_analyzedStrata, fig.width=5.5, fig.height=10, fig.cap="**Trimming Years b/c of Stratum Locations, Stats for the Analyzed Strata** Changes in stratum statistics (number of strata, min/max lon/lat) for strata included in results in paper."
par(mfrow=c(9, 5), mar=c(2,2,1,0.25), cex=1, mgp=c(1,0.15,0), tcl=-0.15, ps=8)
for(r in 1:length(yregs)){
	
	# b <- trim_msom(yregs[r], gridSize=0.5, depthStratum=reg_depthStratum[yregs[r]], tolFraction=reg_tolFraction[yregs[r]], grid_stratum=TRUE, plot=FALSE, cull_show_up=FALSE)
	b <- data_all[reg==yregs[r]]
	
	plot(b[,list("# strata sampled"=trawlData::lu(stratum)),by=c("year")])
	mtext(yregs[r], side=3, line=0, adj=0, font=2)
	abline(v=yr_ablin[[r]])
	plot(b[,list("min latitude"=min(lat)),by=c("year")])
	mtext(yregs[r], side=3, line=0, adj=0, font=2)
	abline(v=yr_ablin[[r]])
	plot(b[,list("max latitude"=max(lat)),by=c("year")])
	mtext(yregs[r], side=3, line=0, adj=0, font=2)
	abline(v=yr_ablin[[r]])
	plot(b[,list("min longitude"=min(lon)),by=c("year")])
	mtext(yregs[r], side=3, line=0, adj=0, font=2)
	abline(v=yr_ablin[[r]])
	plot(b[,list("max longitude"=max(lon)),by=c("year")])
	mtext(yregs[r], side=3, line=0, adj=0, font=2)
	abline(v=yr_ablin[[r]])
}

#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' #Revision Exploration
#' ##Are rare species becoming more common?
#+ rare-more-common, fig.width=7, fig.height=7
blah <- spp_master[present==1, j={
	# list(propSlope=summary(lm(propStrata~year))$coeff[2,1], mu_propStrat=mean(propStrata), ce_categ=ce_categ[1])
	list(propSlope=summary(lm(range_size_mu~year))$coeff[2,1], mu_propStrat=mean(range_size_mu), ce_categ=ce_categ[1])
} ,by=c("reg","spp")]

# par(mfrow=c(3,3))
# blah[,plot(mu_propStrat, propSlope, main=reg[1]),by='reg']

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
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Time Series of Local and Regional Richness
#+ alpha-gamma-TimeSeries, fig.width=7, fig.height=7
localR <- data_all[,list(lR=length(unique(spp))),by=c("reg","stratum","year")]
bothR <- merge(localR, comm_master[,list(reg,year,naive_rich,reg_rich)],by=c("reg","year"))
par(mfrow=c(3,3), cex=1, ps=8, mar=c(2,2,2,2), mgp=c(1,0.2,0), tcl=-0.2)
ureg <- bothR[,unique(reg)]
nreg <- length(ureg)
for(r in 1:nreg){
	bothR[reg==ureg[r],list(mu_lR=mean(lR),naive_rich=naive_rich[1],reg_rich=reg_rich[1]),by=c("reg","year")][,j={
		plot(year, mu_lR, type='l')
		mtext(reg[1], side=3, line=0, adj=0.1, font=2)
		# par(new=TRUE)
		# plot(year, reg_rich, type='l', col='blue', xaxt='n', yaxt='n', xlab='',ylab='')
		par(new=TRUE)
		plot(year, naive_rich, type='l', col='red', xaxt='n', yaxt='n', xlab='',ylab='')
		axis(side=4, col='red')
		NULL
	}]
}

# mtext("Observed Regional (red), MSOM Regional (blue), and Mean Local (black) Richness", side=3, line=-0.75, outer=TRUE, font=2)
mtext("Observed Regional (red) and Mean Local (black) Richness", side=3, line=-0.75, outer=TRUE, font=2)
#' In general, local richness and regional richness are similar. There are differences, possibly in some cases the regional slope might be significant whereas local would not be (I haven't checked, though). However, the trends aren't in opposite directions at the two scales.


#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##All Years All Species Present or Absent
#+ pres-abs-allSpp-time, fig.width=5, fig.height=7
# # dev.new(width=5, height=7)
# # pdf("~/Desktop/pres-abs-allSpp-time.pdf",width=5, height=7)
ureg <- spp_master[,unique(reg)]
for(r in 1:length(ureg)){
	t_table <- spp_master[reg==ureg[r] & present==1, table(spp,year)]
	t_table2 <- t(t_table)#[ncol(t_table):1,]
	t_table3 <- t_table2[,order(colSums(t_table2))]

	par(mar=c(1,5,0.5,1), ps=8, mgp=c(0.75,0.2,0), tcl=-0.15)
	image(t_table3, axes=FALSE)
	# grid(ny=ncol(t_table3)+1, nx=nrow(t_table3)+1)
	abline(h=seq(0,1,length.out=ncol(t_table3)), v=seq(0,1,length.out=nrow(t_table3)), col='gray', lty='dotted', lwd=0.5)
	axis(side=2, at=seq(0,1,length.out=ncol(t_table3)), label=colnames(t_table3), las=1, cex.axis=0.5)
	axis(side=1, at=seq(0,1,length.out=nrow(t_table3)), label=rownames(t_table3))
	text(0.95,0.95, label=ureg[r], font=2)
}
# # # dev.off()
#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Relationships between alpha, beta, and gamma (MSOM) diversity
#+ abgDiversity-MSOM, fig.width=3.5, fig.height=7
localR <- data_all[,list(lR=length(unique(spp))),by=c("reg","stratum","year")]
bothR <- merge(localR, comm_master[,list(reg,year,naive_rich,reg_rich)],by=c("reg","year"))
bothR_mu <- bothR[,list(lR_mu=mean(lR)),by=c("reg","year")]
cm2 <- merge(comm_master, bothR_mu)
eval(figure_setup())
# dev.new(width=3.5, height=7)
par(mfrow=c(3,1), mar=c(1.75,1.5,0.25,0.25),mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
ureg <- cm2[,unique(reg)]
nreg <- length(ureg)
cm2[,j={
	plot(reg_rich, beta_div_mu, col=adjustcolor(pretty_col[reg],0.5), pch=16, xlab="Regional Richness", ylab="Beta Diversity")
	for(r in 1:nreg){
		.SD[reg==ureg[r]][order(reg_rich),j={
			lines(reg_rich,predict(lm(beta_div_mu~reg_rich)),col='black')
		}]
	}
}]
cm2[,j={
	plot(lR_mu, beta_div_mu, col=adjustcolor(pretty_col[reg],0.5), pch=16, xlab="Average Local Richness", ylab="Beta Diversity")
	for(r in 1:nreg){
		.SD[reg==ureg[r]][order(reg_rich),j={
			lines(lR_mu,predict(lm(beta_div_mu~lR_mu)),col='black')
		}]
	}
}]
cm2[,legend("topright",ncol=2,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.01, -0.03), bty='n', x.intersp=0.15, y.intersp=0.65)]
cm2[,j={
	plot(reg_rich, lR_mu, col=adjustcolor(pretty_col[reg],0.5), pch=16, ylab="Average Local Richness", xlab="Regional Richness")
	for(r in 1:nreg){
		.SD[reg==ureg[r]][order(reg_rich),j={
			lines(reg_rich,predict(lm(lR_mu~reg_rich)),col='black')
		}]
	}
}]

#' ##Relationships between alpha, beta, and gamma (Naive) diversity
#+ abgDiversity-Naive, fig.width=3.5, fig.height=7
localR <- data_all[,list(lR=length(unique(spp))),by=c("reg","stratum","year")]
bothR <- merge(localR, comm_master[,list(reg,year,naive_rich,reg_rich)],by=c("reg","year"))
bothR_mu <- bothR[,list(lR_mu=mean(lR)),by=c("reg","year")]
cm2 <- merge(comm_master, bothR_mu)
eval(figure_setup())
# dev.new(width=3.5, height=7)
par(mfrow=c(3,1), mar=c(1.75,1.5,0.25,0.25),mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
ureg <- cm2[,unique(reg)]
nreg <- length(ureg)
cm2[,j={
	plot(naive_rich, beta_div_obs, col=adjustcolor(pretty_col[reg],0.5), pch=16, xlab="Observed Regional Richness", ylab="Observed Beta Diversity")
	for(r in 1:nreg){
		.SD[reg==ureg[r]][order(naive_rich),j={
			lines(naive_rich,predict(lm(beta_div_obs~naive_rich)),col='black')
		}]
	}
}]
cm2[,j={
	plot(lR_mu, beta_div_obs, col=adjustcolor(pretty_col[reg],0.5), pch=16, xlab="Average Local Richness", ylab="Observed Beta Diversity")
	for(r in 1:nreg){
		.SD[reg==ureg[r]][order(naive_rich),j={
			lines(lR_mu,predict(lm(beta_div_obs~lR_mu)),col='black')
		}]
	}
}]
cm2[,legend("topright",ncol=2,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.01, -0.03), bty='n', x.intersp=0.15, y.intersp=0.65)]
cm2[,j={
	plot(naive_rich, lR_mu, col=adjustcolor(pretty_col[reg],0.5), pch=16, ylab="Average Local Richness", xlab="Observed Regional Richness")
	for(r in 1:nreg){
		.SD[reg==ureg[r]][order(naive_rich),j={
			lines(naive_rich,predict(lm(lR_mu~naive_rich)),col='black')
		}]
	}
}]

#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Richness-range size, varying long-term and annual, all spp and transients only
#' By transient I mean species that weren't present in every year.  
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
comm_master[,j={
	plot(range_size_mu_avg_ltAvg, reg_rich, col=adjustcolor(pretty_col[reg],0.5), pch=16)
	for(r in 1:nreg){
		.SD[reg==ureg[r]][order(range_size_mu_avg_ltAvg),j={
			lines(range_size_mu_avg_ltAvg, predict(lm(reg_rich~range_size_mu_avg_ltAvg)))
		}]
	}
}]
comm_master[,j={
	plot(range_size_mu_avg, reg_rich, col=adjustcolor(pretty_col[reg],0.5), pch=16)
	for(r in 1:nreg){
		.SD[reg==ureg[r]][order(range_size_mu_avg),j={
			lines(range_size_mu_avg, predict(lm(reg_rich~range_size_mu_avg)))
		}]
	}
}]
cmNN[,j={
	plot(propStrata_noNeither_avg_ltAvg, reg_rich, col=adjustcolor(pretty_col[reg],0.5), pch=16)
	for(r in 1:nreg){
		.SD[reg==ureg[r]][order(propStrata_noNeither_avg_ltAvg),j={
			lines(propStrata_noNeither_avg_ltAvg, predict(lm(reg_rich~propStrata_noNeither_avg_ltAvg)))
		}]
	}
}]
cmNN[,j={
	plot(propStrata_noNeither_avg, reg_rich, col=adjustcolor(pretty_col[reg],0.5), pch=16)
	for(r in 1:nreg){
		.SD[reg==ureg[r]][order(propStrata_noNeither_avg),j={
			lines(propStrata_noNeither_avg, predict(lm(reg_rich~propStrata_noNeither_avg)))
		}]
	}
}]

#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Range Size Over Time For Community vs Transients
#+ rareExpansion, fig.width=6.5, fig.height=6.5

rE_logic1 <- bquote(reg==rr & ce_categ!="neither" & present==1)
rE_logic2 <- bquote(reg==rr & ce_categ!="neither")
rE_logic3 <- bquote(reg==rr & ce_categ=="neither")
colQ <- bquote(adjustcolor(c("pre_ext"="red","post_col"="blue")[stretch_type],0.5))
colQ2 <- bquote(adjustcolor(c("pre_ext"="red","post_col"="blue")[stretch_type],1))

# dev.new(width=6.5, height=6.5)
par(mfrow=c(3,3), mgp=c(0.85,0.2,0), ps=8, cex=1, mar=c(1.75,1.75,0.5,0.5), tcl=-0.15)
ur <- spp_master[,unique(reg)]
for(r in 1:length(ur)){
	rr <- ur[r]
	rEl3 <- spp_master[order(year)][eval(rE_logic3)]
	rEl2 <- spp_master[order(year)][eval(rE_logic2)]
	rEl1 <- rEl2[order(year)][eval(rE_logic1)]
	u_spp <- rEl1[, unique(spp)]

	rEl1[,j={bp_dat <<- boxplot(range_size_samp~year,plot=FALSE);NULL}]
	bp_ylim <- unlist(bp_dat[c("stats")], use.names=FALSE)
	medRange_pres0 <- rEl2[,list(range_size_samp=median(range_size_samp)),by='year']
	medRange_noNeith <- rEl3[,list(range_size_samp=median(range_size_samp)),by='year']
	rEl1_qylim <- rEl1[,quantile(range_size_samp,c(0.25,0.75)),by='year'][,range(V1)]
	rEl3_qylim <- rEl3[,quantile(range_size_samp,c(0.25,0.75)),by='year'][,range(V1)]
	ylim <- range(c(
		# bp_ylim,
		rEl1_qylim,
		rEl3_qylim,
		# medRange_pres0[,range_size_samp],
		medRange_noNeith[,range_size_samp]#,
	))
	# rEl1[,j={boxplot(range_size_samp~year, add=FALSE, at=unique(year), outline=FALSE, axes=TRUE, ylim=ylim); NULL}]
	
	rEl1[,plot(year, range_size_samp, type='n', ylim=ylim)]
	grid()
	
	r11 <- rEl1[,median(range_size_samp),by='year']
	r12 <- rEl1[,quantile(range_size_samp,0.75),by='year']
	r13 <- rEl1[,quantile(range_size_samp,0.25),by='year']
	# lines(r11, lwd=2, col='red')
	# lines(r12, lwd=1, col='red')
	# lines(r13, lwd=1, col='red')
	poly1y <- c(r12[,V1], r13[,rev(V1)])
	poly1x <- c(r12[,year],r13[,rev(year)])
	
	
	r21 <- rEl3[,median(range_size_samp),by='year']
	r22 <- rEl3[,quantile(range_size_samp,0.75),by='year']
	r23 <- rEl3[,quantile(range_size_samp,0.25),by='year']
	# lines(r21, lwd=2, col='blue')
	# lines(r22, lwd=1, col='blue')
	# lines(r23, lwd=1, col='blue')
	poly2y <- c(r22[,V1], r23[,rev(V1)])
	poly2x <- c(r22[,year],r23[,rev(year)])
	
	
	polygon(poly2x, poly2y, col=adjustcolor('blue',0.15), border=NA)
	polygon(poly1x, poly1y, col=adjustcolor('red',0.15), border=NA)
	lines(r21, lwd=2, col='blue')
	lines(r11, lwd=2, col='red')
	comm_master[reg==rr, lines(year,range_size_samp_avg_ltAvg, col='black')]
	mtext(pretty_reg[rr], line=-0.75, side=3, adj=0.1, font=2)
}

#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Time Series of Range Size, Color is Local Richness
#+ rangeTS-colorAlpha, fig.width=6.5, fig.height=6.5
rE_logic1 <- bquote(reg==rr & ce_categ!="neither" & present==1)
rE_logic2 <- bquote(reg==rr & ce_categ!="neither")
rE_logic3 <- bquote(reg==rr & ce_categ=="neither")
colQ <- bquote(adjustcolor(c("pre_ext"="red","post_col"="blue")[stretch_type],0.5))
colQ2 <- bquote(adjustcolor(c("pre_ext"="red","post_col"="blue")[stretch_type],1))

# dev.new(width=6.5, height=6.5)
par(mfrow=c(3,3), mgp=c(0.85,0.2,0), ps=8, cex=1, mar=c(1.75,1.75,0.5,0.5), tcl=-0.15, oma=c(0.25,0.1,0.1,0.1))
ur <- spp_master[,unique(reg)]
for(r in 1:length(ur)){
	rr <- ur[r]
	rEl3 <- spp_master[order(year)][eval(rE_logic3)]
	rEl2 <- spp_master[order(year)][eval(rE_logic2)]
	rEl1 <- rEl2[order(year)][eval(rE_logic1)]
	u_spp <- rEl1[, unique(spp)]
	
	nCols <- rEl1[,length(unique(year))]
	cols <- viridis(nCols)
	nTrans <- rEl1[,colSums(table(spp,year))]
	nTrans <- nTrans[order(as.integer(names(nTrans)))]
	colVec_ind <- cut(nTrans, breaks=nCols)
	colVec <- cols[colVec_ind]
	
	rEl1[,j={bp_dat <<- boxplot(propStrata~year,plot=FALSE);NULL}]
	bp_ylim <- unlist(bp_dat[c("stats")], use.names=FALSE)
	ylim <- range(c(
		bp_ylim
	))
	rEl1[,j={boxplot(propStrata~year, add=FALSE, at=unique(year), col=colVec, outline=FALSE, axes=TRUE, ylim=ylim); NULL}]
	if(rr=="sa"){
		mapLegend(x=0.3, y=0.78, zlim=range(nTrans),cols=cols)
	}else{
		mapLegend(x=0.05, y=0.78, zlim=range(nTrans),cols=cols)
	}
	mtext(pretty_reg[rr], line=-0.75, side=3, adj=0.1, font=2)
}
mtext("Range Size of Transient Species", side=2, line=-0.75, outer=TRUE, font=2)
mtext("Year", side=1, line=-0.75, outer=TRUE, font=2)

#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Time Series of Range Size, Color is Local Richness (MSOM!)
#+ rangeTS-colorAlpha-MSOM, fig.width=6.5, fig.height=6.5
eval(figure_setup())
rE_logic1 <- bquote(reg==rr & ce_categ!="neither" & present==1)
rE_logic2 <- bquote(reg==rr & ce_categ!="neither")
rE_logic3 <- bquote(reg==rr & ce_categ=="neither")
colQ <- bquote(adjustcolor(c("pre_ext"="red","post_col"="blue")[stretch_type],0.5))
colQ2 <- bquote(adjustcolor(c("pre_ext"="red","post_col"="blue")[stretch_type],1))

# dev.new(width=6.5, height=6.5)
par(mfrow=c(3,3), mgp=c(0.85,0.2,0), ps=8, cex=1, mar=c(1.75,1.75,0.5,0.5), tcl=-0.15, oma=c(0.25,0.1,0.1,0.1))
ur <- spp_master[,unique(reg)]
for(r in 1:length(ur)){
	rr <- ur[r]
	rEl3 <- spp_master[order(year)][eval(rE_logic3)]
	rEl2 <- spp_master[order(year)][eval(rE_logic2)]
	rEl1 <- rEl2[order(year)][eval(rE_logic1)]
	u_spp <- rEl1[, unique(spp)]
	
	nCols <- rEl1[,length(unique(year))]
	cols <- viridis(nCols)
	nTrans <- rEl1[,colSums(table(spp,year))]
	nTrans <- nTrans[order(as.integer(names(nTrans)))]
	colVec_ind <- cut(nTrans, breaks=nCols)
	colVec <- cols[colVec_ind]
	
	rEl1[,j={bp_dat <<- boxplot(range_size_mu~year,plot=FALSE);NULL}]
	bp_ylim <- unlist(bp_dat[c("stats")], use.names=FALSE)
	ylim <- range(c(
		bp_ylim
	))
	rEl1[,j={boxplot(range_size_mu~year, add=FALSE, at=unique(year), col=colVec, outline=FALSE, axes=TRUE, ylim=ylim); NULL}]
	if(rr=="sa"){
		mapLegend(x=0.3, y=0.78, zlim=range(nTrans),cols=cols)
	}else{
		mapLegend(x=0.05, y=0.78, zlim=range(nTrans),cols=cols)
	}
	mtext(pretty_reg[rr], line=-0.75, side=3, adj=0.1, font=2)
}
mtext("Range Size (MSOM) of Transient Species", side=2, line=-0.75, outer=TRUE, font=2)
mtext("Year", side=1, line=-0.75, outer=TRUE, font=2)

#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Time Series of Range Size, Color is Local Richness (Resampled)
#+ rangeTS-colorAlpha-samp, fig.width=6.5, fig.height=6.5
eval(figure_setup())
rE_logic1 <- bquote(reg==rr & ce_categ!="neither" & present==1)
rE_logic2 <- bquote(reg==rr & ce_categ!="neither")
rE_logic3 <- bquote(reg==rr & ce_categ=="neither")
colQ <- bquote(adjustcolor(c("pre_ext"="red","post_col"="blue")[stretch_type],0.5))
colQ2 <- bquote(adjustcolor(c("pre_ext"="red","post_col"="blue")[stretch_type],1))

# dev.new(width=6.5, height=6.5)
par(mfrow=c(3,3), mgp=c(0.85,0.2,0), ps=8, cex=1, mar=c(1.75,1.75,0.5,0.5), tcl=-0.15, oma=c(0.25,0.1,0.1,0.1))
ur <- spp_master[,unique(reg)]
for(r in 1:length(ur)){
	rr <- ur[r]
	rEl3 <- spp_master[order(year)][eval(rE_logic3)]
	rEl2 <- spp_master[order(year)][eval(rE_logic2)]
	rEl1 <- rEl2[order(year)][eval(rE_logic1)]
	u_spp <- rEl1[, unique(spp)]
	
	nCols <- rEl1[,length(unique(year))]
	cols <- viridis(nCols)
	nTrans <- rEl1[,colSums(table(spp,year))]
	nTrans <- nTrans[order(as.integer(names(nTrans)))]
	colVec_ind <- cut(nTrans, breaks=nCols)
	colVec <- cols[colVec_ind]
	
	rEl1[,j={bp_dat <<- boxplot(range_size_samp~year,plot=FALSE);NULL}]
	bp_ylim <- unlist(bp_dat[c("stats")], use.names=FALSE)
	ylim <- range(c(
		bp_ylim
	))
	rEl1[,j={boxplot(range_size_samp~year, add=FALSE, at=unique(year), col=colVec, outline=FALSE, axes=TRUE, ylim=ylim); NULL}]
	if(rr=="sa"){
		mapLegend(x=0.3, y=0.78, zlim=range(nTrans),cols=cols)
	}else{
		mapLegend(x=0.05, y=0.78, zlim=range(nTrans),cols=cols)
	}
	mtext(pretty_reg[rr], line=-0.75, side=3, adj=0.1, font=2)
}
mtext("Range Size (MSOM) of Transient Species", side=2, line=-0.75, outer=TRUE, font=2)
mtext("Year", side=1, line=-0.75, outer=TRUE, font=2)

#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' #Time Series of Tows per Site
#+ towsPerSite-timeSeries, width=6.5, height=6.5
par(mfrow=c(3,3), mar=c(4,3.5,1,1))
ureg <- data_all[reg!='wcann',unique(reg)]
nreg <- length(ureg)
for(r in 1:nreg){
	data_all[reg==ureg[r],list(Kmax=unique(Kmax)),by=c("reg","year","stratum")][,j={boxplot(Kmax~year, main=reg[1]);NULL},by='reg']
}
mtext("Tows per site", side=2, line=-1.25, outer=TRUE)

#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   