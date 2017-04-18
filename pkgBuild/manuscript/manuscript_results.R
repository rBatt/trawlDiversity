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

interactive()

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
#' #Species Richness
#' ##Richness Time Series
#' ###Table. Summary of richness mean and variability
#+ richnessSummary-basic-table, echo=FALSE
kable(
	data.table(
		comm_master[,list(mu=mean(reg_rich), sd=stats::sd(reg_rich)),by='reg'], 
		comm_master[, list(mu_sd=mean(.SD[,stats::sd(reg_rich),by='reg'][,V1]))]
	), 
	caption="Long-term mean and standard deviation in MSOM richness, and average of those sd's. Gmex and Shelf have the highest and lowest long-term averages in species richness (MSOM)."
)

#' ###Figure 1. MSOM richness time series
#+ Richness-ts-fig, fig.height=5, fig.width=3.5, ecal=TRUE, echo=TRUE, fig.cap="**Figure 1.** Time series of MSOM estimates of region richness. Each point is the posterior mean of regional richness in a year. Lines indicate long-term trends from fitted values of linear regression models predicting richness from time."
load("../results/rich_trend_kendall.RData") # load rich trend stats
richness_ts()


#' ###Figure S1. MSOM - naive scatter
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
#' ##Richness Trend
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
#' ###Table S1. Naive tau
#+ naiveTau-table, echo=FALSE
kable(rich_naive_trend_kendall, 
	caption="**Table S1.** Naive richness trends in each region. Estimate is Kendall's Tau_b, BH is the Benjamini-Hochberg corrected p-value, and p.value is the original p-value."
)
#' In the Naive estimates, `r rich_naive_trend_kendall[reg!='wcann', sum(p.value<=0.05)]` regions had significant $\tau_b$.  
#'   
#' ###Table 2. MSOM tau
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
#' #Colonization and Extinction  
#' ##Colonization and Extinction Summary
#' ###Figure S2.  Barplot of categories
#+ Richness-col-ext-barplot, include=TRUE, echo=TRUE, fig.height=3.5, fig.width=3.5, fig.cap="**Figure S2.** Number of species beloning to the categories of both, neither, colonizer, leaver in each region"
categ_barplot()
#+ Richness-col-ext-barplot-table, echo=FALSE
categ_tbl <- t(spp_master[!duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ)])[c(4,1,2,3),]
kable(categ_tbl, caption = "Number of species in each category in each region.")
#' It's the same pattern, whichever way you split it. However, AI is the only region that had more *colonizers* than *both* species. An interesting way to think about some of this is that the average sd in richness was `r comm_master[,stats::sd(reg_rich),by='reg'][,mean(V1)]`, so when the number of *colonizer* or *leaver* species exceed's that region's sd, the impact of those categories, which I consider to be dubious, might start being relevant (though it's not necessarily problematic, nor is this even close to an actual test for the significance of those categories to the trend). EBS and Shelf had significant positive trends in richness and very low numbers in the *colonizer* category. WCTRI and NEWF had similar numbers in the *both* and *colonizer* category.  

#' ###Table Attributing Richness Change to Colonizers (only)
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

att_categ[,reg:=pretty_reg[reg]] # change to pretty region names
setnames(att_categ, old=c("delPred","reg"), c("Richness Change","Region")) # rename columns
att_categ_print <- att_categ[,c("Region", "Richness Change", "colonizer"), with=FALSE] # change order for printing
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
#'   
#' ###Figure NotIncluded. Time series of colonizations and extinctions
#+ Richness-col-ext-ts, include=TRUE, echo=TRUE, fig.height=3.5, fig.width=3.5, fig.cap="**Figure NotIncluded.** Number of colonizations (blue) and extinctions (red) over time in each region."
col_ext_ts()
#' In most regions the differences in colonization and extinction numbers are similar over time. The most obvious exceptions are for the 3 regions that showed large initial spikes in richness; the GOA, GMEX, and AI regions initially have much larger numbers of colonizers than leavers, but this number shrinks rapidly until the two rates are ~equal.  
#'   
#' For the regions with significant positive slopes, there is no visually obvious increase in colonizations relative to extinctions over time. Because the colonization and extinction numbers tend to track each other over the long-term, it it would be difficult to attribute the long-term changes in richness to a change in just colonization or extinction rates.
#'   
#' **Manuscript paragraph:**  
#' A time series of richness can be decomposed into the colonizations and extinctions of individual species over time. We categorized species according to the following colonization extinction patterns: present in all years = neither (536 species), colonized and went extinct = both (263 species),  initially absent but present every year after its colonization = colonizer (61 species), initially present but absent every year after its extinction = leaver (4 species). Most regions had the same overall ranking (neither > both > colonizer > leaver), except in the Northeast US where both was the most common and neither second, and in the Aleutian Islands where colonizer was the second most common and both third (**Figure S2**). In general, changes in richness were not due to permanent departures or introductions of species to the region. Furthermore, colonization and extinction rates did not become more dissimilar over time for any region (**Figure NotIncluded**). Colonizations were initially greater than extinctions in Aleutian Islands, Gulf of Alaska, and Gulf of Mexico, but this difference disappeared in the latter portion of these time series, as evidenced by these regions’ initially rapid increase in richness that later plateaued. The other regions did not show strong trends in the difference between colonizations and extinctions over time, making it difficult to attribute the long-term trends in richness to changes in just colonization or just extinction rates.
#'   
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Richness and Geographic Range
#' ##Richness and Range Density, Range Size
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
#' ###Figure 3. Species richness versus geographic range size
#+ rich-geo-rangeSize, fig.width=3.5, fig.height=3.5, fig.cap="**Figure 3.** Species richness vs geographic range size. Range size is presented as each species' long-term average of the proportion of sites it occupied. Solid lines are linear regressions with MSOM richness as the response and the horizontal axis and an intercept as the predictors."
rich_geoRange(rSMet_ltComm, leg=TRUE, legPan=1, panLab=FALSE)
print(rSMet_ltComm)

#' Range size is a pretty good predictor of species richness. I think I had originally missed the range size relationship b/c I hadn't done the same aggregating procedure. The interpretation I have is that richness is highest when you have a bunch of rare species.  
#'   
#' The goal here is to see if species richness is predicted by the typical range size of community's constituent species. First I'll run different types of models just to explore whether this is true, in general (across regions). Then I'll drill in to each region individually to answer the same question.  
#'   
#' ###Table. Regressions relating richness to range size
#+ rich-rangeSize
range_reg <- make_range_reg(dens="propTow_occ_avg", size="range_size_mu_avg_ltAvg")

# Fit different models to the whole data set
# Models vary in which/ how parameters vary among regions
rSize_mods <- list()
rSize_mods[[1]] <- lm(rich ~ size, data=range_reg) # simple
rSize_mods[[2]] <- lm(rich ~ size*reg, data=range_reg) # slope and intercept ind. among regs
rSize_mods[[3]] <- lme4::lmer(rich ~ size + (1|reg), data=range_reg) # intercept varies randomly among regs
rSize_mods[[4]] <- lme4::lmer(rich ~ size + (size|reg), data=range_reg) # slope & int. vary randomly among regs
rich_size_smry <- smry_modList(rSize_mods) # summarize to get R^2, AIC, p-vals, etc

# Fit simple model to each region separately 
rSize_reg_mods <- list()
ur <- range_reg[,unique(reg)]
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
#' ##Predicting Richness: Range Size or Density?
#' As far as picking one or the other, it doesn't end up mattering much. Range size is a lot better than density in NEUS, and density outperforms size in AI. Otherwise, size as a slight edge over density on average, although both predictors are significant in all regions.  
#'   
#' ##Conclusion for Richness and Geographic Range
#' The two metrics of geographic range are well correlated. Furthermore, richness can be predicted pretty well using regressions with either as a predictor. There are large differences among regions, though. This is probably because richness is not readily comparable among most regions. Regions vary mostly in their intercept values, and they have fairly similar slopes (though they are not identical, and model fits improve when allowing slopes to vary among regions; it's just that the improvement is small compared to allowing intercepts to vary among regions).  
#'  
#' The interpretation of the result that geographic distribution predicts species richness is likely associated with species rarity. When the average range density or range size of a community is low, it means it has a lot of species that are rare (at either spatial scale). It's these rare species that come and go, and form the dynamics of richness that we observe. When that dynamical value is high, it implies that an above-average number of the dynamic species are present. Because those species are transient (dynamic), they are also probably rare.  
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' #Geographic Range and Colonizations and Extinctions
#'   
#' ##Total Colonizations vs Extinctions and Geographic Range
#' The result that richness is predicted by geographic range implied an underlying association between range, colonization/ extinction, and richness itself. Above, richness was explained with range. Later, the results will explain how range changes near a colonization/ extinction event. Here, between the two, the results show how the number of colonizations is related to range.  
#'   
#' ###Figure S7. Total Colonizations/ Extinctions vs Georaphic Range
#+ ceEvents-vs-rangeSize, fig.width=3.5, fig.height=3.5, fig.cap="**Figure S7.** Number of colonizations and extinctions as a function of range size and range density."
ceEventRange("mean_size")
#' Yup, this is definitely a thing. Long-term average range size and range density predict how many colonizations and extinctions a species is likely to have. This will lead nicely into examing how range changes prior to an extinction or after a colonization.  
#'   
#' ##Range Change after Colonization/ before Extinction
#' ###Figure 4. Changes in Geographic Range before Extinction and after Colonization
#+ rangeSize_ColExt, fig.width=6, fig.height=3, fig.cap="**Figure 4.** Geographic range size vs years until extinction (A) and years after colonization (B). For visualization purposes, range size is averaged across species for each unique value on each axis, and a linear model fit through this average. Statistics in main text use unaggregated data. The horizontal axes were formulated as time until (since) the nearest upcoming (previous) absence. Because range size must be zero when either horizontal axis has a value of zero, points at (0,0) were excluded from figures and analyses."
rangeSize_absenceTime(rSMet_base)
print(rSMet_base)
#' ###Figure 4b. Changes in Geographic Range before Extinction and after Colonization
#+ rangeDensity_ColExt, fig.width=6, fig.height=3, fig.cap="**Figure 4b.** Geographic range density vs years until extinction (A) and years after colonization (B). For visualization purposes, range density is averaged across species for each unique value on each axis, and a linear model fit through this average. Statistics in main text use unaggregated data. The horizontal axes were formulated as time until (since) the nearest upcoming (previous) absence. Because range density must be zero when either horizontal axis has a value of zero, points at (0,0) were excluded from figures and analyses."
rangeSize_absenceTime("rangeDensity")
#' Range size declines near an absence much more consistently than does range density; both are (relatively) low just before extinction and just after colonization. However, range density has much more variable intercepts among regions, whereas range size does not.   
#'   
#' This makes sense, at least somewhat, because colonization and extinction events are defined at the site level; though the outcome isn't necessitated by this formulation, because size could drop suddenly. In fact, when a species is absent, both its range size and its range density must be 0 (though, range density is technically calculated for only those sites that are occupied, so I supposed it's technically undefined according to the equations I'm using).  
#'   
#' I think the regressions for range size should omit an intercept, while the regressions for range density should have it. This might be hard to justify fully *a priori* (though see my thinking in previous paragraph), so I'll probably just do a model selection and maybe discuess the difference if one model has an intercept and the other does not.  
#'   
#' ###Table. Regressions of Range Size & Time to Event; All Regions at Once
#+ rangeSize-ColExtTime-data
rangeTimeDT <- make_rangeTime(sizeName=rSMet_base)

#+ rangeSize-ColExtTime-models, results='markup'
# models for range size
sizeCE_mods <- list()
sizeCE_mods[[1]] <- lme4::lmer(size ~ time + (time|spp/reg), data=rangeTimeDT)
sizeCE_mods[[2]] <- lme4::lmer(size ~ time*type + (time|spp/reg), data=rangeTimeDT)

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
#' ###Table. Regressions of Range Size & Time to Event; Each Region Separate Model (Simple)
#+ rangeSizeDensity-ColExtTime-reg-simple
# Set up -- region names and empty named lists
ur <- rangeTimeDT[,unique(reg)]
sTime_reg_mods1 <- structure(vector("list",length(ur)), .Names=ur)
sTime_reg_mods2 <- structure(vector("list",length(ur)), .Names=ur)
sTime_reg_mods3 <- structure(vector("list",length(ur)), .Names=ur)

# Fit 3 models for each region
for(r in 1:length(ur)){
	t_reg <- ur[r]
	t_dat <- rangeTimeDT[reg==t_reg]
	sTime_reg_mods1[[t_reg]] <- lme4::lmer(size ~ time + (time|spp), data=t_dat)
	sTime_reg_mods2[[t_reg]] <- lme4::lmer(size ~ time + type + (time|spp), data=t_dat)
	sTime_reg_mods3[[t_reg]] <- lme4::lmer(size ~ time * type + (time|spp), data=t_dat)
}

# Summarize each model
sTime_reg_smry1 <- smry_modList2(sTime_reg_mods1)
sTime_reg_smry2 <- smry_modList2(sTime_reg_mods2)
sTime_reg_smry3 <- smry_modList2(sTime_reg_mods3)

# Captions for each model
# Also a caption in case I want to report all models for all regions
capS1 <- sTime_reg_smry1[,unique(mod_call)]
capS2 <- sTime_reg_smry2[,unique(mod_call)]
capS3 <- sTime_reg_smry3[,unique(mod_call)]
timeRegMods_cap <- "Summary statistics for fits of predicting range size from years before extinction and years after colonization. These are mixed effect models of the form "

# Compare AIC values from each model fit
compAIC <- rbind(
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
allSmry <- rbind(sTime_reg_smry1, sTime_reg_smry2, sTime_reg_smry3, fill=TRUE) # combine results from all models and regions
bestModels <- allSmry[bestEachOverall, on=c('reg','mod_call')] # for each region, only select models that were best in that region or best overall

# remove columns that are NA for all rows b/c models w/ those parameters were never winners
loserNames <- sapply(bestModels, function(x)all(is.na(x))) # names of columns only pertaining to non-best models
loserNames <- names(loserNames)[loserNames] # okay, get the names for real, not just logic
bestModels[,c(loserNames):=list(NULL)] # drop the loser columns


#+ rangeSizeDensity-ColExtTime-reg-compareModels, echo=FALSE
# kable(bestEach, caption="Which models for predicting range size from time to event were best in each region?")
# kable(compAIC, caption="Comparing fit of models of varying complexity; more complex models test for differences in the slopes or intercepts of pre-extinction and post-colonization trends.")
kable(bestModels, caption="Shows each region's best model, and shows the model that was most often the best for all regions (including those for which it wasn't the best). If all regions had the same best model, only 1 model is shown.")


#+ rangeSize-timeUntil-detailedFigure, fig.width=6.8, fig.height=6.8
par(mfrow=c(3,3), mgp=c(0.85,0.2,0), ps=8, cex=1, mar=c(1.75,1.75,0.75,0.5), oma=c(0.5,0.5,0.1,0.1), tcl=-0.15)
for(r in 1:length(ur)){
	t_reg <- ur[r]
	t_dat <- rangeTimeDT[reg==t_reg]
	t_dat[,spp:=paste0(spp,event)]
	t_dat[type=="pre_ext", c("spp","time"):=list(spp=paste0(spp,type),time=-time)]
	setkey(t_dat, spp, type, time)
	scatterLine(t_dat, x="time", y="size", lineBy="spp", colBy=adjustcolor(pretty_col[t_reg], 1), lwdBy=0.5, type='p', xlab="", ylab="")
	t_dat[,points(mean(time), mean(size), pch=19, col='black', cex=1.2),by=c('time')]
	t_dat[,points(mean(time), mean(size), bg=pretty_col[t_reg], pch=21, col='white'),by=c('time')]
	mtext(pretty_reg[t_reg], line=-0.2, side=3, adj=0.1, font=2)
}
mtext("Time to Event", side=1, line=-0.5, outer=TRUE)
mtext("Range Size", side=2, line=-0.5, outer=TRUE)


#' ##Conclusion for Range Change before Extinction/ after Colonization
#' As stated with the larger (pooled regional data) regressions, range size responds more consistently/ strongly to an approaching extinction/ departure from colonization than does range density. Range size shrinks as extinction approaches, increases as time since colonization increases. In both cases, there was only 1 region for which the direction (into extinction, out of colonization) mattered: for range density, it was the Southeast US (effect of pre-extinction direction = 0.037, p=0.043 [BH correction], Table 12), and for range size it was Scotial Shelf (effect = -0.021, p = 0.002). Time was a significant predictor of range size in all regions; time was *not* a significant predictor of range density in Newfoundland (p = 0.846 [BH]), and the Scotial Shelf (p = 0.916).  
#'  
#' These results support the hypothesis that species rarity is closely associated with proximity to extinction/ colonization, and therefore the probability that the species will contribute to a change in species richness. Furthermore, these relationships have not been directly quantified for entire assemblages. Competing hypotheses exist regarding how range should change approaching extinction and after colonization. We find that the change in range is similar regardless of direction. However, spatial scale was important. Although slopes were similar, the intercept for density was far larger than for range size --- the fraction of tows containing a species in occupied sites may remain high even when extinction is imminent, and may similar be large even if colonization was recent. Range size, on the other hand, was much closer to 0 near an absence.  
#'   
#' Overall, the results suggest that the spatial footprint of individual species is important for understanding changes in species richness. Furthermore, because species contributing most to the dynamics of richness are those that repeatedly colonize and go extinct, it is meaningful to look at a species' long-term rarity in order to gauge whether it is likely to contribute to those long-term richness changes. Determining what drives the geographic range of individual species is probably a powerful way to anticipate richness changes.  
#'   
#' ###Exploring Range size vs absence time as individual regressions
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
o[,j={plot(nTime,Pr); lines(spline(Pr~nTime),col='red',lwd=2)}]
abline(h=0.05, lwd=1.5)
abline(h=0.05, col='white', lwd=0.5)

# slope vs number of years in series
o[,j={plot(nTime,Estimate); lines(spline(Estimate~nTime),col='red',lwd=2)}]
abline(h=0, lwd=1.5)
abline(h=0, col='white', lwd=0.5)

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


#' ##Table of Tows per Site, Sites per Region
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

reg_depthStratum <- trim_msom_settings('depth')
reg_tolFraction <- trim_msom_settings('tolerance')
yr_subs <- trim_msom_settings('years_logic')
yr_ablin <- trim_msom_settings('years_cutoff')
regs <- names(reg_tolFraction)
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
plot_excludeYearsStrata(TRIM=FALSE)


#+ excluding_years_strata_part2_analyzedStrata, fig.width=5.5, fig.height=10, fig.cap="**Trimming Years b/c of Stratum Locations, Stats for the Analyzed Strata** Changes in stratum statistics (number of strata, min/max lon/lat) for strata included in results in paper."
plot_excludeYearsStrata(TRIM=TRUE)

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
#' ##All Years All Species Present or Absent
#+ pres-abs-allSpp-time, fig.width=5, fig.height=7
# pdf("~/Desktop/pres-abs-allSpp-time.pdf",width=5, height=7)
par(mfrow=c(3,3)) # remove this line for separate images; names then readable
ureg <- spp_master[,unique(reg)]
for(r in 1:length(ureg)){
	t_table <- spp_master[reg==ureg[r] & present==1, table(spp,year)]
	t_table2 <- t(t_table)#[ncol(t_table):1,]
	t_table3 <- t_table2[,order(colSums(t_table2))]

	par(mar=c(1,5,0.5,1), ps=8, mgp=c(0.75,0.2,0), tcl=-0.15)
	image(t_table3, axes=FALSE)
	abline(h=seq(0,1,length.out=ncol(t_table3)), v=seq(0,1,length.out=nrow(t_table3)), col='gray', lty='dotted', lwd=0.5)
	axis(side=2, at=seq(0,1,length.out=ncol(t_table3)), label=colnames(t_table3), las=1, cex.axis=0.5)
	axis(side=1, at=seq(0,1,length.out=nrow(t_table3)), label=rownames(t_table3))
	text(0.95,0.95, label=ureg[r], font=2)
}
# dev.off()
#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Alpha, beta, and gamma diversity (Naive)
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


#' ##Alpha, beta, and gamma diversity (MSOM)
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
ptCol <- bquote(adjustcolor(pretty_col[reg],0.5))
scatterLine(Data=cmNN, x="range_size_mu_avg_ltAvg", y="reg_rich", lineBy="reg", col=ptCol, pch=20)
scatterLine(Data=cmNN, x="range_size_mu_avg", y="reg_rich", lineBy="reg", col=ptCol, pch=20)
scatterLine(Data=cmNN, x="propStrata_noNeither_avg_ltAvg", y="reg_rich", lineBy="reg", col=ptCol, pch=20)
scatterLine(Data=cmNN, x="propStrata_noNeither_avg", y="reg_rich", lineBy="reg", col=ptCol, pch=20)

#' 
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' ##Range Size Over Time For Community vs Transients
#+ rareExpansion, fig.width=6.5, fig.height=6.5, results='asis'
eval(neitherPres_subColor())
plot_rangeSize_FullTrans(range_type="range_size_samp")

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


#+ rareExpansion-description, fig.width=6.5, fig.height=6.5, results='markup'
#' Examining how range size changed over time for "core" and "transient" groups of species. Core species are those that are present in all years. Transient species are those that are not present in all years. **`r transExpansionStats[,sum(BH_year<0.05)]`** regions had significant trends in range size for transient species (average slope for regions with significant trend = **`r transExpansionStats[,mean(year[BH_year<0.05])]`**), and **`r coreExpansionStats[,sum(BH_year<0.05)]`** regions had significant trends in range size for core species (**`r coreExpansionStats[,paste(reg[(BH_year>=0.05)], collapse = ", ")]`** did not have significant trends for core species; averge slope for regions with significant trend = **`r coreExpansionStats[,mean(year[BH_year<0.05])]`**).  
#'   
#' Given that many regions showed positive trends in range size for both species groups, it is prudent to compare trends between these groups, to see if they differ. To this end, I analyzed the groups together, and tested the significance of an interaction term that adjusted the a baseline slope if the species was in the transient group.  
#'   
#' Note that transient species are very much expected to have a smaller range size than the core group, particularly at the beginning of the time series; I'm not dynamically reporting the significance of the difference in intercepts here, but at the time of writing all regions showed significant differences in intercept between groups, with the transient group having a smaller intercept than the core group.
#'   
#' Interestingly, **`r ctExpansionStats[,sum(BH_year<0.05)]`** regions had significant trends (main effect), and **`r ctExpansionStats[,paste(reg[(BH_year<0.05 & year <0)], collapse = ", ")]`** had significant negative trends, while **`r ctExpansionStats[,paste(reg[(BH_year<0.05 & year >0)], collapse = ", ")]`** had significant positive trends. **`r ctExpansionStats[,paste(reg[(get("BH_year:categ2")>=0.05)], collapse = ", ")]`** were the regions without significant interactions -- (at the time of this writing) these are not the same regions that lacked significant trends in the core group. If we sum the main effect and transient-interaction together (ignoring significance), we see that **`r ctExpansionStats[,paste(reg[(year+get("year:categ2transient"))<0], collapse = ", ")]`** was/were the only region/s with negative estimated trends for transient species; all others had positive trends in range for transients. Some regions had significant trends for both core and transients, but the direction of these trends differed between the two groups: **`r ctExpansionStats[,paste(reg[(BH_year<0.05 & year <0 & get("BH_year:categ2")< 0.05 & get("year:categ2transient")>0)], collapse = ", ")]`** had core species with trends in range size that were significantly less than 0 but transient species with trends in range size that were significantly greater than 0.



#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#' ##Range Size and Extinction Probability
#' ###Statistics for Range Size Year before Extinction
#+ rangeProbExt-modelFit, results='markup'
std_range_t1 <- spp_master[ce_categ!='neither' & ce_categ!="colonizer", j={
	
	std_range <- .SD[,list(year, ext, std_range=c(scale(get(rSMet_base))), ext_dist=ext_dist),by=c("spp")]
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
(smry_modList(rangeExtProb1, pred_name="std_range"))

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
print(std_range_t1[,{print(summary(glm(ext_dist~std_range, family='poisson')));NULL},by='reg'])
print(std_range_t1[,{print(summary(glmer(ext~std_range+(std_range|reg), family='poisson')));NULL}])

# Set up -- region names and empty named lists
ur <- std_range_t1[,unique(reg)]
rangeExtProb1_2 <- structure(vector("list",length(ur)), .Names=ur)
rangeExtProb2_2 <- structure(vector("list",length(ur)), .Names=ur)
rangeExtProb3_2 <- structure(vector("list",length(ur)), .Names=ur)

# Fit 3 models for each region
for(r in 1:length(ur)){
	t_reg <- ur[r]
	t_dat <- std_range_t1[reg==t_reg]
	rangeExtProb1_2[[t_reg]] <- (glm(ext_dist~std_range, family='poisson', data=t_dat))
	rangeExtProb2_2[[t_reg]] <- (glmer(ext_dist~std_range+(std_range|spp), family='poisson', data=t_dat))
	# rangeExtProb3_2[[t_reg]] <- summary(glm(ext_dist~std_range*(std_range|spp), family='poisson', data=t_dat))
}
(rangeExtProb1_2_smry <- smry_modList(rangeExtProb1_2, pred_name="std_range"))
(rangeExtProb2_2_smry <- smry_modList2(rangeExtProb2_2))

#' ###Figure for Range Size Predicting Years until Extinction
#+ rangeProbExt-modelFit2-fig, fig.width=6, fig.height=6
par(mfrow=c(3,3), mar=c(2.5,2.5,1.5,0.5), ps=10, mgp=c(0.75,0.15,0), tcl=-0.15)
for(r in 1:length(ur)){
	std_range_t1[reg==ur[r],plot(std_range, ext_dist)]
	new_r <- std_range_t1[reg==ur[r],sort(std_range)]
	lines(new_r, predict(rangeExtProb1_2[[ur[r]]], newdata=data.frame(std_range=new_r), type='response'))
	mtext(pretty_reg[ur[r]], side=3, line=0, adj=0.1, font=2)
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
#'   
#' \FloatBarrier  
#'   
#' ***  
#'   
#'   
#' #Time Series of Tows per Site
#+ towsPerSite-timeSeries, width=6.5, height=6.5, caption="Cross-site distribution (boxplot) of the number of tows (per site) over time, for each region."
plot_towsPerSiteTS()


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