#' ---
#' title: "Shifting distributions and long-term changes in marine assemblage richness"
#' author: "Ryan Batt"
#' date: "2015-08-23"
#' abstract: |
#'   These are the results and figures.
#' output:
#'   html_document:
#'     toc: true
#'     toc_depth: 2
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
# 	output_dir='~/Documents/School&Work/pinskyPost/trawlDiversity/pkgBuild/manuscript/manuscript_report',
# )

opts_chunk$set(
	fig.path = 'manuscript_report/', 
	cache.path='manuscript_report/',
	echo=TRUE, 
	include=TRUE, 
	cache=TRUE,
	fig.show="hold",
	fig.lp = if(o_f=="html_document"){"**Figure.**"}else{NULL}
)


# setwd("~/Documents/School&Work/pinskyPost/trawlDiversity/manuscript")
source("../manuscript/manuscript_figures_functions.R")


# ============
# = Richness =
# ============
#+ Richness-basic, include=TRUE, echo=TRUE, eval=TRUE
#' ###Richness Summary
#' Give a basic summary of species richness stats. For example, the regions with the lowest and highest long-term averages in species richness are:  
comm_master[,mean(reg_rich), by='reg'][reg%in%c("gmex","shelf")]

# ---- long-term variability ----
#' Can also measure long-term variability of the different regions. Here is the standard deviation for each region:
comm_master[,stats::sd(reg_rich),by='reg']#[,mean(V1)]
#' And here is the cross-region average of those long-term standard deviations:
comm_master[,stats::sd(reg_rich),by='reg'][,mean(V1)]

#+ Richness-ts-fig.cap, echo=FALSE
richness.ts.fig.cap <- "Time series of MSOM estimates of region richness. Each point is the posterior mean of regional richness in a year. Lines indicate long-term trends from fitted values of linear regression models predicting richness from time."
#' ###Figure 1
#+ Richness-ts-fig, fig.height=5, fig.width=3.5, fig.cap=richness.ts.fig.cap, ecal=TRUE, echo=TRUE
richness_ts()


#+ Richness-trend, include=TRUE, echo=TRUE
#' ###Richness Trend
# ---- long-term trend ----
#' Examine trends for naive richness:
load("pkgBuild/results/rich_naive_trend_kendall.RData")
rich_naive_trend_kendall[reg!="wcann",BH:=p.adjust(taup, method='BH')]
rich_naive_trend_kendall <- rich_naive_trend_kendall[reg!="wcann",list(reg=reg, estimate=tau, BH=BH, p.value=taup)]
print(rich_naive_trend_kendall)

#' Examine trends for MSOM richness:
load("pkgBuild/results/rich_trend_kendall.RData")
rich_trend_kendall[reg!="wcann",BH:=p.adjust(pvalue, method="BH")]
rich_trend_kendall <- rich_trend_kendall[reg!="wcann",list(reg=reg, estimate=tau, BH=BH, p.value=pvalue)]
print(rich_trend_kendall)


#+ 












