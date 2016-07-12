#' trawlDiversity: A package for calculating species richness from bottom trawl surveys
#' 
#' The package is associated with a study aimed at understanding long-term changes in species richness in marine ecosystems around the North American coastline. It employs multispecies occupancy models (MSOM) to estimate richness, as well as a few other methods. Many of the functions are aimed at visualizing data, and conducting supporting statistical analyses to understand the trends reported from MSOM posterior results (MSOM's are Bayesian models).
#' 
#' 
#' 
#' 
#' @importFrom grDevices adjustcolor col2rgb dev.off rainbow recordPlot rgb
#' @importFrom graphics abline barplot layout lines mtext pairs par plot points text
#' @importFrom stats coef cor model.frame model.matrix plogis pnorm predict rbinom residuals runif sd vcov
#' @importFrom utils head tail
#' @importFrom methods as
#' 
#' @docType package
#' @name trawlDiversity
NULL