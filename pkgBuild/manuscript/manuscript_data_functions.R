library('maps')
library('raster')
library('spatstat')
library("spdep")
library('rbLib')
library('trawlDiversity')
library("data.table")

# ========================================================
# = Organize Community Level Richness-Range-Density Data =
# ========================================================
make_range_reg <- function(densName=c("propTow_occ_avg"), sizeName=c("range_size_samp_avg_ltAvg","range_size_mu_avg_ltAvg","propStrata_avg_ltAvg")){
	densName <- match.arg(densName)
	sizeName <- match.arg(sizeName)
	comm_master[,list(reg, year, rich=reg_rich, density=get(densName), size=get(sizeName))]
}


# ==============================================================
# = Organize Species Data on Range, Density, and Time-to-Event =
# ==============================================================
# Does not include years when a species is not present!
# This is to avoid the inevitable size=0 when time=0
# The data is handly for regressions of changes in range size with time to extinction/ time after colonization
make_rangeTime <- function(densName=c("propTow_occ"), sizeName=c("range_size_samp","range_size_mu","propStrata")){
	rangeTimeDT <-  spp_master[!is.na(stretch_type) & propStrata!=0]
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


# =====================================
# = List Trim MSOM Settings On Demand =
# =====================================
# MSOM Data-Trimming Settings
# Settings for trimming data before use in MSOM
# @param type character, indicating what setting to return. Two forms of the 'year' setting available: one for the logic of subsetting, the other convenient for plotting timing of cutoff. Depth is the grid size for depth-based stratum delineation. Tolerance is the proportion of years a stratum can be unsampled and still remain part of the data set.
# 
# @return named vector of length 9 or 10 (some include west coast annual, others don't)
trim_msom_settings <- function(type=c("depth","tolerance","years_logic","years_cutoff")){

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
	
	type <- match.arg(type)
	out <- switch(type,
		depth = reg_depthStratum,
		tolerance = reg_tolFraction,
		years_logic = yr_subs,
		years_cutoff = yr_ablin
	)
	
	return(out)
}


# ========================================================
# = Subsetting to Neither, Present, and Providing Colors =
# ========================================================
neitherPres_subColor <- function(){
	bquote({
		rE_logic1 <- bquote(reg==rr & ce_categ!="neither" & present==1)
		rE_logic2 <- bquote(reg==rr & ce_categ!="neither")
		rE_logic3 <- bquote(reg==rr & ce_categ=="neither")
		colQ <- bquote(adjustcolor(c("pre_ext"="red","post_col"="blue")[stretch_type],0.5))
		colQ2 <- bquote(adjustcolor(c("pre_ext"="red","post_col"="blue")[stretch_type],1))
	})
}

