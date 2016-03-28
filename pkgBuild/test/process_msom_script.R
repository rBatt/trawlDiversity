
library("rstan")
library("trawlDiversity")
library("rbLib")
library("R2jags")
library("maps")

setwd("~/Documents/School&Work/pinskyPost/trawl")

reg_file <- c(
	"msomStatic_norv_1yr_ebs_jags_6kIter_50nZ_start2016-03-22_r1.RData", # ebs 6k iter -- converged, mostly
	# "msomStatic_norv_1yr_ebs_jags_6kIter_89nZ_start2016-03-13_r1.RData", # ebs 6k iter
	# "msomStatic_norv_1yr_ebs_jags_start2016-03-10_16-44-44_r1.RData", # ebs 6k iter, converged
	
	"msomStatic_norv_1yr_ai_jags_12kIter_50nZ_start2016-03-23_r2.RData", # ai 12k iter -- converged
	# "msomStatic_norv_1yr_ai_jags_60kIter_61nZ_start2016-03-14_r2.RData", # ai, 60k iter, didn't converge
	# "msomStatic_norv_1yr_ai_jags_start2016-03-11_16-01-45_r2.RData", #ai 6k iter, didn't converge
	# "msomStatic_norv_1yr_ai_jags_start2016-03-06_19-10-49_r2.RData", # ai 30k iter, didn't converge
	
	"msomStatic_norv_1yr_goa_jags_12kIter_50nZ_start2016-03-23_r3.RData", # goa 12k iter --- converged
	# "msomStatic_norv_1yr_goa_jags_60kIter_62nZ_start2016-03-15_r3.RData", # goa 60k iter, not converged by improving
	# "msomStatic_norv_1yr_goa_jags_start2016-03-11_17-55-00_r3.RData", # goa 6k iter, didn't converge
	# "msomStatic_norv_1yr_goa_jags_start2016-03-07_06-21-58_r3.RData", # goa 30k iter, didn't converge
	
	"msomStatic_norv_1yr_wctri_jags_12kIter_50nZ_start2016-03-23_r4.RData", # wctri 12k iter
	# "msomStatic_norv_1yr_wctri_jags_30kIter_232nZ_start2016-03-16_r4.RData", # wctri 30k iter 200 n0s, didn't converge
	# "msomStatic_norv_1yr_wctri_jags_start2016-03-11_21-14-14_r4.RData", # wctri 6k iter, didn't converge (also too few 0 spp)
	# "msomStatic_norv_1yr_wctri_jags_start2016-03-08_10-53-18_r4.RData", # wctri 30k iter, didn't converge (better than most tho); also, needs more than 100 0 spp
	
	"msomStatic_norv_1yr_wcann_jags_12kIter_50nZ_start2016-03-24_r5.RData", # wcann 12k iter -- mostly converged, except Omega seems to bump up against 1, Omega and alpha1 are negatively correalted, and alpha2 and alpha3 are a bit negatively correlated (correlations for first year; other years don't look as bad)
	# "msomStatic_norv_1yr_wcann_jags_30kIter_216nZ_start2016-03-17_r5.RData", # wcann 30k iter, didn't converge
	# "msomStatic_norv_1yr_wcann_jags_start2016-03-11_23-18-31_r5.RData", # wcann 6k iter, didn't converge, too few 0 sp (50)
	# "msomStatic_norv_1yr_wcann_jags_start2016-03-07_13-25-01_r5.RData", # wcann 30k iter, didn't converge, too few 0 spp
	
	"msomStatic_norv_1yr_gmex_jags_12kIter_50nZ_start2016-03-24_r6.RData", # gmex 12k iter --- mixing, mostly but not quite stationary
	# "msomStatic_norv_1yr_gmex_jags_10kIter_213nZ_start2016-03-18_r6.RData", # gmex 10k iter, didn't really converge (not worst i've seen), too few 0 spp
	# "msomStatic_norv_1yr_gmex_jags_start2016-03-12_02-53-26_r6.RData", # gmex, 6k iter, didn't converge, too few 0 spp
	
	"msomStatic_norv_1yr_sa_jags_12kIter_50nZ_start2016-03-24_r7.RData", # sa 12k iter, some years didn't converge, specifically the 8th year (1997) looks really bad .... but most years look decent
	# "msomStatic_norv_1yr_sa_jags_10kIter_200nZ_start2016-03-18_r7.RData", # sa 10k iter, crap fit b/c trawlData type/error
	# "msomStatic_norv_1yr_sa_jags_start2016-03-12_04-26-53_r7.RData", # sa, 6k iter, almost converged, but too few n0 spp
	# "msomStatic_norv_1yr_sa_jags_start2016-03-06_16-15-08_r7.RData", # sa, 30k iter, very nearly converged, but too few n0 spp
	
	"msomStatic_norv_1yr_neus_jags_12kIter_50nZ_start2016-03-24_r8.RData", # neus, 10k iter -- pretty good convergence, some posterior correlation, e.g.,between Omega and alpha1
	# "msomStatic_norv_1yr_neus_jags_10kIter_230nZ_start2016-03-19_r8.RData", # neus, 10k iter, mostly converged, some years too few n0 spp
	# "msomStatic_norv_1yr_neus_jags_start2016-03-12_06-01-10_r8.RData", # neus, 6k iter, mostly converged, but too few n0 spp
	
	"msomStatic_norv_1yr_shelf_jags_6kIter_50nZ_start2016-03-25_r9.RData", # shelf 6k iter, some years really bad, like 3rd year
	# "msomStatic_norv_1yr_shelf_jags_40kIter_61nZ_start2016-03-19_r9.RData", # shelf 40k iter, mostly converged, but not perfect (oddly, Omega posterior correlated with alpha3)
	# "msomStatic_norv_1yr_shelf_jags_start2016-03-12_12-00-07_r9.RData", # shelf, 6k iter, very nearly converged
	# "msomStatic_norv_1yr_shelf_jags_start2016-03-06_14-03-15_r9.RData", # shelf, 30k iter, better than 6k, but not perfect
	
	"msomStatic_norv_1yr_newf_jags_6kIter_50nZ_start2016-03-25_r10.RData" # newf 6k iter
	# "msomStatic_norv_1yr_newf_jags_6kIter_53nZ_start2016-03-20_r10.RData" # newf 6k iter, good except 2008, when Omega and Beta had problems and were correlated in a chain, and is also the year w/ really high richness :p
	# "msomStatic_norv_1yr_newf_jags_start2016-03-12_13-30-24_r10.RData" # newf, 6k iter, really good mixing, mild lack of stationary
	
)

p <- list()
for(i in 1:length(reg_file)){
	load(paste0("trawlDiversity/pkgBuild/results/", reg_file[i]))
	
	reg_results_ind <- which(sapply(rm_out, function(x)!is.null(x)))
	stopifnot(length(reg_results_ind) == 1)
	
	rm_out <- rm_out[[reg_results_ind]]
	
	
	p[[i]] <- process_msomStatic(rm_out)
	
	rm(list="rm_out")
	gc()
	
}


save(p, file="trawlDiversity/pkgBuild/results/processedMsom/p.RData")




	