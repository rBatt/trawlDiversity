
# nohup R CMD BATCH -cwd --no-save pkgBuild/test/process_msom_script.R process_msom_script.Rout &

library("rstan")
library("trawlDiversity")
library("rbLib")
library("R2jags")
library("maps")

setwd("~/Documents/School&Work/pinskyPost/trawlDiversity")

reg_file <- c(
		"msomStatic_norv_1yr_ebs_jags_12kIter_50nZ_start2016-05-07_r1.RData",
	# "msomStatic_norv_1yr_ebs_jags_12kIter_50nZ_start2016-04-04_r1.RData", # ebs 12k iter, converged
	# "msomStatic_norv_1yr_ebs_jags_6kIter_50nZ_start2016-03-22_r1.RData", # ebs 6k iter -- converged, mostly
	# "msomStatic_norv_1yr_ebs_jags_6kIter_89nZ_start2016-03-13_r1.RData", # ebs 6k iter
	# "msomStatic_norv_1yr_ebs_jags_start2016-03-10_16-44-44_r1.RData", # ebs 6k iter, converged
	
		"msomStatic_norv_1yr_ai_jags_12kIter_50nZ_start2017-04-28_r2.RData", # leopard skate removed
	# "msomStatic_norv_1yr_ai_jags_12kIter_50nZ_start2016-05-07_r2.RData",
	# "msomStatic_norv_1yr_ai_jags_12kIter_50nZ_start2016-04-04_r2.RData", # ai 12k iter, 4th yr didn't converge
	# "msomStatic_norv_1yr_ai_jags_12kIter_50nZ_start2016-03-23_r2.RData", # ai 12k iter -- converged
	# "msomStatic_norv_1yr_ai_jags_60kIter_61nZ_start2016-03-14_r2.RData", # ai, 60k iter, didn't converge
	# "msomStatic_norv_1yr_ai_jags_start2016-03-11_16-01-45_r2.RData", #ai 6k iter, didn't converge
	# "msomStatic_norv_1yr_ai_jags_start2016-03-06_19-10-49_r2.RData", # ai 30k iter, didn't converge
	
		"msomStatic_norv_1yr_goa_jags_12kIter_50nZ_start2016-05-07_r3.RData",
	# "msomStatic_norv_1yr_goa_jags_12kIter_50nZ_start2016-04-05_r3.RData", # goa 12k iter, converged
	# "msomStatic_norv_1yr_goa_jags_12kIter_50nZ_start2016-03-23_r3.RData", # goa 12k iter --- converged
	# "msomStatic_norv_1yr_goa_jags_60kIter_62nZ_start2016-03-15_r3.RData", # goa 60k iter, not converged by improving
	# "msomStatic_norv_1yr_goa_jags_start2016-03-11_17-55-00_r3.RData", # goa 6k iter, didn't converge
	# "msomStatic_norv_1yr_goa_jags_start2016-03-07_06-21-58_r3.RData", # goa 30k iter, didn't converge
	
		"msomStatic_norv_1yr_wctri_jags_12kIter_50nZ_start2016-05-07_r4.RData",
	# "msomStatic_norv_1yr_wctri_jags_12kIter_50nZ_start2016-04-04_r4.RData", # wctri 12k iter, good convergence
	# "msomStatic_norv_1yr_wctri_jags_12kIter_50nZ_start2016-03-23_r4.RData", # wctri 12k iter
	# "msomStatic_norv_1yr_wctri_jags_30kIter_232nZ_start2016-03-16_r4.RData", # wctri 30k iter 200 n0s, didn't converge
	# "msomStatic_norv_1yr_wctri_jags_start2016-03-11_21-14-14_r4.RData", # wctri 6k iter, didn't converge (also too few 0 spp)
	# "msomStatic_norv_1yr_wctri_jags_start2016-03-08_10-53-18_r4.RData", # wctri 30k iter, didn't converge (better than most tho); also, needs more than 100 0 spp
	
		"msomStatic_norv_1yr_wcann_jags_12kIter_50nZ_start2016-05-07_r5.RData",
	# "msomStatic_norv_1yr_wcann_jags_12kIter_50nZ_start2016-04-05_r5.RData", # wcann 12k iter, good convergence, maybe a little Omega close to 1 in some years
	# "msomStatic_norv_1yr_wcann_jags_12kIter_50nZ_start2016-03-24_r5.RData", # wcann 12k iter -- mostly converged, except Omega seems to bump up against 1, Omega and alpha1 are negatively correalted, and alpha2 and alpha3 are a bit negatively correlated (correlations for first year; other years don't look as bad)
	# "msomStatic_norv_1yr_wcann_jags_30kIter_216nZ_start2016-03-17_r5.RData", # wcann 30k iter, didn't converge
	# "msomStatic_norv_1yr_wcann_jags_start2016-03-11_23-18-31_r5.RData", # wcann 6k iter, didn't converge, too few 0 sp (50)
	# "msomStatic_norv_1yr_wcann_jags_start2016-03-07_13-25-01_r5.RData", # wcann 30k iter, didn't converge, too few 0 spp
	
		"msomStatic_norv_1yr_gmex_jags_12kIter_50nZ_start2016-05-07_r6.RData",
	# "msomStatic_norv_1yr_gmex_jags_12kIter_50nZ_start2016-04-04_r6.RData", # gmex 12k iter, converged
	# "msomStatic_norv_1yr_gmex_jags_12kIter_50nZ_start2016-03-24_r6.RData", # gmex 12k iter --- mixing, mostly but not quite stationary
	# "msomStatic_norv_1yr_gmex_jags_10kIter_213nZ_start2016-03-18_r6.RData", # gmex 10k iter, didn't really converge (not worst i've seen), too few 0 spp
	# "msomStatic_norv_1yr_gmex_jags_start2016-03-12_02-53-26_r6.RData", # gmex, 6k iter, didn't converge, too few 0 spp
	
		"msomStatic_norv_1yr_sa_jags_12kIter_50nZ_start2016-05-07_r7.RData",
	# "msomStatic_norv_1yr_sa_jags_12kIter_50nZ_start2016-04-05_r7.RData", # sa 12k iter, some years didn't converge: 8th, and kinda 20th; the 8th year looks a little bit better than it did before
	# "msomStatic_norv_1yr_sa_jags_12kIter_50nZ_start2016-03-24_r7.RData", # sa 12k iter, some years didn't converge, specifically the 8th year (1997) looks really bad .... but most years look decent
	# "msomStatic_norv_1yr_sa_jags_start2016-03-12_04-26-53_r7.RData", # sa, 6k iter, almost converged, but too few n0 spp
	# "msomStatic_norv_1yr_sa_jags_start2016-03-06_16-15-08_r7.RData", # sa, 30k iter, very nearly converged, but too few n0 spp
	
		"msomStatic_norv_1yr_neus_jags_12kIter_50nZ_start2016-05-07_r8.RData",
	# "msomStatic_norv_1yr_neus_jags_12kIter_50nZ_start2016-04-04_r8.RData", # neus 12k iter, converged; yrs 10-12, 16, and 26 Omega almost maxes out; some posterior correlation, modest nonstationarity
	# "msomStatic_norv_1yr_neus_jags_12kIter_50nZ_start2016-03-24_r8.RData", # neus, 10k iter -- pretty good convergence, some posterior correlation, e.g.,between Omega and alpha1
	# "msomStatic_norv_1yr_neus_jags_10kIter_230nZ_start2016-03-19_r8.RData", # neus, 10k iter, mostly converged, some years too few n0 spp
	# "msomStatic_norv_1yr_neus_jags_start2016-03-12_06-01-10_r8.RData", # neus, 6k iter, mostly converged, but too few n0 spp
	
		"msomStatic_norv_1yr_shelf_jags_12kIter_50nZ_start2016-05-07_r9.RData",
	# "msomStatic_norv_1yr_shelf_jags_12kIter_50nZ_start2016-04-04_r9.RData", # shelf 12k iter, convergence way better than before, but a bit spiky in a lot of places.
	# "msomStatic_norv_1yr_shelf_jags_6kIter_50nZ_start2016-03-25_r9.RData", # shelf 6k iter, some years really bad, like 3rd year
	# "msomStatic_norv_1yr_shelf_jags_40kIter_61nZ_start2016-03-19_r9.RData", # shelf 40k iter, mostly converged, but not perfect (oddly, Omega posterior correlated with alpha3)
	# "msomStatic_norv_1yr_shelf_jags_start2016-03-12_12-00-07_r9.RData", # shelf, 6k iter, very nearly converged
	# "msomStatic_norv_1yr_shelf_jags_start2016-03-06_14-03-15_r9.RData", # shelf, 30k iter, better than 6k, but not perfect
	
		"msomStatic_norv_1yr_newf_jags_12kIter_50nZ_start2016-05-07_r10.RData"
	# "msomStatic_norv_1yr_newf_jags_12kIter_50nZ_start2016-04-05_r10.RData" # newf 12k iter
	# "msomStatic_norv_1yr_newf_jags_6kIter_50nZ_start2016-03-25_r10.RData" # newf 6k iter
	# "msomStatic_norv_1yr_newf_jags_6kIter_53nZ_start2016-03-20_r10.RData" # newf 6k iter, good except 2008, when Omega and Beta had problems and were correlated in a chain, and is also the year w/ really high richness :p
	# "msomStatic_norv_1yr_newf_jags_start2016-03-12_13-30-24_r10.RData" # newf, 6k iter, really good mixing, mild lack of stationary
	
)

p <- list()
for(i in 1:length(reg_file)){
	cat(i)
	
	load(paste0("pkgBuild/results/", reg_file[i]))
	
	reg_results_ind <- which(sapply(rm_out, function(x)!is.null(x)))
	stopifnot(length(reg_results_ind) == 1)
	
	rm_out <- rm_out[[reg_results_ind]]
	
	#' #TODO: Remove if I ever re-run MSOM for all regions; I added a check on inputData to make sure it doesn't include species that aren't in data_all .... 
	# sort(unique(unlist(sapply(rm_out, function(x){(dimnames(x$inputData$X)$spp)[!grepl("Unknown", dimnames(x$inputData$X)$spp)]})))) # the unique species from all years (ignoring 'unknowns')	
	fix2names <- function(x){
		nms <- dimnames(x$inputData$X)$spp
		fixInd1 <- nms=="Reinhardtius stomias"
		fixed1 <- "Atheresthes stomias"
		fixInd2 <- nms=="Cross papposus"
		fixed2 <- "Crossaster papposus"
		dimnames(x$inputData$X)$spp[fixInd1] <- fixed1
		dimnames(x$inputData$X)$spp[fixInd2] <- fixed2
		return(x)
	}
	for(fi in 1:length(rm_out)){
		rm_out[[fi]] <- fix2names(rm_out[[fi]])
	}
	
	p[[i]] <- process_msomStatic(rm_out)
	
	t_reg <- p[[i]]$processed[,una(reg)]
	t_X <- data_all[reg==t_reg]
	msom_yrs <- p[[i]]$processed[,sort(una(year))]
	p_obs <- process_obsRich(X=t_X, msom_yrs=msom_yrs)
	p[[i]] <- c(p[[i]],p_obs)
	
	t_proc <- p[[i]]$processed
	t_proc <- merge(t_proc, p_obs$colonization$n_cep, by="year", all=TRUE)
	t_proc <- merge(t_proc, p_obs$bt[,list(bt_ann=mean(bt,na.rm=TRUE)), by="year"], by="year", all=TRUE)
	t_proc <- merge(t_proc, p[[i]]$naive_rich, by="year", all=TRUE)
	t_proc <- merge(t_proc, p[[i]]$local_rich_obs, by="year", all=TRUE)
	t_proc <- merge(t_proc, p[[i]]$local_rich_samp, by="year", all=TRUE)
	p[[i]]$processed <- t_proc
	
	rm(list="rm_out")
	gc()
	
}


save(p, file="pkgBuild/results/processedMsom/p.RData")
