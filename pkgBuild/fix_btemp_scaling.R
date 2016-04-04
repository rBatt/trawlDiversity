


run_msom_get_btemp_sd <- function(reg = c("ai", "ebs", "gmex", "goa", "neus", "newf", "ngulf", "sa", "sgulf", "shelf", "wcann", "wctri"), regX.a1, params_out=c("params","params_main","params_random","params_latent","custom"), custom_params=NULL, model_type=c("Dynamic", "Static"), n0=50, chains=4, cores=parallel::detectCores()/2, iter, thin=max(1, floor((iter/2)/200)), language=c("JAGS", "Stan"), test=FALSE, test_sub=list(stratum=4, year=3, spp=10), seed=1337, pre_save=FALSE, save_dir=".", model_dir=file.path(system.file(package="trawlDiversity"), tolower(language)), compiled_model=NULL, ...){
	
	model_type <- match.arg(model_type)
	language <- match.arg(language)
	reg <- match.arg(reg)
	params_out <- match.arg(params_out, several.ok=TRUE)
	
	
	requireNamespace("parallel", quietly=TRUE)
	requireNamespace("trawlData", quietly=TRUE)
	requireNamespace("rbLib", quietly=TRUE)
	if(language=="Stan"){requireNamespace("rstan", quietly=TRUE)}
	if(language=="JAGS"){requireNamespace("R2jags", quietly=TRUE)}
	if(language=="JAGS"){library("R2jags")}
		
	if(missing(iter)){
		if(test){
			iter <- ifelse(language=="Stan", 50, 500)
		}else{
			iter <- ifelse(language=="Stan", 200, 5E3)
		}
	}

	
	
	# check files/ paths
	if(!file.exists(save_dir)){
		stop("save directory (", save_dir, ") does not exist")
	}
	if(!file.exists(model_dir) & ((is.null(compiled_model) & language=="Stan") | language=="JAGS")){
		stop("model directory (", model_dir, ") does not exist")
	}


	# ======================
	# = Subset Data to Use =
	# ======================
	if(missing(regX.a1)){
		regX.a1 <- trim_msom(reg, gridSize=1, plot=FALSE)
	}
	

	if(test){
		set.seed(seed)
		ind <- mpick(regX.a1, p=c(stratum=test_sub$stratum, year=test_sub$year), weight=TRUE, limit=60)
		logic <- expression(
			spp%in%spp[ind]
			& stratum%in%stratum[ind]
			& year%in%year[ind]
		)
		regX.a1 <- regX.a1[eval(logic)][pick(spp, test_sub$spp, w=FALSE)]
	}
	
	regX.a2 <- regX.a1[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, depth=depth, doy=yday(datetime))]


	# ==================
	# = Set Covariates =
	# ==================
	cov.vars_use <- c(bt="bt",bt2="bt2",yr="yr", yr2="yr2", doy="doy", doy2="doy2", depth="depth", depth2="depth2")


	# ======================================
	# = Aggregate and Transform Covariates =
	# ======================================
	# rename columns for shorthand
	setnames(regX.a2, c("btemp"), c("bt"))
	# setnames(regX.a2, c("stemp"), c("st"))
	# regX.a2[,yr:=scale(as.integer(year))] # fails with 1 year; also, mu and sd not weighted to unique years, so kinda weird
	regX.a2[,yr:=as.numeric(year)]
	regX.a2[,year:=as.character(year)]

	# aggregate and transform (^2) btemp stemp
	btemp.mu <- regX.a2[,mean(bt, na.rm=TRUE)]
	btemp.sd <- regX.a2[,sd(bt, na.rm=TRUE)]
	
	return(btemp.sd)
	
}


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
	
	"msomStatic_norv_1yr_neus_jags_12kIter_50nZ_start2016-03-24_r8.RData", # neus, 10k iter
	# "msomStatic_norv_1yr_neus_jags_10kIter_230nZ_start2016-03-19_r8.RData", # neus, 10k iter, mostly converged, some years too few n0 spp
	# "msomStatic_norv_1yr_neus_jags_start2016-03-12_06-01-10_r8.RData", # neus, 6k iter, mostly converged, but too few n0 spp
	
	"msomStatic_norv_1yr_shelf_jags_6kIter_50nZ_start2016-03-25_r9.RData", # shelf 6k iter
	# "msomStatic_norv_1yr_shelf_jags_40kIter_61nZ_start2016-03-19_r9.RData", # shelf 40k iter, mostly converged, but not perfect (oddly, Omega posterior correlated with alpha3)
	# "msomStatic_norv_1yr_shelf_jags_start2016-03-12_12-00-07_r9.RData", # shelf, 6k iter, very nearly converged
	# "msomStatic_norv_1yr_shelf_jags_start2016-03-06_14-03-15_r9.RData", # shelf, 30k iter, better than 6k, but not perfect
	
	"msomStatic_norv_1yr_newf_jags_6kIter_50nZ_start2016-03-25_r10.RData" # newf 6k iter
	# "msomStatic_norv_1yr_newf_jags_6kIter_53nZ_start2016-03-20_r10.RData" # newf 6k iter, good except 2008, when Omega and Beta had problems and were correlated in a chain, and is also the year w/ really high richness :p
	# "msomStatic_norv_1yr_newf_jags_start2016-03-12_13-30-24_r10.RData" # newf, 6k iter, really good mixing, mild lack of stationary
	
)






library("trawlDiversity")
library("raster")
library("sp")
library("rstan")
library("R2jags")

Sys.time()
sessionInfo()

regs <- c("ebs", "ai", "goa", "wctri", "wcann", "gmex", "sa", "neus", "shelf", "newf")


stan_folder <- file.path(system.file(package="trawlDiversity"), tolower("Stan"))
model_location <- file.path(stan_folder, "msomStatic_norv_1yr.stan")
compiled_stan_model <- stan_model(model_location)


reg_n0_pad <- c(
	"ebs" = 50,
	"ai" = 50,
	"goa" = 50,
	"wctri" = 50,
	"wcann" = 50,
	"gmex" = 50,
	"sa" = 50,
	"neus" = 50,
	"shelf" = 50,
	"newf" = 50
)

reg_iter <- c(
	"ebs" = 6E3,
	"ai" = 12E3,
	"goa" = 12E3,
	"wctri" = 12E3, 
	"wcann" = 12E3, 
	"gmex" = 12E3, 
	"sa" = 12E3, 
	"neus" = 12E3, 
	"shelf" = 6E3, 
	"newf" = 6E3
)


for(r in 1:length(regs)){
	
	rm_out <- vector("list", length(regs)) # yes, this reset the contents of the list. Saving all regions together is too big
	
	t_reg <- regs[r]
	
	data_in_all0 <- trim_msom(t_reg, gridSize=0.5, depthStratum=100, tolFraction=0.15, grid_stratum=TRUE, plot=FALSE)
	data_in_all <- data_in_all0
	setkey(data_in_all, year, stratum, haulid, spp)
	
	u_yrs <- data_in_all[,unique(year)]
	n_spp <- data_in_all[,list(n_spp=lu(spp)), by="year"]
	
	S <- data_in_all[,lu(spp)]
	annual_n0 <- (S + reg_n0_pad[regs[r]]) - n_spp[,n_spp]
	
	rm_out[[r]] <- vector("list", length(u_yrs)) 
	
	# ---- load old results ----
	load(paste0("trawlDiversity/pkgBuild/results/", reg_file[r]))
	reg_results_ind <- which(sapply(rm_out, function(x)!is.null(x)))
	stopifnot(length(reg_results_ind) == 1)
	
	for(i in 1:length(u_yrs)){
		t_data <- data_in_all[year==u_yrs[i]]
	
		msg_reg <- toupper(t_data[,unique(reg)])
		msg_yr_id <- paste0("Year = ",u_yrs[i])
		msg_yr_cnt <-paste0("(", i, " of ", length(u_yrs), ")")
		msg_progress <- paste(msg_reg, msg_yr_id, msg_yr_cnt)
		cat(paste("\n\n\n", msg_progress, "\n"))
		print(Sys.time())
		
		# ---- replace bad btemp scaling ----
		rm_out[[r]][[i]]$inputData$scaling["btemp.sd"] <- run_msom_get_btemp_sd(
			reg = t_reg,
			regX.a1 = t_data,
			params_out = c("params"),
			language="JAGS", 
			model_type = "Static", 
			compiled_model = compiled_stan_model,
			cores = 4, chains = 4,
			test=FALSE, n0=annual_n0[i], iter=reg_iter[regs[r]], pre_save=FALSE, save_warmup=FALSE
		)
		
		print(Sys.time())
		cat("\n\n")
		
	}
	
	
	save(rm_out, file=paste0("trawlDiversity/pkgBuild/results/", reg_file[r]), compress="xz")
	
	cat(paste0("\n\n\n\n\n",  paste(rep("=", 50), collapse=""), "\n\n\n"))
	
}






