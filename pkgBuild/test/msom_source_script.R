
# run on amphiprion: 
# nohup R CMD BATCH -cwd --no-save trawlDiversity/pkgBuild/test/msom_source_script.R &

# nohup R CMD BATCH -cwd --no-save trawlDiversity/pkgBuild/test/msom_source_script.R msom_source_script_AllRegs_annual_r9-10.Rout &


# ========
# = Load =
# ========
library("trawlDiversity")
library("raster")
library("sp")
library("rstan")
library("R2jags")

Sys.time()
sessionInfo()

n0_pad <- 100
regs <- c("ebs", "ai", "goa", "wctri", "wcann", "gmex", "sa", "neus", "shelf", "newf")


stan_folder <- file.path(system.file(package="trawlDiversity"), tolower("Stan"))
model_location <- file.path(stan_folder, "msomStatic_norv.stan")
compiled_stan_model <- stan_model(model_location)

# for(r in length(regs):1){
# for(r in 1:2){ # ebs and ai
# for(r in 3:4){ # goa and wctri
# for(r in 5:6){ # wcann and gmex
# for(r in 7:8){ # sa and neus
for(r in 9:10){ # shelf and newf
	
	rm_out <- vector("list", length(regs)) # yes, this reset the contents of the list. Saving all regions together is too big
	
	t_reg <- regs[r]
	
	data_in_all0 <- trim_msom(t_reg, gridSize=0.5, depthStratum=100, tolFraction=0.15, grid_stratum=TRUE, plot=FALSE)

	data_in_all <- data_in_all0
	setkey(data_in_all, year, stratum, haulid, spp)
	u_yrs <- data_in_all[,unique(year)]
	n_spp <- data_in_all[,list(n_spp=lu(spp)), by="year"]
	
	S <- data_in_all[,lu(spp)]
	annual_n0 <- (S + n0_pad) - n_spp[,n_spp]
	
	rm_out[[r]] <- vector("list", length(u_yrs))
	
	for(i in 1:length(u_yrs)){
		t_data <- data_in_all[year==u_yrs[i]]
	
		msg_reg <- toupper(t_data[,unique(reg)])
		msg_yr_id <- paste0("Year = ",u_yrs[i])
		msg_yr_cnt <-paste0("(", i, " of ", length(u_yrs), ")")
		msg_progress <- paste(msg_reg, msg_yr_id, msg_yr_cnt)
		cat(paste("\n\n\n", msg_progress, "\n"))
		print(Sys.time())
	
		rm_out[[r]][[i]] <- tryCatch(run_msom(
			reg = t_reg,
			regX.a1 = t_data,
			params_out = c("params"),
			language="Stan",
			model_type = "Static", 
			compiled_model = compiled_stan_model,
			cores = 4, chains = 4,
			test=FALSE, n0=annual_n0[i], iter=100, pre_save=FALSE, save_warmup=FALSE
		), error=function(cond){message(paste("**Run failed for:", msg_progress));NA})
		
		cat("\n\n")
		
	}
	
	append_r <- paste0("_r", r, ".RData")
	save_name <- gsub("\\.RData", append_r, rm_out[[r]][[i]][[3]]["save_path"])
	save_name <- gsub("^\\./", "./trawlDiversity/pkgBuild/results/", save_name)
	save(rm_out, file=save_name, compress="xz")
	
	cat(paste0("\n\n\n\n\n",  paste(rep("=", 50), collapse=""), "\n\n\n"))
	
}


Sys.time()
sessionInfo()