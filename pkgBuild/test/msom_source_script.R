
# run on amphiprion: 
# nohup R CMD BATCH -cwd --no-save trawlDiversity/pkgBuild/test/msom_source_script.R &

# nohup R CMD BATCH -cwd --no-save trawlDiversity/pkgBuild/test/msom_source_script.R msom_source_script_AllRegs_annual.Rout &


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
# n0_pad <- 50
regs <- c("ebs", "ai", "goa", "wctri", "wcann", "gmex", "sa", "neus", "shelf", "newf")
rm_out <- vector("list", length(regs))

for(r in 1:length(regs)){
	
	t_reg <- regs[r]
	
	data_in_all0 <- trim_msom(t_reg, gridSize=1, grid_stratum=TRUE, plot=FALSE)

	data_in_all <- data_in_all0
	setkey(data_in_all, year, stratum, haulid, spp)
	u_yrs <- data_in_all[,unique(year)]
	n_spp <- data_in_all[,list(n_spp=lu(spp)), by="year"]
	# max_n_spp <- n_spp[,max(n_spp)]
	S <- data_in_all[,lu(spp)]
	annual_n0 <- (S + n0_pad) - n_spp[,n_spp]
	# annual_n0 <- (max_n_spp + n0_pad) - n_spp[,n_spp]

	
	rm_out[[r]] <- vector("list", length(u_yrs))
	
	for(i in 1:length(u_yrs)){
		t_data <- data_in_all[year==u_yrs[i]]
	
		cat(paste0("\n\n\n", toupper(t_data[,unique(reg)]), " Year = ",u_yrs[i], " (", i, " of ", length(u_yrs), ")\n"))
		print(Sys.time())
	
		rm_out[[r]][[i]] <- run_msom(
			reg = t_reg,
			regX.a1 = t_data,
			params_out = c("params"),
			language="Stan",
			model_type = "Static",
			cores = 4, chains=4,
			test=FALSE, n0=annual_n0[i], iter=60, pre_save=FALSE, save_warmup=FALSE
		)
		
		cat("\n\n")
	
	}
	
	append_r <- paste0("_r1-", r, ".RData")
	save_name <- gsub("\\.RData", append_r, rm_out[[r]][[i]][[3]]["save_path"])
	save_name <- gsub("^\\./", "./trawlDiversity/pkgBuild/results/", save_name)
	save(rm_out, file=save_name, compress="xz")
	
	cat(paste0("\n\n\n\n\n",  paste(rep("=", 50), collapse=""), "\n\n\n"))
	
}





Sys.time()
sessionInfo()

# save.image(rbLib::renameNow(rm_out[[1]][[3]]["save_path"]), compress="xz")
