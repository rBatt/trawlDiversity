
# run on amphiprion: 
# nohup R CMD BATCH -cwd --no-save trawlDiversity/pkgBuild/test/msom_source_script.R &

# nohup R CMD BATCH -cwd --no-save trawlDiversity/pkgBuild/test/msom_source_script.R msom_source_script_AllRegs_annual_jags_noCulling_r6-7.Rout &


# amphiprion process id's
# ebs: 2041
# ai & goa: 2138
# wctri & wcann: 2307


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


# ===========
# = Options =
# ===========
lang_used <- c("JAGS", "Stan")[1]

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
	"ebs" = 12E3,
	"ai" = 12E3,
	"goa" = 12E3,
	"wctri" = 12E3, 
	"wcann" = 12E3, 
	"gmex" = 12E3, 
	"sa" = 12E3, 
	"neus" = 12E3, 
	"shelf" = 12E3, 
	"newf" = 12E3
)


# =========
# = Setup =
# =========
regs <- c("ebs", "ai", "goa", "wctri", "wcann", "gmex", "sa", "neus", "shelf", "newf")


if(lang_used=="Stan"){
	stan_folder <- file.path(system.file(package="trawlDiversity"), tolower("Stan"))
	model_location <- file.path(stan_folder, "msomStatic_norv_1yr.stan")
	compiled_stan_model <- stan_model(model_location)
}else{
	compiled_stan_model <- NULL
}


# =====================
# = Loop and Run MSOM =
# =====================
# for(r in 1:length(regs)){
# for(r in 1:1){ # ebs
# for(r in 2:3){ # ai & goa
# for(r in 4:5){ # wctri & wcann
for(r in 6:7){ # gmex & sa
# for(r in 8:8){ # neus
# for(r in 9:10){ # shelf and newf

	rm_out <- vector("list", length(regs)) # yes, this reset the contents of the list. Saving all regions together is too big
	
	t_reg <- regs[r]
	data_in_all <- data_all[reg==t_reg]
	
	u_yrs <- data_in_all[,unique(year)]
	n_spp <- data_in_all[,list(n_spp=lu(spp)), by="year"]
	
	S <- data_in_all[,lu(spp)]
	annual_n0 <- (S + reg_n0_pad[regs[r]]) - n_spp[,n_spp]
	
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
			language=lang_used, 
			model_type = "Static", 
			compiled_model = compiled_stan_model,
			cores = 4, chains = 4,
			test=FALSE, n0=annual_n0[i], iter=reg_iter[regs[r]], pre_save=FALSE, save_warmup=FALSE
		), error=function(cond){message(paste("**Run failed for:", msg_progress));NA})
		
		rm_out[[r]][[i]]$info['year'] <- u_yrs[i]
		
		print(Sys.time())
		cat("\n\n")
		
	}
	
	append_r <- paste0("_r", r, ".RData")
	save_name <- gsub("\\.RData", append_r, rm_out[[r]][[1]][[3]]["save_path"])
	save_name <- gsub("^\\./", "./trawlDiversity/pkgBuild/results/", save_name)
	save_name <- gsub("[0-9]+(?=nZ)",reg_n0_pad[regs[r]],save_name, perl=TRUE)
	save(rm_out, file=save_name, compress="xz")
	
	cat(paste0("\n\n\n\n\n",  paste(rep("=", 50), collapse=""), "\n\n\n"))
	
}


Sys.time()
sessionInfo()