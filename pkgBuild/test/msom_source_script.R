


# run on amphiprion: 
# nohup R CMD BATCH -cwd --no-save trawlDiversity/pkgBuild/test/msom_source_script.R &

# ========
# = Load =
# ========
library("devtools")
library("trawlData")
library("trawlDiversity")
library("rbLib")
library("raster")
library("rstan")
library("R2jags")
library("parallel")



# ---- alternative approach to doing stan model, maybe ----
# probably requires installing StanHeaders and rstan from source on development GitHub
# msom_model <- stan_model(file=model_path, auto_write=TRUE)
#
# do_stan <- function(x, object, data, iter){
# 	sampling(
# 			object=object,
# 			data=inputData,
# 			control=stan_control,
# 			pars= params_out,
# 			chains=1, iter=iter, seed=1337, cores=1, verbose=F, refresh=1, chain_id=x
# 	)
# }
#
#
# t_out <- sampling(
# 		object=msom_model,
# 		data=inputData,
# 		control=stan_control,
# 		pars= params_out,
# 		chains=2, iter=4, seed=1337, cores=1, verbose=F, refresh=1
# )
#
# out <- do_stan(1, msom_model, data=inputData, iter=iter)
# # read1 <- read_stan_csv("ebs_msom_samples1")
#
# out <- mclapply(1:chains, do_stan, data=inputData, object=msom_model, iter=iter, mc.cores=cores, mc.preschedule=FALSE)

sessionInfo()




save.image(renameNow(save_path), compress="xz")
