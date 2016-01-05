
# run on amphiprion: 
# nohup R CMD BATCH -cwd --no-save trawlDiversity/pkgBuild/test/msom_source_script.R &

# ========
# = Load =
# ========
library("trawlDiversity")
library("rstan")

Sys.time()
sessionInfo()

data_in_all <- trim_msom("ebs", gridSize=5, grid_stratum=FALSE, plot=FALSE)

# data_in_year <- data_in_all[year<=2000]

rm_out <- run_msom(
	reg = "ebs", 
	regX.a1 = data_in_all,
	params_out = c("params"), 
	language="Stan",
	model_type = "Static",
	test=FALSE, n0=50, iter=50, pre_save=F
)

Sys.time()
sessionInfo()

save.image(rbLib::renameNow(rm_out[[2]]["save_path"]), compress="xz")
