
# run on amphiprion: 
# nohup R CMD BATCH -cwd --no-save trawlDiversity/pkgBuild/test/msom_source_script.R &

# ========
# = Load =
# ========
library("trawlDiversity")

sessionInfo()

rm_out <- run_msom(
	reg = "ebs", 
	params_out = c("params"), 
	# custom_params = c("Z", "Psi"), 
	test=TRUE, n0=10, iter=200, pre_save=TRUE
)

save.image(renameNow(rm_out[[2]]["save_path"]), compress="xz")
