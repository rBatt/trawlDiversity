
# run on amphiprion: 
# nohup R CMD BATCH -cwd --no-save trawlDiversity/pkgBuild/test/msom_source_script.R &

# ========
# = Load =
# ========
library("trawlDiversity")

Sys.time()
sessionInfo()

rm_out <- run_msom(
	reg = "ebs", 
	params_out = c("params"), 
	test=FALSE, n0=50, iter=100, pre_save=TRUE
)

Sys.time()
sessionInfo()

save.image(rbLib::renameNow(rm_out[[2]]["save_path"]), compress="xz")
