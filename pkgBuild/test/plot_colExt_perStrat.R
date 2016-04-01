plot_colExt_perStrat <- function(prn, Figures){
	unpack_p(prn)
	
	if(missing(Figures)){
		Figures <- list()
	}
	
	r_lon <- colonization$n_spp_col_weighted_tot[,range(lon)]
	r_lat <- colonization$n_spp_col_weighted_tot[,range(lat)]
	
	lay_logic <- diff(r_lon) > diff(r_lat)
	if(lay_logic){
		fig7_mfr <- c(2,1*3)
		fig7_h <- 2.5#*3
		fig7_w <- (diff(r_lon) * fig7_h / diff(r_lat)) * (fig7_mfr[2]/fig7_mfr[1])
	}else{
		# fig7_mfr <- c(1*3,2)
# 		fig7_w <- 2.5#*3
# 		fig7_h <- (diff(r_lat) * fig7_w / diff(r_lon)) / (fig7_mfr[2]/fig7_mfr[1])
		fig7_mfr <- c(2,1*3)
		fig7_h <- 7#*3
		fig7_w <- (diff(r_lon) * fig7_h / diff(r_lat)) * (fig7_mfr[2]/fig7_mfr[1])
		
		lay_logic <- !lay_logic
	}

	fig7_name <- paste0("colonizations_per_stratum_", reg, ".png")
	fig7_dim <- c(fig7_w, fig7_h)
	dev.new(width=fig7_w, height=fig7_h)
	if(lay_logic){
		par(mfcol=fig7_mfr)
	}else{
		par(mfrow=fig7_mfr)
	}
	par(mar=c(1.25,1.25,0.1,0.1), oma=c(0.1,1,1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=8, cex=1)

	# site-specific colonizations from data
	colonization$n_spp_col_weighted_tot[,plot_space(lon, lat, n_spp_col_weighted, TRUE, pch=19, ylab="", xlab="")]
	map(add=TRUE, fill=TRUE, col="white")
	mtext("Colonizations (C)", side=ifelse(lay_logic, 3, 1), line=ifelse(lay_logic, 0.25, 2))
	# smoothed map for convex hull of site observations
	colonization$n_spp_col_weighted_tot[,plot_space(lon, lat, n_spp_col_weighted, pch=19, ylab="", xlab="")]
	map(add=TRUE, fill=TRUE, col="white")
	
	
	# site-specific extinctions from data
	colonization$n_spp_ext_weighted_tot[,plot_space(lon, lat, n_spp_ext_weighted, TRUE, pch=19, ylab="", xlab="")]
	map(add=TRUE, fill=TRUE, col="white")
	mtext("Extinctions (E)", side=ifelse(lay_logic, 3, 1), line=ifelse(lay_logic, 0.25, 2))
	# smoothed map for convex hull of site observations
	colonization$n_spp_ext_weighted_tot[,plot_space(lon, lat, n_spp_ext_weighted, pch=19, ylab="", xlab="")]
	map(add=TRUE, fill=TRUE, col="white")
	
	# site-specific extinctions from data
	cre <- merge(colonization$n_spp_col_weighted_tot, colonization$n_spp_ext_weighted_tot, by=c("stratum","lon","lat","depth"), all=TRUE)
	# cre <- data.table(cre[,list(stratum, n_spp_ext_weighted, n_spp_col_weighted)], cre[,strat2lld(stratum)])
# 	cre[is.na(n_spp_ext_weighted), n_spp_ext_weighted:=0]
# 	cre[is.na(n_spp_col_weighted), n_spp_col_weighted:=0]
	cre[,rel_col_ext:=(n_spp_col_weighted - n_spp_ext_weighted) ]
	cre[,plot_space(lon, lat, rel_col_ext, TRUE, pch=19, ylab="", xlab="")]
	map(add=TRUE, fill=TRUE, col="white")
	mtext("C - E", side=ifelse(lay_logic, 3, 1), line=ifelse(lay_logic, 0.25, 2))
	# smoothed map for convex hull of site observations
	cre[,plot_space(lon, lat, rel_col_ext, pch=19, ylab="", xlab="")]
	map(add=TRUE, fill=TRUE, col="white")
	
	

	Figures[[reg]][['Figure7']][["figure"]] <- recordPlot()
	Figures[[reg]][['Figure7']][["name"]] <- fig7_name
	Figures[[reg]][['Figure7']][["dim"]] <- fig7_dim

	return(Figures)
	
	# ---- Figure 8 ----
	# ---- Number of extinctions per stratum ----
	# fig8_name <- paste0("extinctions_per_stratum_", reg, ".png")
# 	fig8_dim <- c(fig7_w, fig7_h)
# 	dev.new(width=fig7_w, height=fig7_h)
# 	par(mfrow=fig7_mfr, mar=c(1.5,1.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=8, cex=1)

	# # site-specific extinctions from data
# 	colonization$n_spp_ext_weighted_tot[,plot_space(lon, lat, n_spp_ext_weighted, TRUE, pch=19)]
# 	map(add=TRUE, fill=TRUE, col="white")
# 	# smoothed map for convex hull of site observations
# 	colonization$n_spp_ext_weighted_tot[,plot_space(lon, lat, n_spp_ext_weighted, pch=19)]
# 	map(add=TRUE, fill=TRUE, col="white")

	# Figures[[reg]][['Figure8']][["figure"]] <- recordPlot()
	# Figures[[reg]][['Figure8']][["name"]] <- fig8_name
	# Figures[[reg]][['Figure8']][["dim"]] <- fig8_dim
	# dev.off()
	
	
	# ---- Figure 8.5 ----
	# ---- Colonizations relative to extinctions per stratum ----
	# fig8.5_name <- paste0("ColRelExt_per_stratum_", reg, ".png")
	# fig8.5_dim <- c(fig7_w, fig7_h)
	# dev.new(width=fig7_w, height=fig7_h)
	# par(mfrow=fig7_mfr, mar=c(1.5,1.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=8, cex=1)
	#
	# strat2lld <- function(x){
	# 	s <- strsplit(x, split=" ")
	# 	lon <- sapply(s, function(x)x[1])
	# 	lat <- sapply(s, function(x)x[2])
	# 	depth_interval <- sapply(s, function(x)x[3])
	# 	data.table(lon=as.numeric(lon), lat=as.numeric(lat), depth_interval=as.numeric(depth_interval))
	# } 

	# # site-specific extinctions from data
# 	cre <- merge(colonization$n_spp_col_weighted_tot, colonization$n_spp_ext_weighted_tot, by=c("stratum","lon","lat","depth"), all=TRUE)
# 	# cre <- data.table(cre[,list(stratum, n_spp_ext_weighted, n_spp_col_weighted)], cre[,strat2lld(stratum)])
# # 	cre[is.na(n_spp_ext_weighted), n_spp_ext_weighted:=0]
# # 	cre[is.na(n_spp_col_weighted), n_spp_col_weighted:=0]
# 	cre[,rel_col_ext:=(n_spp_col_weighted - n_spp_ext_weighted) ]
# 	cre[,plot_space(lon, lat, rel_col_ext, TRUE, pch=19)]
# 	map(add=TRUE, fill=TRUE, col="white")
# 	# smoothed map for convex hull of site observations
# 	cre[,plot_space(lon, lat, rel_col_ext, pch=19)]
# 	map(add=TRUE, fill=TRUE, col="white")
}