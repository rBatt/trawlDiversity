library('maps')
library('raster')
library('spatstat')
library("spdep")
library('rbLib')
library('trawlDiversity')
library("data.table")

figure_setup <- function(){
	bquote({
		regs <- comm_master[,una(reg)] #sapply(p, function(x)x$processed[,una(reg)])
		pretty_reg <- c("ebs"="E. Bering Sea", "ai"="Aleutian Islands", "goa"="Gulf of Alaska", "wctri"="West Coast US", "gmex"="Gulf of Mexico", "sa"="Southeast US", "neus"="Northeast US", "shelf"="Scotian Shelf", "newf"="Newfoundland")
		pretty_col <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','navy','#a65628','salmon','#999999')
		names(pretty_col) <- names(pretty_reg)
	})
}

panLab <- function(){
	pm <- par("mfg")
	nmat <- matrix(1:(prod(pm[3:4])), nr=pm[3], nc=pm[4])
	panelLabel <- LETTERS[nmat[pm[1], pm[2]]]
	return(panelLabel)
}


# =============
# = Main Text =
# =============
# ---- Richness Time Series ----
richness_ts <- function(){
	eval(figure_setup())
	par(mar=c(1.75,1.5,0.25,0.25),mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
	comm_master[,plot(year, reg_rich, col=pretty_col[reg], xlab="Year", ylab="Estimated Richness", pch=20)]
	for(r in 1:length(regs)){
		comm_master[reg==regs[r],j={lines(year, reg_rich, lwd=0.5, col=adjustcolor(pretty_col[reg], 0.5))}]
	}
	comm_master[,lines(year,fitted(lm(reg_rich~year))),by='reg']
	comm_master[,legend("topleft",ncol=1,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.085,-0.01), bty='n')]
	invisible(NULL)
}

# ---- Prevalence Time to Absence ----
rangeSize_absenceTime <- function(pred_var=c("rangeSize","rangeDensity")){
	eval(figure_setup())
	avg_prev_abs <- spp_master[!is.na(ext_dist) & ext_dist!=0,list(rangeSize=mean(propStrata), rangeDensity=mean(propTow_occ)),by=c("reg","ext_dist","stretch_type")]
	# avg_prev_abs <- spp_master[!is.na(ext_dist) & ext_dist!=0,list(rangeSize=mean(propStrata), rangeDensity=mean(propTow_occ)),by=c("reg","ext_dist","stretch_type","spp")]
	
	pred_var <- match.arg(pred_var, several.ok=TRUE) #c("rangeSize","rangeDensity")
	par(mfrow=c(length(pred_var),2), mar=c(1.65,1.65,0.25,0.25), oma=c(0.1,0.1,0.1,0.1), cex=1, ps=8, mgp=c(0.75,0.1,0), tcl=-0.1, ylbias=0.35)
	counter <- 0
	for(v in 1:length(pred_var)){
		pv <- pred_var[v]
		for(st in 1:2){
			counter <- counter + 1L
			t_st <- c("pre_ext","post_col")[st]
			t_xlab <- c("Years before extinction", "Years after colonization")[st]
			t_panel <- LETTERS[counter]
			leg_log <- c(TRUE, FALSE)[st]
	
			xlim <- avg_prev_abs[stretch_type==t_st, range(ext_dist, na.rm=TRUE)]
			if(st==1){
				xlim <- rev(xlim)
			}
			avg_prev_abs[stretch_type==t_st,plot(ext_dist, eval(s2c(pv))[[1]], col=adjustcolor(pretty_col[(reg)], 0.5), xlim=xlim, pch=16, xlab=t_xlab, ylab=c(rangeSize="Range size",rangeDensity="Range density")[pv])]
			avg_prev_abs[stretch_type==t_st,j={
				lines(ext_dist, fitted(lm(eval(s2c(pv))[[1]]~ext_dist)), col=pretty_col[(reg[1])], lwd=1.5)
			},by=c('reg')]

			# avg_prev_abs[stretch_type==t_st,plot(eval(s2c(pv))[[1]], ext_dist, col=adjustcolor(pretty_col[(reg)], 0.5), pch=16, ylab=t_xlab, xlab=c("Range Size","Range Density")[v])]
			# avg_prev_abs[stretch_type==t_st,j={
			# 	lines(eval(s2c(pv))[[1]], fitted(lm(ext_dist ~ eval(s2c(pv))[[1]])), col=pretty_col[(reg[1])], lwd=1.5)
			# },by=c('reg')]
	
			# legend("topleft", legend=t_panel, inset=c(-0.075, -0.03), bty='n', text.font=2)
			mtext(t_panel, side=3, line=-1.5, adj=c(0.05,0.95)[st], font=2, cex=1.25)
			if(counter == 1){
				comm_master[,legend("topright",ncol=1,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.02, -0.03), bty='n', x.intersp=1, y.intersp=0.65)]
			}	
		}
	}
	
	
	invisible(NULL)
}

# ---- Colonization Rate ----
ceRate_map <- function(ce=c("colonization","extinction")){
	ce <- match.arg(ce)
	eval(figure_setup())
	map_layout <- trawl_layout()
	par(mar=c(0.9,0.9,0.25,0.25), mgp=c(0.5,0.075,0), tcl=-0.1, ps=8, cex=1, oma=c(0.85,0.6,1,0.1))
	layout(map_layout)

	map_col <- grDevices::colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
	toRast <- function(p){
		crs_orig <- CRS("+proj=longlat")
		
		x0 = rep(p$xcol, length(p$yrow))
		y0 = rep(p$yrow, each=length(p$xcol))
		# xy <- mapproject(x0, y0, proj="albers", parameters=extent(r)[3:4], orientation=c(mean(y0), mean(x0), 0))
		z = c(t(as.matrix(p)))
		
		r0 <- raster::rasterFromXYZ(cbind(x0, y0, z), crs=crs_orig)
		
		# crs_new <- CRS(paste0("+proj=aea +lat_1=", min(y0), " +lat_2=", max(y0), " +lat_0=", mean(y0), " +lon_0=", mean(x0), " +x_0=0 +y_0=0 +ellps=WGS84 +datum=NAD83 +units=m +no_defs"))
		# proj_to <- projectExtent(r0, crs=crs_new)
		# r <- projectRaster(from=r0, to=proj_to)
		# plot(r, asp=1)
		# map(add=T, proj="albers", parameters=extent(r)[3:4], orientation=c(mean(y0), mean(x0), 0))
		r <- r0
		return(r)
	}

	u_regs <- mapDat[,unique(reg)]
	rs <- mapDat[,una(reg)]
	nr <- length(rs)
	mapPPP_ce <- list()
	for(r in 1:nr){
		td <- mapDat[reg==rs[r]]
		if(ce=="colonization"){
			mapPPP_ce[[r]] <- spatstat::ppp(x=td[,lon], y=td[,lat], marks=td[,n_spp_col_weighted], window=mapOwin[[r]]) # /avgRich
			# mapPPP_ce[[r]] <- spatstat::ppp(x=td[,lon], y=td[,lat], marks=td[,n_spp_col_unique], window=mapOwin[[r]]) # /avgRich
		}else if(ce=="extinction"){
			mapPPP_ce[[r]] <- spatstat::ppp(x=td[,lon], y=td[,lat], marks=td[,n_spp_ext_weighted], window=mapOwin[[r]])
		}
		
		t_idw <- spatstat::Smooth(mapPPP_ce[[r]], hmax=1)
		z <- toRast(t_idw)
		raster::image(z, col=map_col, xlab="", ylab="", asp=1)
		map(add=TRUE, fill=TRUE, col="lightgray")
	
		zl <- range(values(z)*10, na.rm=TRUE)
		switch(rs[r],
			ebs = mapLegend(x=0.05, y=0.25, h=0.375, w=0.025, zlim=zl, cols=map_col, lab.cex=1),
			ai = mapLegend(x=0.985, y=0.3, w=0.02, h=0.75, zlim=zl, cols=map_col, lab.cex=1),
			goa = mapLegend(x=0.985, y=0.15, w=0.02,  zlim=zl, cols=map_col, lab.cex=1),
			wctri = mapLegend(x=0.1, y=0.120, w=0.2, h=0.15, zlim=zl, cols=map_col, lab.cex=1),
			gmex = mapLegend(x=0.95, y=0.2, h=0.375, zlim=zl, cols=map_col, lab.cex=1),
			sa = mapLegend(x=0.95, y=0.15, zlim=zl, cols=map_col, lab.cex=1),
			neus = mapLegend(x=0.95, y=0.15, zlim=zl, cols=map_col, lab.cex=1),
			shelf = mapLegend(x=0.95, y=0.15, zlim=zl, cols=map_col, lab.cex=1),
			newf = mapLegend(x=0.05, y=0.15, h=0.20, zlim=zl, cols=map_col, lab.cex=1)
		)
		switch(rs[r],
			ebs = legend("topright", legend="A", bty='n', text.font=2, inset=c(-0.02,-0.15), cex=1.25, text.col='black'),
			ai = legend("topleft", legend="C", bty='n', text.font=2, inset=c(-0.065,-0.45), cex=1.25, xpd=T),
			goa = legend("topleft", legend="B", bty='n', text.font=2, inset=c(-0.065,-0.06), cex=1.25),
			wctri = legend("top", legend="E", bty='n', text.font=2, inset=c(0,0.05), cex=1.25, text.col='black'),
			gmex = legend("topleft", legend="G", bty='n', text.font=2, inset=c(-0.175,-0.12), cex=1.25, text.col='black'),
			sa = legend("topleft", legend="I", bty='n', text.font=2, inset=c(-0.15,-0.075), cex=1.25, text.col='black'),
			neus = legend("topleft", legend="H", bty='n', text.font=2, inset=c(-0.125,-0.05), cex=1.25, text.col='black'),
			shelf = legend("topleft", legend="D", bty='n', text.font=2, inset=c(-0.1,-0.125), cex=1.25, text.col='black'),
			newf = legend("topright", legend="F", bty='n', text.font=2, inset=c(-0.01,-0.05), cex=1.25)
		)
	}
	if(ce=="colonization"){
		mtext(bquote(Colonization~Rate~(C[w]~~decade^-1)), side=3, outer=TRUE, font=2, line=-0.3)
	}else if(ce=="extinction"){
		mtext(bquote(Extinction~Rate~(E[w]~~decade^-1)), side=3, outer=TRUE, font=2, line=-0.3)
	}
	mtext(bquote(Longitude~(phantom()*degree*E)), side=1, line=0.15, outer=TRUE)
	mtext(bquote(Latitude~(phantom()*degree*N)), side=2, line=-0.4, outer=TRUE)
	invisible(NULL)
}

# ---- Richness vs Range Density ----
# ---- Richness vs Range Size ----
#' Plot richness amd geographic range
#' Plots species richness against two measures of geographic distribution: range size and range density. Size is the proportion of sites a species is present, density is the proportion of tows in occupied sites. 
#' 
#' @param gR0 character vector (can be length 1 or 2) specifying metrics to plot
#' @param leg logical, plot a legend? Default is to plot the legend associating each region with a color.
#' @param legPan integer, which panel should contain the legend? Default is 2, and if \code{gR0} is default, the second panel which is range density (because it fits better). 
#' @param panLab logical, add letter for panel label? Default is TRUE, to add A in bottom-left of first panel, B in second
#' 
#' @details
#' All values taken from \code{\link{comm_master}}.
#' 
#' The measures are species- and year-specific; for this figure they are first averaged across all years for each species-region combination. This first average yields a value that is intended to be "characteristic" of that species in that region. In each year, the characteristic values of each species are then averaged to form a community descriptor. That community descriptor can change among years because the identity of the species present can change, thus the species characteristics being averaged can change.
#' 
#' @note The above description of these metrics is subject to change. Check manuscript_data.R to be sure.
#' 
#' @return NULL returned invisibly
#' 
#' @seelaso \code{\link{comm_master}} \code{\link{spp_master}}
#' 
#' @examples
#' rich_geoRange()
#' 
#' @export
rich_geoRange <- function(gR0=c("size", "density"), leg=TRUE, legPan=2, panLab=TRUE){
	eval(figure_setup())
	
	gR0 <- match.arg(gR0, several.ok=TRUE)
	if(length(gR0)>1){
		mfr <- rbLib::auto.mfrow(length(gR0), tall=TRUE)
		par(mfrow=mfr, mar=c(1.75,1.5,0.25,0.25),mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
	}else{
		par(mar=c(1.75,1.5,0.25,0.25),mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
	}
	for(g in 1:length(gR0)){
		gR <- switch(gR0[g],
			density = "propTow_occ_avg",
			size = "propStrata_avg_ltAvg"
		)	
		xlab <- switch(gR0[g],
			density = "Community-level range density",
			size = "Community-level range size"
		)
		comm_master[,plot(eval(s2c(gR))[[1]], reg_rich, col=adjustcolor(pretty_col[reg],0.5), xlab=xlab, ylab="Species richness", pch=16)]
		ur <- comm_master[,unique(reg)]
		mod_expr <- bquote(fitted(lm(reg_rich~eval(s2c(gR))[[1]])))
		for(r in 1:length(ur)){
			comm_master[reg==ur[r],lines(eval(s2c(gR))[[1]],eval(mod_expr))]
		}
		if(leg & g==legPan){
			if(gR0[g]=="density"){
				comm_master[,legend("topright",ncol=2,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.02, -0.02), bty='n')]
			}else if(gR0[g]=="size"){
				comm_master[,legend("topleft",ncol=2,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(0.03, -0.02), bty='n', x.intersp=0.15, y.intersp=0.65)]
			}
			
		}
		if(panLab){
			# xy <- par()$usr[c(1,3)]
			# xy <- xy + c(sign(xy[1]), sign(xy[2]))*c(0.1, 0.15)*xy
			# text(xy[1], xy[2], labels=c("A","B")[g], font=2)
			mtext(c("A","B")[g], side=1, line=-1.5, adj=0.05, font=2, cex=1.25)
		}	
	}
	invisible(NULL)
}





# ==============
# = Supplement =
# ==============
# ---- MSOM Richness vs Naive Richness Scatterplot ----
naive_msom_scatter <- function(){
	eval(figure_setup())
	par(mfrow=c(3,3), mar=c(1.25,1.0,0.5,0.1), oma=c(0.35,0.5,0.1,0.1), mgp=c(0.25,0.1,0), tcl=-0.1, ps=8, cex=1)
	for(r in 1:length(regs)){
		comm_master[reg==regs[r],j={
			plot(naive_rich, reg_rich, main=pretty_reg[una(reg)], type='p', xlab="", ylab="", cex.main=1)
			abline(a=0, b=1)
		}]
	}

	mtext("MSOM Richness", side=2, line=-0.2, outer=TRUE)
	mtext("Naive Richness", side=1, line=-0.5, outer=TRUE)
	invisible(NULL)
}

# ---- Category Barplot ----
categ_barplot <- function(){
	eval(figure_setup())
	categ_table <- t(spp_master[!duplicated(paste(reg,ce_categ,spp)), table(reg, ce_categ)])[c(4,1,2,3),]
	colnames(categ_table) <- pretty_reg[colnames(categ_table)]
	colnames(categ_table) <- gsub("^(.*) (.*)$", "\\1\n\\2", colnames(categ_table))
	par(cex=1, mar=c(3,2,1,0.1), ps=8)
	bp <- barplot(categ_table, beside=T, legend=T, names.arg=rep("",ncol(categ_table)), args.legend=list(bty='n'))
	text(colMeans(bp)-1, -11, labels=colnames(categ_table), srt=45, xpd=TRUE)
	invisible(NULL)
}

# ---- Colonization / Extinction Time Series ----
col_ext_ts <- function(){
	eval(figure_setup())
	par(mfrow=c(3,3), mar=c(1.25,1.0,0.5,0.1), oma=c(0.35,0.5,0.1,0.1), mgp=c(0.25,0.1,0), tcl=-0.1, ps=8, cex=1)
	for(r in 1:length(regs)){
		comm_master[reg==regs[r],j={
			ylim=range(c(n_col,n_ext));
			plot(year, n_col, main=pretty_reg[una(reg)], type='l', col='blue', ylim=ylim, xlab="", ylab="", cex.main=1)
			lines(year, n_ext, col='red')
		}]
	}
	mtext("Colonizations or Extinctions", side=2, line=-0.2, outer=TRUE)
	mtext("Year", side=1, line=-0.5, outer=TRUE)
	invisible(NULL)
}


# ---- Neighborhood used in Moran's I ----
nb_moranI <- function(ce=c("colonization", "extinction")){
	eval(figure_setup())
	map_layout <- trawl_layout()
	par(mar=c(0.25,0.25,0.25,0.25), mgp=c(0.25,0.075,0), tcl=-0.1, ps=8, cex=1, oma=c(0.1,0.1,0.1,0.1))
	layout(map_layout)
	pretty_reg <- c("ebs"="E. Bering Sea", "ai"="Aleutian Islands", "goa"="Gulf of Alaska", "wctri"="West\nCoast\nUS", "gmex"="Gulf of Mexico", "sa"="Southeast US", "neus"="Northeast US", "shelf"="Scotian Shelf", "newf"="Newfoundland")

	rs <- mapDat[,una(reg)]
	nr <- length(rs)
	for(r in 1:nr){
		t_lac <- localAC[[ce]][[rs[r]]]
		plot(mapOwin[[rs[r]]], coords=t_lac$I[,list(lon,lat)], add=FALSE, main="")
		box()
		if(rs[r]=='wctri'){
			mtext(pretty_reg[rs[r]], side=3, line=-3)
		}else{
			mtext(pretty_reg[rs[r]], side=3, line=-1)
		}
		plot(t_lac$nb, t_lac$I[,list(lon,lat)], add=TRUE, col=adjustcolor('black', 0.5), cex=0.8, lwd=0.5)
		sig_lac <- t_lac$I[,lI_pvalue]<0.05
		if(any(sig_lac)){
			zl <- range(t_lac$I[sig_lac,Ii], na.rm=TRUE)
			t_col <- rbLib::zCol(256,t_lac$I[sig_lac,Ii])
			map_col <- rbLib::zCol(6, 1:6)
			points(x=t_lac$I[sig_lac,lon], y=t_lac$I[sig_lac,lat], bg=t_col, pch=21, cex=1.1)
			switch(rs[r],
				ebs = mapLegend(x=0.05, y=0.25, h=0.375, w=0.025, zlim=zl, cols=map_col, lab.cex=1),
				ai = mapLegend(x=0.985, y=0.3, w=0.02, h=0.75, zlim=zl, cols=map_col, lab.cex=1),
				goa = mapLegend(x=0.985, y=0.15, w=0.02,  zlim=zl, cols=map_col, lab.cex=1),
				wctri = mapLegend(x=0.1, y=0.120, w=0.2, h=0.15, zlim=zl, cols=map_col, lab.cex=1),
				gmex = mapLegend(x=0.95, y=0.2, h=0.375, zlim=zl, cols=map_col, lab.cex=1),
				sa = mapLegend(x=0.95, y=0.15, zlim=zl, cols=map_col, lab.cex=1),
				neus = mapLegend(x=0.95, y=0.15, zlim=zl, cols=map_col, lab.cex=1),
				shelf = mapLegend(x=0.95, y=0.15, zlim=zl, cols=map_col, lab.cex=1),
				newf = mapLegend(x=0.05, y=0.15, h=0.20, zlim=zl, cols=map_col, lab.cex=1)
			)
		}else{
			# plot(x=locs[,1], y=locs[,2], col='blue')
		}
	}
	invisible(NULL)
}


# ---- community range density vs range size ----
rangeSizeDens <- function(){
	eval(figure_setup())
	
	par(mfrow=c(2,1), mar=c(1.75,1.5,0.25,0.25),mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
	range_reg_spp <- spp_master[,list(reg, year, density=propTow_occ, size=propStrata)]
	range_reg_spp[,plot(size, density, col=adjustcolor(pretty_col[reg],0.1), xlab="Range size", ylab="Range density", pch=16)]
	ur <- range_reg_spp[,unique(reg)]
	for(r in 1:length(ur)){
		td <- range_reg_spp[reg==ur[r] & !is.na(density) & !is.na(size)]
		setorder(td, size)
		lf <- td[,fitted(loess(density~size))]
		td[, lines(size, lf, col="white", lwd=4)]
		td[, lines(size, lf, col=pretty_col[reg], lwd=2)]
	}
	mtext(c("A"), side=1, line=-1.5, adj=0.95, font=2, cex=1.25)
	
	range_reg <- comm_master[,list(reg, year, rich=reg_rich, density=propTow_occ_avg, size=propStrata_avg_ltAvg)]
	range_reg[,plot(size, density, col=adjustcolor(pretty_col[reg],0.5), xlab="Community range size", ylab="Community range density", pch=16)]
	comm_master[,legend("topright",ncol=2,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.02, -0.02), bty='n')]
	mtext(c("B"), side=1, line=-1.5, adj=0.95, font=2, cex=1.25)

	invisible(NULL)
	
}

ceEventRange <- function(pred_vars = c("mean_size", "mean_density")){
	eval(figure_setup())
	pred_vars <- match.arg(pred_vars, several.ok=TRUE)
	par(mfrow=c(length(pred_vars),1), mar=c(1.75,1.5,0.25,0.25),mgp=c(0.85,0.1,0), tcl=-0.1, cex=1, ps=8)
	
	range_ceEvents <- spp_master[present==1,.SD[,list(mean_size=mean(propStrata), mean_density=mean(propTow_occ), total_colExt=sum(col+ext)),by='spp'],by='reg']
	ur <- range_ceEvents[,unique(reg)]
	
	for(v in 1:length(pred_vars)){
		pv <- pred_vars[v]
		
		range_ceEvents[,plot(eval(s2c(pv))[[1]], total_colExt, col=adjustcolor(pretty_col[reg],0.5), xlab=c("mean_size"="Range Size", "mean_density"="Range Density")[pv], ylab="Total Colonizations and Extinctions", pch=16)]
		for(r in 1:length(ur)){
			td <- range_ceEvents[reg==ur[r] & !is.na(mean_density) & !is.na(mean_size)]
			setorderv(td, pv)
			lf <- td[,fitted(loess(total_colExt~eval(s2c(pv))[[1]]))]
			td[, lines(eval(s2c(pv))[[1]], lf, col="white", lwd=4)]
			td[, lines(eval(s2c(pv))[[1]], lf, col=pretty_col[reg], lwd=2)]
		}
		if(v==1 | length(pred_vars)==1){
			range_ceEvents[,legend("topright",ncol=2,legend=pretty_reg[una(reg)],text.col=pretty_col[una(reg)], inset=c(-0.02, -0.02), bty='n')]
		}	
		if(length(pred_vars)>1){
			mtext(c("A","B")[v], side=1, line=-2, adj=0.95, font=2, cex=1.25)
		}
	}
	invisible(NULL)
}


