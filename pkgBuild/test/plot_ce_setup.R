plot_ce_setup <- function(spp2use, width.max=12, height.max=7, max_spp_columns=15, nPlots=3){
	ncolspp <- max_spp_columns
	
	
	if(length(spp2use) > ncolspp){
		nrowspp <- ceiling(length(spp2use)/ncolspp)
		lmat0 <- matrix(rep(1:(ncolspp*nrowspp),each=nPlots), nrow=nrowspp*nPlots, ncol=ncolspp)
		lmat <- lmat0
		lmat[] <- 1:prod(dim(lmat))
		lmat <- lmat[,apply(lmat < (length(spp2use)*nPlots), 2, function(x)any(x))]
	
		fwc_ldim <- dim(lmat)
		fwc_w <- min(max(1*fwc_ldim[2], 1.5), width.max) # width
		fwc_h <- min(max(fwc_w * fwc_ldim[1] / fwc_ldim[2] * 1.1, 3.0), height.max)
		fwc_dim <- c("width"=fwc_w, "height"=fwc_h)
		dev.new(width=fwc_dim[1], height=fwc_dim[2])
		layout(mat=lmat, heights=rep(c(1.5,1,1), nrowspp))
		par(mar=c(1.5,1,0.25,0.1) ,oma=c(0.1,0.1,0.1,0.1), ps=6, mgp=c(0.5, 0.01, 0), tcl=-0.01, cex=1)
			
	}else{
		fwc_mfr <- c(nPlots, length(spp2use))
		fwc_w <- min(max(1*fwc_mfr[2], 1.5), width.max) # width
		fwc_h <- min(max(fwc_w * fwc_mfr[1] / fwc_mfr[2] * 1.1, 3.0), height.max)
		fwc_dim <- c(width=fwc_w, height=fwc_h)
		dev.new(width=fwc_dim[1], height=fwc_dim[2])
		par(mfcol=fwc_mfr, mar=c(1,1,1,0.1) ,oma=c(0.1,0.1,1.5,0.1), ps=6, mgp=c(0.6, 0.1, 0), tcl=-0.1, cex=1)
	}
	
	return(c(fig_dim=fwc_dim))

}