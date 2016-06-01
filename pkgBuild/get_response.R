
get_response <- function(p, n_grid=20){
	alpha_unscale<-p$alpha_unscale
	alpha_unscale <- alpha_unscale[,.SD[sample(1:nrow(.SD),20)],by=c("spp","year")]
	nIter <- alpha_unscale[,.N,by=c("spp","year")][,una(N)]

	bt <- p$bt

	# n_grid <- 20

	d_vals <- bt[,seq(min(depth), max(depth), length.out=n_grid)]
	b_vals <- bt[,seq(min(bt, na.rm=TRUE), max(bt, na.rm=TRUE), length.out=n_grid)]
	vals <- expand.grid(b=b_vals, d=d_vals)

	pred_resp <- function(b, d, alpha){
		a <- matrix(alpha, nrow=5)
		x <- as.matrix(cbind(1, b, b^2, d, d^2))
		x%*%a
	}

	pr <- pred_resp(vals[,"b"], vals[,"d"], t(alpha_unscale[,list(a1,a2,a3,a4,a5)]))

	ac <- as.character
	spp_year0 <- alpha_unscale[,paste(spp,year),by=c("spp","year")]
	spp_year <- spp_year0[,V1]
	dn <- list(btemp=ac(b_vals),depth=ac(d_vals),iter=NULL, sppYear=spp_year)

	pr2 <- array(pr, dim=c(n_grid, n_grid, nIter, ncol(pr)/nIter), dimnames=dn)
	pr3 <- apply(pr2, c(1,2,4), mean)

	spp_3 <- spp_year0[,spp]
	year_3 <- spp_year0[,year]
	
	return(list(responses=pr3, spp=spp_3, year=year_3))
}

# alpha_unscale<-p[[1]]$alpha_unscale
# alpha_unscale <- alpha_unscale[,.SD[sample(1:nrow(.SD),20)],by=c("spp","year")]
# nIter <- alpha_unscale[,.N,by=c("spp","year")][,una(N)]
#
bt <- p[[9]]$bt
#
# n_grid <- 20
#
# d_vals <- bt[,seq(min(depth), max(depth), length.out=n_grid)]
# b_vals <- bt[,seq(min(bt, na.rm=TRUE), max(bt, na.rm=TRUE), length.out=n_grid)]
# vals <- expand.grid(b=b_vals, d=d_vals)
#
# pred_resp <- function(b, d, alpha){
# 	a <- matrix(alpha, nrow=5)
# 	x <- as.matrix(cbind(1, b, b^2, d, d^2))
# 	x%*%a
# }
#
# pr <- pred_resp(vals[,"b"], vals[,"d"], t(alpha_unscale[,list(a1,a2,a3,a4,a5)]))
#
# ac <- as.character
# spp_year0 <- alpha_unscale[,paste(spp,year),by=c("spp","year")]
# spp_year <- spp_year0[,V1]
# dn <- list(btemp=ac(b_vals),depth=ac(d_vals),iter=NULL, sppYear=spp_year)
#
# pr2 <- array(pr, dim=c(n_grid, n_grid, nIter, ncol(pr)/nIter), dimnames=dn)
# pr3 <- apply(pr2, c(1,2,4), mean)
#
#
# spp_3 <- spp_year0[,spp]
# year_3 <- spp_year0[,year]

resp_out <- get_response(p[[9]], 20)
pr3 <- resp_out$responses
spp_3 <- resp_out$spp
year_3 <- resp_out$year

u_spp <- unique(spp_3)

for(i in 1:length(u_spp)){
	t_spp <- u_spp[i]

	sy_which <- which(spp_3==t_spp)

	dev.new()
	par(mfrow=auto.mfrow(length(sy_which)+1), oma=c(1,1,0.1,0.1))
	
	par(mar=c(0.1,0.1,1,0.1), ps=6, cex=1)
	tsi <- sppImg(t_spp)
	if(is.null(tsi)){plot(1, type="n")}

	for(j in 1:length(sy_which)){
		t_sy <- sy_which[j]
		t_yr <- year_3[t_sy]
		z <- plogis(pr3[,,t_sy])
		x <- as.numeric(rownames(pr3))
		y <- as.numeric(colnames(pr3))
	
		par(mar=c(2,2,1.75,0.1), ps=8, cex=1, mgp=c(0.5,0.1,0), tcl=-0.1)
		fields::image.plot(z=z, x=x, y=y, main=t_yr, xlab="", ylab="")
		par(mar=c(2,2,1.75,0.1), ps=8, cex=1, mgp=c(0.5,0.1,0), tcl=-0.1)
		bt[(year)==t_yr,points(bt, depth)]

	}
	mtext("Depth (m)", side=2, outer=TRUE, line=0)
	mtext("Temperature (C)", side=1, outer=TRUE, line=0)
}

