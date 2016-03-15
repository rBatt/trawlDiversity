library("rstan")
library("trawlDiversity")
library("rbLib")
library("R2jags")
library("maps")

reg_file <- c(
	"msomStatic_norv_1yr_ebs_jags_start2016-03-10_16-44-44_r1.RData", # ebs 6k iter, converged
	# "msomStatic_norv_1yr_ai_jags_start2016-03-11_16-01-45_r2.RData", #ai 6k iter, didn't converge
	"msomStatic_norv_1yr_ai_jags_start2016-03-06_19-10-49_r2.RData", # ai 30k iter, didn't converge
	# "msomStatic_norv_1yr_goa_jags_start2016-03-11_17-55-00_r3.RData", # goa 6k iter, didn't converge
	"msomStatic_norv_1yr_goa_jags_start2016-03-07_06-21-58_r3.RData", # goa 30k iter, didn't converge
	# "msomStatic_norv_1yr_wctri_jags_start2016-03-11_21-14-14_r4.RData", # wctri 6k iter, didn't converge (also too few 0 spp)
	"msomStatic_norv_1yr_wctri_jags_start2016-03-08_10-53-18_r4.RData", # wctri 30k iter, didn't converge (better than most tho); also, needs more than 100 0 spp
	# "msomStatic_norv_1yr_wcann_jags_start2016-03-11_23-18-31_r5.RData", # wcann 6k iter, didn't converge, too few 0 sp (50)
	"msomStatic_norv_1yr_wcann_jags_start2016-03-07_13-25-01_r5.RData", # wcann 30k iter, didn't converge, too few 0 spp
	"msomStatic_norv_1yr_gmex_jags_start2016-03-12_02-53-26_r6.RData", # gmex, 30k iter, didn't converge, too few 0 spp
	# "msomStatic_norv_1yr_sa_jags_start2016-03-12_04-26-53_r7.RData", # sa, 6k iter, almost converged, but too few n0 spp
	"msomStatic_norv_1yr_sa_jags_start2016-03-06_16-15-08_r7.RData", # sa, 30k iter, very nearly converged, but too few n0 spp
	"msomStatic_norv_1yr_neus_jags_start2016-03-12_06-01-10_r8.RData", # neus, 6k iter, mostly converged, but too few n0 spp
	# "msomStatic_norv_1yr_shelf_jags_start2016-03-12_12-00-07_r9.RData", # shelf, 6k iter, very nearly converged
	"msomStatic_norv_1yr_shelf_jags_start2016-03-06_14-03-15_r9.RData", # shelf, 30k iter, better than 6k, but not perfect
	"msomStatic_norv_1yr_newf_jags_start2016-03-12_13-30-24_r10.RData" # newf, 6k iter, really good mixing, mild lack of stationary
	
)

p <- list()
# for(i in 1:length(reg_file)){
for(i in 10:10){	
	load(paste0("trawlDiversity/pkgBuild/results/", reg_file[i]))
	
	reg_results_ind <- which(sapply(rm_out, function(x)!is.null(x)))
	stopifnot(length(reg_results_ind) == 1)
		
	rm_out <- rm_out[[reg_results_ind]]
	
	
	p[[i]] <- process_msomStatic(rm_out)
	
	
}



# ===================
# = Unpack a region =
# ===================
rd <- p[[10]]$rd
processed <- p[[10]]$processed
bt <- p[[10]]$bt
colonization <- p[[10]]$colonization
param_iters <- p[[10]]$param_iters
ab <- p[[10]]$ab

reg <- processed[,una(reg)]

naive_rich <- processed[,naive_rich, by='year']
reg_rich <- processed[,reg_rich, by='year']
bt_ann <- bt[,list(bt_ann=mean(bt)), by='yr']



# ===========
# = Figures =
# ===========
fig1_name <- paste0("richness_bt_timeSeries_", reg, ".png")
# png(file.path("trawlDiversity/pkgBuild/figures",fig1_name), width=3.5, height=6, units="in", res=200)
dev.new(width=3.5, height=6)
par(mfrow=c(3,1), mar=c(1.75,1.5,0.25,0.25), oma=c(0.1,0.1,0.75,0.1), mgp=c(0.75,0.1,0), tcl=-0.1, ps=8, cex=1)
plot(naive_rich, type="o", ylab="Naive Region Richness", xlab="Year")
plot(reg_rich, type="o", xlab="Year", ylab="MSOM Region Richness")
plot(bt_ann, type="o", xlab="Year", ylab="Annual Mean Bottom Temperature")
mtext(reg, side=3, line=0, outer=TRUE, font=2)
# dev.off()

fig2_name <- paste0("richness_bt_scatter_", reg, ".png")
# png(file.path("trawlDiversity/pkgBuild/figures",fig2_name), width=3.5, height=5, units="in", res=200)
dev.new(width=3.5, height=5)
par(mfrow=c(2,1), mar=c(1.75,1.5,0.25,0.25), oma=c(0.1,0.1,0.75,0.1), mgp=c(0.75,0.1,0), tcl=-0.1, ps=8, cex=1)
plot(processed[,list(bt_ann,naive_rich)], type="p", ylab="Naive Region Richness", xlab="Annual Mean Bottom Temperature")
plot(processed[,list(bt_ann,reg_rich)], type="p", ylab="MSOM Region Richness", xlab="Annual Mean Bottom Temperature")
mtext(reg, side=3, line=0, outer=TRUE, font=2)
# dev.off()

fig3_name <- paste0("btempMap_", reg, ".png")
# png(file.path("trawlDiversity/pkgBuild/figures",fig3_name), width=8, height=3, units="in", res=200)
dev.new(width=8, height=3)
par(mfrow=auto.mfrow(bt[,lu(year)]), oma=c(0.1,0.1, 1,0.1), mar=c(1,1,0.1,0.1), mgp=c(0.75,0.1,0), tcl=-0.15, cex=1, ps=8)
bt[,j={
	plot(lon, lat, type="n")
	map(add=TRUE)
	points(lon, lat, col=bt_col, pch=20)
	mtext(unique(year), side=3, adj=0.1, line=-0.75, font=2)
}, by="yr"]
mtext(paste(reg, "Bottom Temperature"), outer=TRUE, side=3, line=0.25, font=2)
# dev.off()

fig4_name <- paste0("traceplot_", reg, ".png")
# png(file.path("trawlDiversity/pkgBuild/figures",fig4_name), width=12, height=6, units="in", res=200)
dev.new(width=12, height=6)
n_yrs <- param_iters[,lu(year)]
par(mfrow=c(length(pars_trace), n_yrs), oma=c(1,1, 1,0.1), mar=c(0.5,0.5,0.1,0.1), mgp=c(0.25,0.1,0), tcl=-0.1, cex=1, ps=6)
for(h in 1:length(pars_trace)){
	for(i in 1:n_yrs){
		t_yr <- param_iters[,una(year)][i]
		t_iters <- param_iters[year==t_yr]
		mytrace(t_iters, pars=pars_trace[h], lang=lang, xaxt='n')
		if(i == 1){
			mtext(pars_trace[h], side=2, line=0.75)
		}
	}
}
# dev.off()


# ---- look for correlation in posterior distribution of parameteres ----
# i.e., do parameter values covary
dev.new()
pairs(param_iters[year==param_iters[,una(year)][1], eval(s2c(pars_trace))])


# ---- Examine Colonizations[t] vs unobserved_spp[t-1] ----
dev.new()
processed[,plot(unobs_rich[-length(unobs_rich)], n_col[-1], xlab="Unobserved species present last year", ylab="Species colonizing this year")]
abline(a=0, b=1)


# ---- Number of Colonizations per Stratum ----
dev.new(width=6, height=3)
par(mfrow=c(1,2), mar=c(1.5,1.5,0.1,0.1), mgp=c(1,0.1,0), tcl=-0.1, ps=8, cex=1)

# site-specific colonizations from data
colonization$n_spp_col_weighted_tot[,plot_space(lon, lat, n_spp_col_weighted, TRUE, pch=19)]
map(add=TRUE, fill=TRUE, col="white")

# smoothed map for convex hull of site observations
colonization$n_spp_col_weighted_tot[,plot_space(lon, lat, n_spp_col_weighted, pch=19)]
map(add=TRUE, fill=TRUE, col="white")



# ---- Function to plot time series of species parameter ----
plot_ab <- function(X, t_spp){
	
	si <- sppImg(t_spp)
	if(!is.null(si)){
		par(new=T)
	}
	
	plot(X[,year], X[,value], col=adjustcolor('gray', 0.01), cex=0.5, pch=21, bg=adjustcolor('white',0.01))
	mu <- X[,list(mu=mean(value)),by="year"]
	mu[,lines(year, mu, lwd=2, col='gray')]
	mu[,lines(year, mu, lwd=1, col='white')]
	if(is.null(si)){
		common_name <- spp.key[spp==t_spp, una(common)]
		mtext(paste(t_spp, common_name, sep="\n"), side=3)
	}
	
}

# ---- Time Series of Species Detection Parameter (Beta) ----
dev.new()
par(mfrow=auto.mfrow(ab[,lu(spp)]), mar=c(1,1,1,0.1), ps=6, mgp=c(0.6,0.1,0), tcl=-0.1, cex=1)
yc <- ab[,zCol(lu(year), una(year))]
names(yc) <- ab[,una(year)]

u_spp <- ab[,una(spp)]
n_spp <- ab[,lu(spp)]
for(s in 1:n_spp){
	t_s <- u_spp[s]
	t_dat <- ab[par=="beta" & (spp)==t_s]
	plot_ab(X=t_dat, t_spp=t_s)
}


# ---- Time Series of Presence Parameters (alpha) ----
dev.new()
par(mfrow=auto.mfrow(ab[,lu(spp)]), mar=c(1,1,1,0.1), oma=c(0.1,0.1,1.5,0.1), ps=6, mgp=c(0.6,0.1,0), tcl=-0.1, cex=1)
u_spp <- ab[,una(spp)]
n_spp <- ab[,lu(spp)]
for(s in 1:n_spp){
	t_s <- u_spp[s]
	t_dat <- ab[par=="alpha" & ab_ind==1 & (spp)==t_s]
	plot_ab(X=t_dat, t_spp=t_s)
}
mtext(paste0(reg, "Alpha[1] (presence intercept)"))



	