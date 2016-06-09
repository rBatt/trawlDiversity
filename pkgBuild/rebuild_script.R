
library(devtools)
setwd("~/Documents/School&Work/pinskyPost/trawl")

# =========================
# = Uninstall & Reinstall =
# =========================
# shoudl be repeated after script is run, or just done at end
# putting here so easy to reinstall w/o rebuilding
remove.packages('trawlDiversity')
install('trawlDiversity', upgrade_dependencies=FALSE)
library(trawlDiversity)


# ==================
# = Options/ Setup =
# ==================
show_up_choice <- FALSE

reg_depthStratum <- c(
	"ebs" = 500,
	"ai" = 100,
	"goa" = 500,
	"wctri" = 100, 
	"wcann" = 500, 
	"gmex" = 500, 
	"sa" = 500, 
	"neus" = 500, 
	"shelf" = 500, 
	"newf" = 500
)

reg_tolFraction <- c(
	"ebs" = 0,
	"ai" = 0.15,
	"goa" = 0,
	"wctri" = 0.15, 
	"wcann" = 0.15, 
	"gmex" = 0.15, 
	"sa" = 0.15, 
	"neus" = 0.15, 
	"shelf" = 0.15, 
	"newf" = 0.15
)

sus_reg <- list()
data_all_list <- list()

# ===================
# = Show Up Species =
# ===================
poss_regs <- c("ebs", "ai", "goa", "wctri", "wcann", "gmex", "sa", "neus", "shelf", "newf")
for(r in 1:length(poss_regs)){
	t_reg <- poss_regs[r]
	t_dat <- trim_msom(t_reg, gridSize=0.5, depthStratum=reg_depthStratum[t_reg], tolFraction=reg_tolFraction[t_reg], grid_stratum=TRUE, plot=FALSE, cull_show_up=FALSE)
	sus_reg[[r]] <- get_show_up_spp(t_dat, n_yrs_start=2, n_yrs_remain=5, n_yrs_miss_after=0)
	
	if(!show_up_choice){
		# If not culling show up spp, can re-use t_dat for data_all creation
		data_all_list[[r]] <- t_dat
	}
}

show_up_spp <- rbindlist(sus_reg)
save(show_up_spp, file="trawlDiversity/data/show_up_spp.RData")


# ============
# = Data All =
# ============
# only need to loop through if choosing to cull show up spp
# otherwise can just use the same result of trim_msom used to create show up spp
if(show_up_choice){
	for(r in 1:length(poss_regs)){
		t_reg <- poss_regs[r]
		data_all_list[[r]] <- trim_msom(t_reg, gridSize=0.5, depthStratum=reg_depthStratum[t_reg], tolFraction=reg_tolFraction[t_reg], grid_stratum=TRUE, plot=FALSE, cull_show_up=show_up_choice)
	}
}
data_all <- rbindlist(data_all_list, fill=TRUE)
setkey(data_all, reg, year, stratum, haulid, spp)
save(data_all, file="trawlDiversity/data/data_all.RData")
# dev.new(); par(mfrow=c(5,2), mar=c(1.85, 1.85, 1, 0.1),mgp=c(0.5,0.1,0), tcl=-0.1, cex=1, ps=8); data_all[,plot(.SD[,list(x=una(year), y=lu(spp)),by='year'][,list(x,y)], type="o", main=una(reg)),by='reg']




