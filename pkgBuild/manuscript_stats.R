
library('trawlDiversity')
library('rbLib')
library("EnvStats")

# ========================
# = Richness Time Series =
# ========================
# ---- mean richness ----
# smallest and largest long-term averages
comm_master[,mean(reg_rich), by='reg'][reg%in%c("gmex","shelf")]

# long-term variability
comm_master[,sd(reg_rich),by='reg']#[,mean(V1)]
comm_master[,sd(reg_rich),by='reg'][,mean(V1)]





