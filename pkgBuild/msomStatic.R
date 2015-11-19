load_all("trawlData")
load_all("trawl/trawlDiversity")

ebs.a2 <- ebs.a[,list(year=year, spp=spp, stratum=stratum, K=K, abund=abund, btemp=btemp, doy=yday(datetime))]

Yc <- trawlCast(x=ebs.a2, formula="spp~K~stratum~year", valueName="abund", grandNamesOut=c("s","k","j","t"))
dimnames(Yc) <- lapply(dimnames(Yc), function(x)1:length(x))
Ym0 <- reshape2::melt(Yc)
Ym <- Ym0[complete.cases(Ym0),]



