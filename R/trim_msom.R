#' Trim trawl data for msom
#' 
#' Performs the \code{trawlTrim} function the clean region data, then trims down to species (not genera) that have been observed at least 10 times. Does aggregating within a haul/ spp using \code{trawlAgg}. Adds an abundance column, which is actually just 0 if wtcpue == 0, and 1 if wtcpue > 0.
#' 
#' @param reg Character name of region or data.table to be passed to \code{\link{trawlTrim}}.
#' @param gridSize numeric grid size to be passed to \code{\link{ll2strat}}
#' @param grid_stratum logical, default TRUE; whether or not the strata should be defined on a grid, or used original stratum definition
#' @param depthStratum The depth range to use to define stratum; the default does not use depth at all. If using, a reasonable value might be 100.
#' @param tolFraction The fraction of years that a stratum can not be sampled yet still be included in output. Default is 1/3, indicating that a stratum will be included in output so long as it is sampled for at least 2/3 of years. Actual number of years tolerated is rounded up. This value determines the strat_tol argument in \code{\link{check_strat}}
#' @param plot Logical, default FALSE; for \code{check_strat}, should the stratum tolerance be plotted?
#' @param cull_show_up Logical, default \code{FALSE}; should spp in \code{\link{show_up_spp}} be removed when trimming?
#' 
#' @return
#' A data.table 
#' 
#' @export
trim_msom <- function(reg, gridSize=1, grid_stratum=TRUE, depthStratum=NULL, tolFraction=1/3, plot=FALSE, cull_show_up=FALSE){
		#
	# reg = 'neus'
	# gridSize = 0.5
	# grid_stratum = TRUE
	# depthStratum = 100
	# tolFraction = 0.15
	# plot=TRUE
	
	# notes on depthStratum
	# neus doesn't need it, or set to 500
	# shelf doesn't need it, or set to 500
	# goa doesn't need it, or set to 500 (200 isn't bad either)
	# gmex doesn't change much between 1000 or 100 ...
	# ai loses about 20 strata going from 100 to 500 (keep at 100)
	# ebs is basically same; loses 3 strata going from 500 to 100
	# wctri loses almost half strata going from 100 to 500 (keep at 100)
	# wcann loses 12 strata and gets gaps in coverage going from 500 to 100 (recommend 500)
	# sa doesn't change at all between 100 and 500
	# newf loses about 80 strata going from 500 to 100 and gets gaps in coverage (recommend 500)
	
	if(reg == "shelf"){ # have to do this here b/c trimming drops
		
	}
	
	# use default trimming
	if(reg == "shelf"){
		X.t <- trawlTrim(reg, c.add=c("CODE"))
	}else{
		X.t <- trawlTrim(reg)
	}
	
	
	
	if(grid_stratum){
		# define stratum on gridSizeº grid
		X.t[,stratum:=ll2strat(lon, lat, gridSize=gridSize)]
		if(!is.null(depthStratum)){
			X.t[,stratum:=paste(stratum, roundGrid(depth, depthStratum))]
		}
	}
	
	
	# drop strata that weren't sampled every year
	if(reg == "ai"){
		X.t <- X.t[(year)>1900,]
	}
	if(reg == "ebs"){
		X.t <- X.t[(year)>1983,]
	}
	if(reg == "gmex" | reg == "neus"){
		X.t <- X.t[(year)!=2015,]
	}
	if(reg == "gmex"){
		X.t <- X.t[(year)>1983 & (year)<2001]
	}
	if(reg == "neus"){
		X.t <- X.t[(year)>1981 & (year)<2014]
	}
	if(reg == "newf"){
		X.t <- X.t[(year)>1995,]
	}
	if(reg == "sa"){
		X.t <- X.t[(year)>=1990]
	}
	if(reg == "shelf"){
		X.t <- X.t[(year)!=2011 & (year)>1950,]
	}
	if(reg == "wcann"){
		X.t <- X.t[(year)>2003,]
	}
	
	# ---- constratin/ standardize time of year (day of year) sampling occurred ----
	yd <- X.t[,yday(datetime)]
	X.t <- switch(reg,
		ebs = X.t[yd <= 220],
		# ai = X.t[yd >= 170 & yd <= 225], # sampling periods really overlap between 175 and 225
		# goa = X.t[yd <= 210],
		# wc = X.t[yd >= 170 & yd <= 240],
		# wctri = X.t[yd >= 170 & yd <= 240],
		gmex = X.t,
		sa = X.t[yd >= 150 & yd <= 250],
		neus = X.t,
		shelf = X.t[yd >=  180 & yd <= 215],
		newf = X.t[ yd >= 250 | yd <= 28],

		X.t # just returns X.t if reg isn't matched (e.g., if reg was wcann or sgulf)

	)

	# ---- Cut out methodologically suspicious species ----
	if(reg == "ebs"){
		bad_spp_ebs <- c(
			"Lepidopsetta polyxystra", # showed up in 1996 and was prolific ever since; also line 182 in Kotwicki & Lauth
			"Bathyraja parmifera", # line 182 Kotwicki & Lauth draft MS
			# "Chrysaora melanaster", # jellyfish, showed up in 2003 and slightly increasingly prevalent; more gradual in ai; but not mentioned as of variable detectabilty, and actually rated as a 1 (best detectability) in every year
			
			# species that have variable detectability according to data providers
			# "Bathyraja parmifera", # I already spotted this one
			# "Lepidopsetta polyxystra", # I already spotted this one
			"Acantholumpenus mackayi", "Atheresthes evermanni", "Atheresthes stomias", 
			"Bathyraja aleutica", "Bathyraja interrupta", "Bathyraja maculata", 
			 "Bathyraja taranetzi", "Boreogadus saida", 
			"Careproctus rastrinus", "Ctenodiscus crispatus", "Elassochirus cavimanus", 
			"Elassochirus tenuimanus", "Eurymen gyrinus", "Gymnocanthus detrisus", 
			"Gymnocanthus galeatus", "Gymnocanthus pistilliger", "Gymnocanthus tricuspis", 
			"Hemilepidotus papilio", "Isopsetta isolepis", "Labidochirus splendescens", 
			"Lepidopsetta bilineata", "Leptychaster anomalus", 
			"Lumpenella longirostris", "Lumpenus fabricii", "Lumpenus medius", 
			"Lumpenus sagitta", "Myoxocephalus jaok", "Myoxocephalus polyacanthocephalus", 
			"Myoxocephalus verrucosus", "Octopus dofleini", "Pagurus aleuticus", 
			"Pagurus beringanus", "Pagurus brandti", "Pagurus capillatus", 
			"Pagurus confragosus", "Pagurus cornutus", "Pagurus dalli", "Pagurus kennerlyi", 
			"Pagurus ochotensis", "Pagurus rathbuni", "Pagurus tanneri", 
			"Pagurus trigonocheirus", "Podothecus veternus", "Poroclinus rothrocki", 
			"Raja binoculata", "Sebastes aleutianus", "Sebastes melanostictus", 
			"Stichaeus punctatus"
			
		)
		X.t <- X.t[!spp%in%bad_spp_ebs,]
	}
	
	bad_spp_ai_goa <- c("Lethasterias nanimensis", "Elassochirus tenuimanus", "Sebastes variabilis", "Sebastes ciliatus", "Lepidopsetta bilineata", "Lepidopsetta polyxystra")
	
	if(reg == "ai"){
		man_remove_ai <- c("Aphrocallistes vastus", "Aplidium soldatovi", "Ceramaster japonicus", "Ceramaster patagonicus", "Chlamys albida", "Chrysaora melanaster", "Clathrina blanca", "Crossaster papposus", "Cucumaria fallax", "Diplopteraster multipes", "Elassochirus cavimanus", "Fanellia compressa", "Geodia lendenfeldi", "Halichondria panicea", "Halichondria sitiens", "Halocynthia aurantium", "Henricia aspera", "Henricia leviuscula", "Hippasteria phrygiana", "Homaxinella amphispicula", "Isodictya rigida", "Lebbeus groenlandicus", "Leptasterias coei", "Muriceides nigra", "Mycale loveni", "Ophiura sarsii", "Pagurus brandti", "Pagurus confragosus", "Pagurus trigonocheirus", "Plakina tanaga", "Pododesmus macrochisma", "Porella compressa", "Pseudarchaster parelii", "Pteraster marsippus", "Pteraster tesselatus", "Sebastes melanostictus", "Stegophiura ponderosa", "Strongylocentrotus polyacanthus", "Styela rustica", "Suberites ficus", "")
		bad_spp_ai <- c(
			man_remove_ai,
			bad_spp_ai_goa,
			"Paragorgia arborea", # coral, shows up in 3rd year at ~15% strata
			"Lepidopsetta polyxystra", # flatfish, shows up in 1997 for 80% strata (also in ebs)
			"Pteraster militaris", # seastar, show sup at 40% in 1994
			"Ophiopholis aculeata" # sea star, shows up in 1994 at ~35%
		)
		
		X.t <- X.t[!spp%in%bad_spp_ai,]
	}
	if(reg == "goa"){
		man_remove_goa <- c("Actinauge verrilli", "Alcyonidium pedunculatum", "Aphrodita negligens", "Asterias amurensis", "Buccinum plectrum", "Buccinum scalariforme", "Ceramaster japonicus", "Ceramaster patagonicus", "Cheiraster dawsoni", "Chrysaora melanaster", "Clathrina blanca", "Crossaster borealis", "Crossaster papposus", "Diplopteraster multipes", "Dipsacaster borealis", "Elassochirus cavimanus", "Elassochirus gilli", "Evasterias troscheli", "Fanellia compressa", "Halichondria panicea", "Halichondria sitiens", "Halipteris willemoesi", "Halocynthia dumosa", "Halocynthia igaboja", "Henricia leviuscula", "Hyas lyratus", "Laqueus californianus", "Metridium farcimen", "Mycale loveni", "Myxilla incrustans", "Neoesperiopsis infundibula", "Neptunea amianta", "Ophiopholis aculeata", "Ophiura sarsii", "Orthasterias koehleri", "Pagurus capillatus", "Pagurus confragosus", "Pagurus kennerlyi", "Pagurus trigonocheirus", "Phacellophora camtschatica", "Pododesmus macrochisma", "Pseudarchaster parelii", "Pseudostichopus mollis", "Pteraster militaris", "Pteraster tesselatus", "Ptilosarcus gurneyi", "Pyrulofusus harpa", "Sebastes melanostictus", "Solaster dawsoni", "Solaster endeca", "Stegophiura ponderosa", "Strongylocentrotus polyacanthus", "Styela rustica", "Suberites domuncula", "Suberites ficus", "Synallactes challengeri", "Terebratalia transversa", "")
		bad_spp_goa <- c(
			man_remove_goa,
			bad_spp_ai_goa,
			"Sebastes aleutianus", # shows up at 45% in 2007 and stays high (previously grouped w/ "Sebastes melanostictus")
			"Lepidopsetta bilineata", # rock sole, shows up at 60% in 1996 and stays pretty high
			"Lepidopsetta polyxystra" # rock sole/ flatfish that shows up at ~50% in 1996 and stays high
		)
		X.t <- X.t[!spp%in%bad_spp_goa,]
	}
	
	if(reg == "wctri"){
		man_remove_wctri <- c("Apostichopus leukothele", "Brisaster latifrons", "Lepidopsetta bilineata", "Liponema brevicorne", "Metridium farcimen", "Sardinops sagax", "Stylasterias forreri", "Tritonia diomedea", "")
		bad_spp_wctri <- c(
			man_remove_wctri,
			""
		) 
		X.t <- X.t[!spp%in%bad_spp_wctri,]
	}
	
	if(reg == "wcann"){
		man_remove_wcann <- c("")
		bad_spp_wcann <- c(
			man_remove_wcann,
			""
		) 
		X.t <- X.t[!spp%in%bad_spp_wcann,]
	}
	
	if(reg == "gmex"){
		man_remove_gmex <- c("Astropecten cingulatus", "Astropecten duplicatus", "Aurelia aurita", "Etropus cyclosquamus", "Luidia clathrata", "Ogcocephalus declivirostris", "Pareques iwamotoi", "Pitar cordatus", "Reilla mulleri", "")
		bad_spp_gmex <- c(
			man_remove_gmex,
			"Ophiolepis elegans" # "from 1982-1988, Ophiolepis elegans probably was identified only to the order level of Ophiuroidea" -- Jeff Rester, email to rbatt and mpinsky 2016-04-18
		) 
		X.t <- X.t[!spp%in%bad_spp_gmex,]
	}
	
	
	if(reg == "sa"){
		bad_spp_sa <- c(
			"Limulus polyphemus", # horseshoe crab weren't previously ID'd
			"Libinia dubia", # disappeared b/c they started aggregating to genus
			"Libinia emarginata", # aggreagting to genus
			"Anchoa hepsetus", # aggregating to genus; removed beforehand in trawlData clean.trimRow
			"Anchoa lyolepis", # aggregating to genus; removed beforehand in trawlData clean.trimRow
			"Anchoa mitchilli", # aggregating to genus; removed beforehand in trawlData clean.trimRow
			"Anchoa cubana", # aggregating to genus; removed beforehand in trawlData clean.trimRow
			# "Sciaenops ocellatus", # long line red drum
			"Stomolophus meleagris" # cannonball jelly; started ID'ing late
		)
		X.t <- X.t[!spp%in%bad_spp_sa,]
	}
	if(reg == "newf"){
		bad_spp_newf <- c("Thalarctos maritimus")
		X.t <- X.t[!spp%in%bad_spp_newf,]
	}
	if(reg == "shelf"){
		bad_spp_shelf <- c("Strongylocentrotus droebachiensis","Cucumaria frondosa") # see email from Don Clark on May 3, 2016
		bad_ref_shelf <- X.t[CODE>=1E3, una(ref)] # see email from Don Clark on May 3, 2016
		X.t <- X.t[!ref%in%bad_ref_shelf,]
		X.t <- X.t[!spp%in%bad_spp_shelf,]
		X.t[,CODE:=NULL] # delete this column
	}
	
	if(cull_show_up){
		o_reg <- reg
		bad_spp <- show_up_spp[(reg)==o_reg,una(spp)]
		X.t <- X.t[!spp%in%bad_spp,]
	}
	
	
	# manually removed species (eye-balled and based on confidence guides)
	
	
	
	# ---- Deal with consistent strata ----
	strat_tol <- X.t[,floor(lu(year)*tolFraction)]
	check_strat(X.t, reg, gridSize=gridSize, strat_tol=strat_tol, append_keep_strat=TRUE, plot=plot)
	
	X.t <- X.t[(keep_strat)]
	X.t[,keep_strat:=NULL]
	
	# get the names of the species that were observed at least 5 times
	lots_spp <- X.t[,apply(table(spp, stratum,year)>1, 1, sum)]>=10
	names_lots_spp <- names(lots_spp[lots_spp])
	
	# drop taxa that aren't species or weren't observed at least 10 times
	X.t <- X.t[(spp%in%names_lots_spp) & (taxLvl=="species" | taxLvl=="subspecies")]
	
	# X.t[,hist(rowSums(table(stratum, year)>0))]
	# X.t[,plot(colSums(table(stratum, year)>0), xlab="year", ylab='nstrat', type='o')]
	# X.t[, lu(stratum)]
	# X.t[, lu(spp)]

	# aggregate to make sure a haul
	# doesn't have duplicates of individuals
	# this aggregation uses the biological sum
	
	if(reg== 'neus'){
		X.t <- trawlAgg(
			X = X.t,
			bioFun = una,
			envFun = meanna,
			bio_lvl = "ref",
			space_lvl = "haulid",
			time_lvl = "haulid",
			bioCols= c("wtcpue"),
			envCols = c("stemp","btemp","depth","lon","lat"),
			metaCols = c("reg","datetime","season","year","lon","lat","stratum","common","spp","species"),
			meta.action = "unique1"
		)
	}
	
	X.agg1 <- trawlAgg(
		X = X.t,
		bioFun = sumna,
		envFun = meanna,
		bio_lvl = "spp",
		space_lvl = "haulid",
		time_lvl = "haulid",
		bioCols= c("wtcpue"),
		envCols = c("stemp","btemp","depth","lon","lat"),
		metaCols = c("reg","datetime","season","year","lon","lat","stratum","common","species"),
		meta.action = "unique1"
	)
	# X.agg1[nAgg!=1,j={par(mar=c(0,0,3,0),mfrow=rbLib::auto.mfrow(lu(spp))); for(i in 1:lu(spp)){trawlData::sppImg(unique(spp)[i],unique(common)[i])}}]

	# drop new "time" column
	# change name of agg to indicate step
	# add a 1 and 1/10 degree grid
	X.agg1[,time_lvl:=NULL]
	setnames(X.agg1, "nAgg", "nAgg1")
	# X.agg1[,aggStrat:=ll2strat(lon, lat, gridSize=0.1)]

	# ========================
	# = Give Abundance Value =
	# ========================
	# X.agg1[,abund:=as.integer(wtcpue>0)]
	# X.agg1[,abund:=as.integer(!is.na(wtcpue))]
	X.agg1[,abund:=1]
	
	
	# ============
	# = Define K =
	# ============
	X.agg1[,K:=as.integer(as.factor(haulid)), by=c("reg","year","stratum")]
	X.agg1[,Kmax:=max(K), by=c("reg","year","stratum")]
	
	return(X.agg1)
	
}