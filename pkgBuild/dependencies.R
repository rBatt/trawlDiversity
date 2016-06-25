
update_dependencies <- function(){
	devtools::use_package("data.table", type="Depends") # Basis for handling all data sets
	devtools::use_package("trawlData", type="Depends") # Meets basic requirements for data content and format
	
	devtools::use_package("rbLib", type="Imports")
	
	devtools::use_package("R2jags", type="Suggests")
	
	devtools::use_package("rstan", type="Suggests")
	
	devtools::use_package("reshape2", type="Suggests")
	
	devtools::use_package("vegan", type="Suggests")
	
	devtools::use_package("ade4", type="Suggests")
	
	devtools::use_package("zoo", type="Suggests")
	
	devtools::use_package("spatstat", type="Suggests")
	
	devtools::use_package("spdep", type="Suggests")
	
	devtools::use_package("raster", type="Suggests")
	
	devtools::use_package("Kendall", type="Suggests")
	
	devtools::use_package("forecast", type="Suggests")
	
}



