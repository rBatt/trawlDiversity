
update_dependencies <- function(){
	devtools::use_package("data.table", type="Depends") # Basis for handling all data sets
	devtools::use_package("trawlData", type="Suggests") # Meets basic requirements for data content and format
}



