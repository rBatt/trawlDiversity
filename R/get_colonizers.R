#' Get Colonizers
#' 
#' Identify species colonizing, going extinct, and persisting between years
#' 
#' @param d An input data.table with columns spp, year (may have others, but they aren't used)
#' @param smry Logical, when TRUE (default), the results are summarized in terms of number of species in each category; if FALSE, species-specific categorizations are return
#' 
#' @return
#' A data.table with a column for year and columns relating to colonization, extinction, and persistence. If \code{smry} is TRUE, remaining columns are \code{n_col} (number of colonizations), \code{n_ext} (number of extinctions), and \code{n_pers} (number of species persisting). If \code{smry} is FALSE, columns are year, spp, and value. If value is 1 that species colonized (is present now, wasn't last year), if -1 it went extinct (not present now, but was last year), and if 0 it persisted (present both last year and this year).
#' 
#' @export
get_colonizers <- function(d, smry=TRUE){
	
	tbl <- d[,table(spp, year)]
	col_ext <- apply(tbl>0, 1, diff)
	col_ext_dt <- data.table(reshape2::melt(col_ext, varnames=c("year","spp")), key=c("year","spp"))
	
	n_cep <- col_ext_dt[,list(n_col=.SD[value==1, sum(value)], n_ext=.SD[value==-1, abs(sum(value))], n_pers=.SD[value==0, lu(spp)]), by="year"]
	
	if(smry){
		return(n_cep)
	}else{
		return(col_ext_dt)
	}
}