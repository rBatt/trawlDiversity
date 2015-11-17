#' ebs.a
#' 
#' Aggregated data from Eastern Berring Sea, small subset for testing
#' 
#' @format
#' A dim = 184 x 21 data.table data.frame:
#' \tabular{rlll}{
#' [,1] \tab aggStrat \tab character \tab stratum used in spatial aggregation\cr
#' [,2] \tab time_lvl \tab character \tab the temporal factor\cr
#' [,3] \tab spp \tab character \tab taxonomic identifier (species name)\cr
#' [,4] \tab wtcpue \tab numeric \tab weight per unit effort\cr
#' [,5] \tab abund \tab integer \tab number of aggrated hauls with wtcpue >= 0\cr
#' [,6] \tab stemp \tab numeric \tab surface tempeature\cr
#' [,7] \tab btemp \tab numeric \tab bottom temperature\cr
#' [,8] \tab depth \tab numeric \tab depth\cr
#' [,9] \tab lon \tab numeric \tab longitude\cr
#' [,10] \tab lat \tab numeric \tab latitude\cr
#' [,11] \tab reg \tab character \tab region name\cr
#' [,12] \tab datetime \tab c("POSIXct", "POSIXt") \tab date and time\cr
#' [,13] \tab season \tab character \tab season of year\cr
#' [,14] \tab year \tab character \tab year\cr
#' [,15] \tab stratum \tab character \tab stratum, site; based on 1 degree grid\cr
#' [,16] \tab common \tab character \tab species common name\cr
#' [,17] \tab species \tab character \tab species name\cr
#' [,18] \tab nAgg2 \tab integer \tab number of aggregations involved for this species-place-time\cr
#' [,19] \tab haulid \tab character \tab most specific space-time identifier\cr
#' [,20] \tab K \tab integer \tab replicate visits to a stratum during the time_lvl period\cr
#' [,21] \tab Kmax \tab integer \tab maximum number of visits across strata during a given year\cr
#' }
#' @seealso \code{\link{ebs.agg2}}
#' @import data.table
#' @export
"ebs.a"

#' ebs.agg2
#' 
#' #' Aggregated data from Eastern Berring Sea
#' 
#' @format
#' A dim = 254654 x 21 data.table data.frame:
#' \tabular{rlll}{
#' [,1] \tab aggStrat \tab character \tab insert_description_here\cr
#' [,2] \tab time_lvl \tab character \tab insert_description_here\cr
#' [,3] \tab spp \tab character \tab insert_description_here\cr
#' [,4] \tab wtcpue \tab numeric \tab insert_description_here\cr
#' [,5] \tab abund \tab integer \tab insert_description_here\cr
#' [,6] \tab stemp \tab numeric \tab insert_description_here\cr
#' [,7] \tab btemp \tab numeric \tab insert_description_here\cr
#' [,8] \tab depth \tab numeric \tab insert_description_here\cr
#' [,9] \tab lon \tab numeric \tab insert_description_here\cr
#' [,10] \tab lat \tab numeric \tab insert_description_here\cr
#' [,11] \tab reg \tab character \tab insert_description_here\cr
#' [,12] \tab datetime \tab c("POSIXct", "POSIXt") \tab insert_description_here\cr
#' [,13] \tab season \tab character \tab insert_description_here\cr
#' [,14] \tab year \tab character \tab insert_description_here\cr
#' [,15] \tab stratum \tab character \tab insert_description_here\cr
#' [,16] \tab common \tab character \tab insert_description_here\cr
#' [,17] \tab species \tab character \tab insert_description_here\cr
#' [,18] \tab nAgg2 \tab integer \tab insert_description_here\cr
#' [,19] \tab haulid \tab character \tab insert_description_here\cr
#' [,20] \tab K \tab integer \tab insert_description_here\cr
#' [,21] \tab Kmax \tab integer \tab insert_description_here\cr
#' }
#' @seealso \code{\link{ebs.a}}
#' @import data.table
#' @export
"ebs.agg2"