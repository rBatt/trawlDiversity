#' Unpack P
#' 
#' Unpacks the object containing the processed MSOM results
#' 
#' @param p_rn The p object subset to a region
#' 
#' @details
#' Make sure to only supply p after it's been subset (e.g., \code{p[[1]]}). The structure of p is such that the first level is the regions, the second is the years, then the third contains the varius elements associated with the processed output of a single model run.
#' 
#' @export
unpack_p <- function(p_rn){
	rd <<- p_rn$rd
	processed <<- p_rn$processed
	bt <<- p_rn$bt
	colonization <<- p_rn$colonization
	param_iters <<- p_rn$param_iters
	ab <<- p_rn$ab
	naive_rich <<- p_rn$naive_rich

	reg <<- processed[,una(reg)]
	lang <<- "JAGS"

	if(lang == "Stan"){
		pars_trace <<- c("Omega","alpha_mu[1]", "alpha_mu[2]", "alpha_mu[3]", "alpha_mu[4]", "alpha_mu[5]", "beta_mu[1]")
	}else{
		pars_trace <<- c("Omega","alpha_mu[1]", "alpha_mu[2]", "alpha_mu[3]", "alpha_mu[4]", "alpha_mu[5]", "beta_mu")
	}

	naive_rich <<- naive_rich[,naive_rich, keyby='year']
	reg_rich <<- processed[,reg_rich, keyby='year']
	bt_ann <<- bt[,list(bt_ann=mean(bt, na.rm=TRUE)), keyby='year']

	n_pars <<- length(pars_trace)
	n_yrs <<- param_iters[,lu(year)]
	n_spp <<- rd[,lu(spp)]
	
	spp_col <<- colonization[[1]][value==1, una(spp)]
	spp_ext <<- colonization[[1]][value==-1, una(spp)]
	spp_col_only <<- spp_col[!spp_col %in% spp_ext]
	spp_ext_only <<- spp_ext[!spp_ext %in% spp_col]
	spp_col_and_ext <<- spp_col[spp_col%in%spp_ext]
	spp_neither <<- rd[!spp%in%spp_col & !spp%in%spp_ext, una(spp)]
}