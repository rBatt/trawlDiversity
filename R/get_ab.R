#' Get Alpha and Beta
#' 
#' Get alpha and beta parameters given output of an MSOM model run in Stan or JAGS.
#' 
#' @param iD inputData element from output of \code{\link{run_msom}}
#' @param o out element from output of \code{\link{run_msom}}
#' @param yr the year corresponding to input
#' 
#' @details
#' Gets the alpha and beta parameters from the msom. Alpha parameters are related to true presence of species, beta parameters are related to detection.
#' 
#' @return
#' A data.table with iterations of all alpha and beta parameters, which are species-specific. Columns for chain, parameter, value (value is a sample of the posterior), par (alpha or beta), ab_ind (which alpha/ beta parameter, which is the first dimension in parameter), spp_id (species index, second dimension in parameter), spp (species name), and year.
#' 
#' @export
get_ab <- function(iD, o, yr){
	
	# ---- put together spp dimension of parameter names ----
	cnx <- colnames(iD$X)
	cnx_spp <- !grepl("Unknown_[0-9]*", cnx)
	t_spp <- cnx[cnx_spp]
	spp_brack <- paste0(which(cnx_spp), "]")
	n_spp <- length(spp_brack)
	
	# ---- put together alpha and beta dimension of parameter names ----
	sl <- o$BUGSoutput$sims.list
	n_alpha <- ncol(sl$alpha_mu)
	n_beta <- ncol(sl$beta_mu)
	dnames <- dimnames(o$BUGSoutput$sims.array)[[3]]
	
	# alpha_brack <- paste0("alpha[",1:n_alpha,",")
	use_alpha_comma <- grepl("alpha[[0-9]*,.*]", dnames[grepl("alpha",dnames)])
	if(n_alpha==1 & !any(use_alpha_comma)){
		alpha_brack <- "alpha["
	}else{
		alpha_brack <- paste0("alpha[",1:n_alpha,",")
	}
	
	# beta_brack <- paste0("beta[",1:n_beta, ",")
	use_beta_comma <- grepl("beta[[0-9]*,.*]", dnames[grepl("beta",dnames)])
	if(n_beta==1 & !any(use_beta_comma)){
		beta_brack <- "beta["
	}else{
		beta_brack <- paste0("beta[",1:n_beta, ",")
	}
	
	# ---- alpha and beta parameter names to grab ----
	alpha_params <- paste0(rep(alpha_brack, each=n_spp), rep(spp_brack, n_alpha))
	beta_params <- paste0(rep(beta_brack, each=n_spp), rep(spp_brack, n_beta))
	
	# ---- keep track of what's what ----
	alpha_id <- rep(1:n_alpha, each=n_spp)
	beta_id <- rep(1:n_beta, each=n_spp)
	spp_alpha_id <- rep(1:n_spp, n_alpha)
	spp_beta_id <- rep(1:n_spp, n_beta)
	
	# ---- get iters ----
	gi <- get_iters(X=o, pars=c(alpha_params, beta_params), lang="JAGS", FALSE)
	n_iter <- nrow(gi)
	
	# ---- output data.table ----
	ab <- data.table:::melt.data.table(gi, id.vars="chain", variable.name="parameter", "value.name"="value")
	
	a_name <- rep("alpha", length(alpha_id)*n_iter)
	b_name <- rep("beta", length(beta_id)*n_iter)
	ab[,par:=c(a_name, b_name)] # indicate whether each value is an alpha or beta
	ab[,ab_ind:=c(rep(alpha_id, each=n_iter), rep(beta_id, each=n_iter))] # indicate which alpha/ beta

	ab[,spp_id:=rep(c(spp_alpha_id,spp_beta_id), each=n_iter)] # spp index
	ab[,spp:=t_spp[(spp_id)]] # spp name
	ab[,year:=as.integer(yr)]
	
	# ---- return ----
	return(ab)
}
