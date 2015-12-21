#' Print Names of MSOM Parameters
#' 
#' Given a specified language and model type (as in \code{\link{run_msom}}), and the dimensions of the data (e.g., \code{\link{msomData}}), list the potential parameter names as well as the total number of parameters
#' 
#' @param language The language of the MSOM model (JAGS or Stan)
#' @param model_type The type of model (Dynamic or Static)
#' @param nT the number of years
#' @param Jmax the total number of unique sites ever visited
#' @param nU the number of presence covariates
#' @param nV the number of detection covariates
#' @param nS the number of species (including 0's)
#' @param noIndex logical, default FALSE; if TRUE, strip indices from parameter names and return unique values (this can return a much smaller object, and greatly reduce computation time)
#' 
#' @details
#' Currently does not include reporting covariates specified as random variables, nor do these extra parameters factor into the total parameter count. Z is only available in JAGS, as is w. JAGS is parameterized with tau's, Stan with sd's. Only dynamic models include parameters for persistence (phi) and colonization (gamma). Even though they aren't strictly parameters, all models and languages can report parameters for probability of presence (Psi) and probability of detection (Theta). Most hyperparameters are available reporting; however, parameters used strictly to generate non-centered distrubtions (used in Stan only) are not available. In addition to those already mentioned, parameters available include the following (and the pertinent variations on these parameters): alpha, beta, Omega.
#' 
#' Parameter names are semi-redundant: e.g., the 'w' parameter in a JAGS model is species-specific. In the output list of all possible parameters, 'w' will be listed, as will 'w[1]', 'w[2]', ..., 'w[nS]'. The interpretation here is that 'w' implies all of its indices, but that both languages accept both specifications, so both versions are included. Note that this introduces a small asymmetry between the lengths reported in \code{n} and the lengths in the other elements of the output list (see Examples).
#' 
#' @return
#' A named list with the following elements:  
#' \tabular{ll}{
#' n \tab a named vector with integer elements (see below)\cr
#' params \tab the names of all potential parameters\cr
#' params_main \tab the names of all main-effect parameters\cr
#' params_random \tab the names of all random-effect parameters\cr
#' params_latent \tab the names of all latent paramters\cr
#' }
#' 
#' \code{n} contains the following elements
#' \tabular{ll}{
#' nP \tab total parameters\cr
#' n_main \tab total main-effect parameters, including hyperparameters\cr
#' n_random \tab number of random-effect parameters\cr
#' n_latent \tab number of random-effect parameters\cr
#' }
#' 
#' @examples
#' lang <- "JAGS"
#' mt <- "Static"
#' out <- msom_params(lang, mt, nT=3, Jmax=4, nU=3, nV=1, nS=5)
#' 
#' # In reference to the asymetry mentioned in Details
#' out[["n"]][2:4] == sapply(out[3:5], length) # symmetrical
#' out[["n"]][1] < length(out[["params"]]) # asymmetric
#' 
#' @export
msom_params <- function(language, model_type, nT, Jmax, nU, nV, nS, noIndex=FALSE){
	
	if(!noIndex){
		p_ind_s <- paste0("[", 1:nS, "]")
		p_ind_u <- paste0("[", 1:nU, "]")
		p_ind_v <- paste0("[", 1:nV, "]")
		
		p_grid_tj <- expand.grid(1:nT, 1:Jmax)
		p_ind_tj <- paste0("[", paste(p_grid_tj[,1], p_grid_tj[,2], sep=","), "]")
		
		p_grid_us <- expand.grid(1:nU, 1:nS)
		p_ind_us <- paste0("[", paste(p_grid_us[,1], p_grid_us[,2], sep=","), "]")
		
		p_grid_vs <- expand.grid(1:nV, 1:nS)
		p_ind_vs <- paste0("[", paste(p_grid_vs[,1], p_grid_vs[,2], sep=","), "]")
		
		p_grid_tjs <- expand.grid(1:nT, 1:Jmax, 1:nS)
		p_ind_tjs <- paste0("[", paste(p_grid_tjs[,1], p_grid_tjs[,2], p_grid_tjs[,3], sep=","), "]")
	}
	
	if(language=="JAGS"){
		
		params <- c(
			"Omega", "w",
			"alpha_mu", "alpha_tau", 
			"beta_mu", "beta_tau",
			"alpha", "beta", "Z"
			# , "Psi", "Theta"
		)
		if(model_type == "Dynamic"){
			params <- c(params, "phi_mu_logit", "phi_tau_logit", "gamma_mu_logit", "gamma_tau_logit")
		}
		
		if(!noIndex){
			params <- c(
				params,
				paste0("alpha_mu",p_ind_u), "alpha_tau", paste0("alpha_tau",p_ind_u),
				paste0("beta_mu",p_ind_v), "beta_tau", paste0("beta_tau",p_ind_v),
				paste0("w",p_ind_s), paste0("alpha", p_ind_us), paste0("beta", p_ind_vs)
				# , paste0("Psi",p_ind_tjs), paste0("Theta", p_ind_tjs), paste0("Z", p_ind_tjs)
			)
			
			if(model_type == "Dynamic"){
				params <- c(params, paste0("phi", p_ind_s), paste0("gamma", p_ind_s))
			}
			
		}
		
	}else if (language=="Stan"){

		params <- c(
			"Omega",
			"alpha_mu", "alpha_sd", 
			"beta_mu", "beta_sd",
			"alpha", "beta", "logit_psi", "logit_theta"
		)
		if(model_type == "Dynamic"){
			params <- c(params, "phi_mu_logit", "phi_sd_logit", "gamma_mu_logit", "gamma_sd_logit")
		}
		
		if(!noIndex){
			params <- c(
				params,
				paste0("alpha_mu",p_ind_u), "alpha_sd", paste0("alpha_sd",p_ind_u),
				paste0("beta_mu",p_ind_v), "beta_sd", paste0("beta_sd",p_ind_v),
				paste0("w",p_ind_s), paste0("alpha", p_ind_us), paste0("beta", p_ind_vs), 
				paste0("logit_psi",p_ind_tjs), paste0("logit_theta", p_ind_tjs)
			)
			
			if(model_type == "Dynamic"){
				params <- c(params, paste0("phi", p_ind_s), paste0("gamma", p_ind_s))
			}
			
		}
	}
	
	
	params_index_irrelevant_ind <- grepl("Omega|phi_mu|phi_sd|gamma_mu|gamma_sd|phi_tau|gamma_tau", params)
	params_main_ind <- grepl("Omega|(alpha|beta|phi|gamma)_(mu|tau|sd)", params)
	params_latent_ind <- grepl("[Pp]si|[Tt]heta|Z", params)
	params_random_ind <- !params_main_ind & !params_latent_ind & !params_index_irrelevant_ind
	params_just_indices_ind <- grepl("\\[", params)
	

	
	unique_with_ind <- ifInd_strip_noInd(params)
	params_main <- ifInd_strip_noInd(params[params_main_ind])
	params_random <- ifInd_strip_noInd(params[params_random_ind])
	params_latent <- ifInd_strip_noInd(params[params_latent_ind])
	
	# ---- count parameters ----
	n_main <- (1 + nU*2 + nV*2) # Omega, alpha_mu/tau, beta_mu/tau
	n_random <- (nU*nS + nV*nS) # alpha beta
	n_latent <- (nT*Jmax*nS*2) # Psi/ Theta
	if(language=="JAGS"){
		n_random <- n_random + nS # w
		n_latent <- n_latent + nT*Jmax*nS # Z
	}
	if(model_type=="Dynamic"){
		n_main <- n_main + 4 # phi/gamma mu/tau
		n_random <- n_random + 2*nS # phi/gamma
	}
	
	nP <- n_main + n_random + n_latent
	
	out <- list(
		n = c("nP"=nP, "n_main"=n_main, "n_random"=n_random, "n_latent"=n_latent),
		params = params,
		params_main = params_main,
		params_random = params_random,
		params_latent = params_latent
	)
	
	
	return(out)
	
}
