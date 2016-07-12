globalVariables(c("V1","spp","stratum","K","abund","btemp","stemp","depth","datetime","yr","doy","thin","lon", "lat", "keep_strat","spp", "haulid", "taxLvl", "time_lvl", "abund", "wtcpue", "K", "Kmax", "A", "sigmaA", "b", "a", "strat_mean_6yr_mean", "strat_mean_9yr_mean", "strat_mean_6yr_sd", "strat_mean_9yr_sd", "ab_ind", "spp_ind", "value", "n_spp_col", "col_logic", "n_strat_col", "n_spp_col_weighted", "chain", "reg", "year_added", "years_elapsed", "years_remaining", "years_seen_after", "bt_col", "extract", "Omega", "iter", "variable", "bt", "nK", "richness", "spp_id", "strat_mean", "a1", "a2", "a3", "a4", "a5", "ext", "CODE", "ab", "bt_ann", "bt_max", "bt_mean", "bt_min", "btemp.mu", "btemp.sd", "common", "depth.mu", "depth.sd", "lang", "latP", "lonP", "n_col", "n_pars", "n_spp_ext_weighted", "n_yrs", "na.pass", "naive_rich", "ref", "reg_rich", "rel_col_ext", "sppID", "spp_col", "spp_col_and_ext", "spp_col_only", "spp_ext", "spp_ext_only", "spp_neither", "status", "unobs_rich", "value_mu", "x", "y", "n"))

# utils::globalVariables(c("rd", "processed", "bt", "colonization", "param_iters", "ab", "naive_rich", "reg", "lang", "reg_rich", "bt_ann", "n_pars", "n_yrs", "n_spp", "spp_col", "spp_ext", "spp_col_only", "spp_ext_only", "spp_col_and_ext", "spp_neither"))
# to_assign <- c("rd", "processed", "bt", "colonization", "param_iters", "ab", "naive_rich", "reg", "lang", "reg_rich", "bt_ann", "n_pars", "n_yrs", "n_spp", "spp_col", "spp_ext", "spp_col_only", "spp_ext_only", "spp_col_and_ext", "spp_neither")
# suppressBindingNotes <- function(variablesMentionedInNotes) {
#     for(variable in variablesMentionedInNotes) {
# 				# if(exists(variable, envir = parent.frame())) {assign(variable,NULL, envir = parent.frame())}else{assign(variable,NULL, envir = .GlobalEnv)}
# 				assign(variable,NULL, envir = .GlobalEnv)
#     }
# }
# suppressBindingNotes(to_assign)