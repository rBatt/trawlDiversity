library('maps')
library('raster')
library('spatstat')
library("spdep")
library('rbLib')
library('trawlDiversity')
library("data.table")


# ===================
# = Model Summary 1 =
# ===================
# This is a function that'll help summarize model fit and coeffs and parameter significance:  
mod_smry <- function(m, pred_name=c("density","size","time","type","time:type")){
	pred_name <- match.arg(pred_name)
	sc <- sem.coefs(m)
	sc[,c("estimate","p.value")] <- lapply(sc[,c("estimate","p.value")], signif, 4)
	mod_call <- switch(class(m), lmerMod=m@call, lm=m$call)
	mod_call <- as.character(mod_call)[2]
	fits <- sem.model.fits(m)
	fits[,c("Marginal","Conditional")] <- lapply(fits[,c("Marginal","Conditional")], signif, 4)
	out <- cbind(
		mod_call = mod_call,
		sc[sc[,"predictor"]==pred_name,],
		fits, 
		AIC=round(AIC(m), getOption("digits"))
	)
	out[,c(
		# "std.error","N",
		"Class","mod_call","predictor","estimate","p.value","Marginal","Conditional","AIC"
	)]
}


# ===================
# = Model Summary 2 =
# ===================
# summary function used in tables
mod_smry2 <- function(m){
	mod_call <- switch(class(m), lmerMod=m@call, lm=m$call)
	mod_call <- as.character(mod_call)[2]
	fits <- tryCatch(sem.model.fits(m), error=function(cond)NA)
	if(all(is.na(fits))){
		warning("Error in sem.model.fits")
		fits <- data.frame(Class=class(m), Family="gaussian", Link="identity", N=as.numeric(nobs(m)),"Marginal"=NA_real_,"Conditional"=NA_real_) # family always gaussian for lmer and lm
	}else{
		fits[,c("Marginal","Conditional")] <- lapply(fits[,c("Marginal","Conditional")], signif, 4)
	}
	anova_sum <- car::Anova(m)
	
	out <- data.frame(
		mod_call = mod_call,
		data.frame(predictor=rownames(anova_sum), anova_sum)[,c("predictor","Pr..Chisq.")],
		fits, 
		AIC=AIC(m)#signif(AIC(m), getOption("digits"))
	)
	out[,c("Class","mod_call", "predictor", "Pr..Chisq.","Marginal","Conditional","AIC")]
}


# =============
# = Get Coefs =
# =============
getCoefs <- function(mList){
	outList <- list()
	for(r in 1:length(ur)){
		outList[[ur[r]]] <- rbindlist(lapply(
			coef(mList[[r]]), function(x)as.list(colMeans(x))
		), idcol=TRUE)
		setnames(outList[[ur[r]]], c(".id"), c("randGrp"))
		mod_call <- switch(class(mList[[r]]), lmerMod=mList[[r]]@call, lm=mList[[r]]$call)
		mod_call <- as.character(mod_call)[2]
		outList[[ur[r]]][,mod_call:=mod_call]
	}
	outList <- rbindlist(outList, idcol=TRUE)
	setnames(outList, c(".id"), c("reg"))
	# setcolorder(outList, c("reg", "mod_call", "randomGroup", "(Intercept)", "time", "typepre_ext"))
	outList
}


# ============================
# = Summarize List of Models =
# ============================
# Wrapper for Applying Model Summary to list of Models
smry_modList <- function(ml, pred_name="size"){
	rbindlist(lapply(ml, mod_smry, pred_name="size"))
}

smry_modList2 <- function(ml){
	o <- rbindlist(lapply(ml, mod_smry2), idcol=TRUE) # reg=names(ml), 
	setnames(o, old=c(".id","Pr..Chisq.","Marginal","Conditional"), new=c('reg',"pval","MargR2","CondR2"))
	setkey(o, reg, Class, mod_call, predictor)
	o[, BH:=round(p.adjust(pval, "BH"), 3), by=c("mod_call", "predictor")]# adjust p-values for multiple tests
	o <- dcast(o, reg+Class+mod_call+MargR2+CondR2+AIC~predictor, value.var=c("pval", "BH"))# rearrange so each model on 1 line
	
	# combine fit stats with pvals and coefficients
	o <- merge(o, getCoefs(ml), by=c('reg','mod_call'))
	
	# shorten width of output
	setnames(o, old=c("(Intercept)"), new=c("Int")) # shorter names
	o[,mod_call:=gsub(" ", "", mod_call)] # remove spaces to shorter
	o[,Class:=NULL] # all the same when using mod_smry2
	return(o)
}


# ==================================================
# = Get Estimate, P-value, R Squared from LM Model =
# ==================================================
getEPR <- function(x){
	col_select <- c("Estimate","Pr(>|t|)")
	sx <- summary(x)
	EP <- sx$coefficients[2,col_select]
	names(EP) <- c("Estimate","Pr")
	R <- sx$r.squared
	data.table(as.data.table(as.list(EP)), Rsquared=R)
}



