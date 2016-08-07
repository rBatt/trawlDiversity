#' Kendall's Tau with Serial Dependence
#' 
#' Use Kendall's Tau to test for a significant trend in a time series, while accounting for potential serial dependence.
#' 
#' @param x numeric vector of times (dates, step number, etc)
#' @param y numeric vector of value to be inspected for a trend
#' @param tauOnly Logical, if TRUE (FALSE is default) only the Kendall's Tau is returned
#' @param ... arguments passed to \code{\link{cor}}
#' 
#' @details
#' The general approach has 3 components: estimate the linear trend, estimate the serial dependence, then estimate Kendall's Tau. 
#' 
#' The simplest way to do this is to estimate and remove the linear trend from the data, model and remove first order autocorrelation, add the trend back into the AR residuals, then fit Kendall's Tau to the Trend+Residuals series. This is the approach described in the studies in the references section, and is implemented in the \href{https://cran.r-project.org/web/packages/zyp/index.html}{zyp package}.
#' 
#' This function improves upon the above approach by fitting the linear trend and the serial dependence at the same time, and by considering more complex models of serial dependence (specifically, up to ARMA(2,2) models). The fitting of the ARMA model is done using the \code{forecast::auto.arima} function; if one already knows the order of the ARMA model, an equivalent result could be achieved through \code{\link{arima}}. The key to this type of model is that it is fit as a regression with ARMA residuals (see the Hyndman 2004 reference below for more details, or \href{https://www.otexts.org/fpp/9/1}{or this chapter in the free online text}). The order of the model is limited to second order autoregressive (AR(2) component) and second order moving average (MA(2) component). Differencing, drift, and seasonality are not considered. This is partially done for speed, and partially to protect against the model soaking up too much of any nonlinear trend (because the trend we're explicitly modeling is linear, but Kendall's Tau wouldn't be needed if that was the actual trend was linear!).
#' 
#' \code{Kendall::Kendall} is used. The variance of the numerator in Kendall's Tau is more complex when there are ties, and changes the p-value. \code{cor(..., method="kendall")} does not make the correction.
#' 
#' @return
#' If \code{tauOnly} is TRUE, returns a length-1 numeric vector with name \code{"tau"} that is Kendall's Tau. Otherwise, returns a length 8 list:
#' \tabular{ll}{
#' \code{b} \tab the slope of the linear trend, as estimated alongisde the ARMA model \cr
#' \code{bp} \tab P-value associated with the linear trend \cr
#' \code{tau} \tab Kendall's Tau \cr
#' \code{S} \tab score \cr
#' \code{varS} \tab variance of score \cr
#' \code{tauZ} \tab Z statistic associated with Kendall's Tau (either S/sqrt(varS) or qnorm(taup))
#' \code{taup} \tab p-value of tau \cr
#' \code{adjusted} \tab character indicating whether S or Z had to be adjusted to match p-value reported by Kendall::Kendall 
#' }
#' 
#' @references
#' Yue, S., and C. Y. Wang. 2002. Applicability of prewhitening to eliminate the influence of serial correlation on the Mann-Kendall test. Water resources research 38:4-1-4-7.
#' 
#' Yue, S., and C. Wang. 2004. The Mann-Kendall Test Modified by Effective Sample Size to Detect Trend in Serially Correlated Hydrological Series. Water Resources Management 18:201-218.
#' 
#' Hyndman, R. J., and G. Athanasopoulos. 2014. Forecasting: principles and practice: OTexts.
#' 
#' @examples
#' n <- 10
#' x <- 1:n
#' slope <- 0.5
#' ar_coeff <- 0.85
#' y <- slope*x+arima.sim(model=list(ar=c(ar_coeff)), n=n)
#' tsTau(x, y)
#' 
#' @seealso \code{\link{post_trend}}
#' 
#' @export
tsTau <- function(x, y, tauOnly=FALSE, ...){
	requireNamespace("forecast", quietly = TRUE)
	requireNamespace("Kendall", quietly = TRUE)
	
	mod_stat <- forecast::auto.arima(y, xreg=x, d=0, max.p=2, max.q=2, seasonal=FALSE, allowdrift=FALSE, approximation=TRUE)
	b <- unname(coef(mod_stat)['x'])
	y_stat <- residuals(mod_stat) + b*x[!is.na(y)]
	
	mod <- Kendall::Kendall(x, y_stat)
	tau <- mod[['tau']][1]
	
	if(!tauOnly){
		bp_se <- sqrt(vcov(mod_stat)['x','x'])
		bp_tail <- c(FALSE,TRUE)[(sign(b)==-1)+1L]
		bp <- pnorm(b, 0, sd=bp_se, lower.tail=bp_tail)*2
	}
	
	if(tauOnly){
		return(c(tau=tau))
	}else{
		taup <- mod[['sl']][1]
		varS <- mod[['varS']][1]
		S <- mod[['S']][1]
		Z <- qnorm(taup/2, lower.tail=(S<0)) # infer Z score from Kendall::Kendall pvalue
		
		# Check that the output can reproduce the Kendall::Kendall p-value
		# Assumes Kendall::Kendall p-value is "right"
		# First, I calculate a p-value using Kendall:Kendall S and varS (S/varS=Z)
		# A lack of match for that ^ manual p-value might be due to absence of correction in S (continuity)
		# Also, mismatch can be due to rounding error [e.g., very small p], which means I need to recalc Z
		# ---- Check p-value match, due continuity correction if needed ----
		manual_p <- pnorm(S/sqrt(varS), lower.tail=(S<0))*2 # calculate Z score & pvalue from Kendall::Kendall output
		p_check <- isTRUE(all.equal(taup, manual_p, tolerance=1E-6)) # Logical, do they match?
		has_ties <- any(duplicated(y_stat))
		adjusted <- ''
		if(!p_check & has_ties){ # one reason p_check would be FALSE is if there are ties
			S <- S - sign(S) # if p-values don't match, make continuity correction (made if ties)
			adjusted <- c(adjusted, "S")
		}
		# ---- Check that finite Z obtained from qnorm(p), recalc if not ----
		if(!is.finite(Z)){ # if Kendall::Kendall pvalue implies non-finite Z, use Sadj to calculate my own Z
			Z <- (S)/sqrt(varS) # recalc Z, but using Kendall::Kendall S
			adjusted <- c(adjusted, "Z")
		}
		adjusted <- as.character(paste(adjusted, collapse=""))
		return(list(b=b, bp=bp, tau=tau, S=S, varS=varS, tauZ=Z, taup=taup, adjusted=adjusted))
	}
}


#' Trend through Posterior Distributions
#' 
#' Fit a trend to a time series of posterior estimates. Ends up being similar to bootstrapping.
#' 
#' @param x time step
#' @param Y matrix of posterior samples. Each column is a time step, each row is a sample of the posterior
#' @param nSamp number of times each set of posterior samples should be resampled to fit trend. I.e., number of bootstrap iterations.
#' 
#' @details 
#' Uses Kendall's Tau and accounts for serial correlation in the time series (up to ARMA(2,2)). Uses Kendall::Kendall. Accounts for ties. 
#' 
#' @return
#' A named numeric vector of length 3, containing Kendall's Tau, the associated Z statistic, and the p-value
#' 
#' @examples
#' n <- 20
#' x <- 1:n
#' slope <- 0.5
#' n_post <- 50
#' ar_coeff <- 0.85
#' y <- c(slope*x+arima.sim(model=list(ar=c(ar_coeff)), n=n))
#' y_post <- sapply(y, function(x)rnorm(n_post, mean=x, sd=1))
#' post_trend(x, y_post, 10)
#' 
#' @seealso \code{\link{tsTau}}
#' 
#' @export
post_trend <- function(x, Y, nSamp=100){
	
	nPost <- nrow(Y)
	n <- length(x)
	# D <- n*(n-1)/2
	
	post_tauOut <- list() #matrix(0, nrow=nSamp, ncol=1, dimnames=list(NULL, c("tau")))
	blank_tsTau <- function(){
		list(b=NA_real_, bp=NA_real_, tau=NA_real_, S=NA_real_, varS=NA_real_, tauZ=NA_real_, taup=NA_real_, adjusted=NA_character_)
	}
	for(i in 1:nSamp){
		ty <- Y[cbind(sample(seq_len(nPost), n, replace=TRUE), seq_len(n))]
		post_tauOut[[i]] <- tryCatch(tsTau(x, ty, tauOnly=FALSE), error=function(cond)blank_tsTau())
	}
	
	post_tauOut <- rbindlist(post_tauOut)
	post_tauZ <- post_tauOut[!is.na(tau), tau]	
	tau_mu <- mean(post_tauZ, na.rm=TRUE)

	# # S <- post_tauOut[!is.na(tau),S] #post_tauZ*D
	# S <- post_tauOut[!is.na(tau),Sadj]
	# S_mu <- mean(S)
	# S_sd <- post_tauOut[!is.na(tau),sqrt(varS)] #sqrt(D*(2*n+5))/3
	# S_sd_mu <- mean(S_sd)
	# Z <- S_mu/S_sd_mu # DONT DO CONTINUITY CORRECTION HERE: already performed in tsTau
	
	Z <- post_tauOut[!is.na(tau), mean(tauZ, na.rm=TRUE)]
	pvalue <- pnorm(Z, lower.tail=(Z<0))*2
	
	return(c(tau=tau_mu, Z=Z, pvalue=pvalue))
}