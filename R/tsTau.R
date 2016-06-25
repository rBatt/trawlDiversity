#' Kendall's Tau with Serial Dependence
#' 
#' Use Kendall's Tau to test for a significant trend in a time series, while accounting for potential serial dependence.
#' 
#' @param x numeric vector of times (dates, step number, etc)
#' @param y numeric vector of value to be inspected for a trend
#' @param tauOnly Logical, if FALSE (default) only the Kendall's Tau is returned
#' @param checkTies Logical, if FALSE (default) \code{\link{cor}} is used to compute Kendall's Tau. Otherwise, \code{Kendall::Kendall} is used.
#' @param ... arguments passed to \code{\link{cor}}
#' 
#' @details
#' The general approach has 3 components: estimate the linear trend, estimate the serial dependence, then estimate Kendall's Tau. 
#' 
#' The simplest way to do this is to estimate and remove the linear trend from the data, model and remove first order autocorrelation, add the trend back into the AR residuals, then fit Kendall's Tau to the Trend+Residuals series. This is the approach described in the studies in the references section, and is implemented in the \href{https://cran.r-project.org/web/packages/zyp/index.html}{zyp package}.
#' 
#' This function improves upon the above approach by fitting the linear trend and the serial dependence at the same time, and by considering more complex models of serial dependence (specificall, up to ARMA(2,2) models). The fitting of the ARMA model is done using the \code{forecast::auto.arima} function; if one already knows the order of the ARMA model, an equivalent result could be achieved through \code{\link{arima}}. The key to this type of model is that it is fit as a regression with ARMA residuals (see the Hyndman 2004 reference below for more details, or \href{https://www.otexts.org/fpp/9/1}{or this chapter in the free online text}). The order of the model is limited to second order autoregressive (AR(2) component) and second order moving average (MA(2) component). Differencing, drift, and seasonality are not considered. This is partially done for speed, and partially to protect against the model soaking up too much of any nonlinear trend (b/c the trend we're explicitly modeling is linear, but Kendall's Tau wouldn't be needed if that was the actual trend was linear!).
#' 
#' When \code{checkTies} is TRUE, \code{Kendall::Kendall} is used instead of \code{cor} to estimate Kendall's Tau. This is because the variance of the numerator in Kendall's Tau is more complex when there are ties, and the former handles this case whereas the latter does not.
#' 
#' @return
#' If \code{tauOnly} is TRUE, returns a length-1 numeric vector with name \code{"tau"} that is Kendall's Tau. Otherwise, returns a data.frame with 1 row and 5 numeric columns:
#' \tabular{ll}{
#' \code{b} \tab the slope of the linear trend, as estimated alongisde the ARMA model \cr
#' \code{tau} \tab Kendall's Tau \cr
#' \code{bp} \tab P-value associated with the linear trend \cr
#' \code{taup} \tab P-value associated with Kendall's Tau \cr
#' \code{tauZ} \tab Z statistic associated with Kendall's Tau \cr
#' }
#' 
#' @references
#' Yue, S., and C. Y. Wang. 2002. Applicability of prewhitening to eliminate the influence of serial correlation on the Mann-Kendall test. Water resources research 38:4–1–4–7.
#' 
#' Yue, S., and C. Wang. 2004. The Mann-Kendall Test Modified by Effective Sample Size to Detect Trend in Serially Correlated Hydrological Series. Water Resources Management 18:201–218.
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
tsTau <- function(x, y, tauOnly=FALSE, checkTies=FALSE, ...){
	requireNamespace("forecast", quietly = TRUE)
	
	mod_stat <- forecast::auto.arima(y, xreg=x, d=0, max.p=2, max.q=2, seasonal=FALSE, allowdrift=FALSE, approximation=TRUE)
	b <- unname(coef(mod_stat)['x'])
	y_stat <- residuals(mod_stat) + b*x[!is.na(y)]
	
	if(checkTies){
		requireNamespace("Kendall", quietly = TRUE)
		mod <- Kendall::Kendall(x, y_stat)
		tau <- mod[['tau']][1]
	}else{
		tau <- cor(x, y_stat, method="kendall", ...)
	}
	
	if(!tauOnly){
		bp_se <- sqrt(vcov(mod_stat)['x','x'])
		bp_tail <- c(FALSE,TRUE)[(sign(b)==-1)+1L]
		bp <- pnorm(b, 0, sd=bp_se, lower.tail=bp_tail)*2
	}
	
	if(checkTies){
		if(tauOnly){
			varS <- mod$varS
			return(c(tau=tau))
		}else{
			taup <- mod[['sl']][1]
			varS <- mod$varS
			S <- mod[['S']][1]
			Z <- S/sqrt(varS)
			return(data.frame(b=b, tau=tau, bp=bp, taup=taup, tauZ=Z))
		}
	}else{
		if(tauOnly){
			return(c(tau=tau))
		}else{
			D <- n*(n-1)/2
			S <- tau*D
			S_sd <- sqrt(D*(2*n+5))/3
			Z <- S/S_sd
			Z_tail <- c(FALSE,TRUE)[(sign(b)==-1)+1L]
			taup <- pnorm(Z,lower.tail=Z_tail)*2
			return(data.frame(b=b, tau=tau, bp=bp, taup=taup, tauZ=Z))
		}
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
#' Uses Kendall's Tau and accounts for serial correlation in the time series (up to ARMA(2,2)). Currently assumes that there are no ties in the resampled Y. If there are, the p-value may be incorrect.
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
	D <- n*(n-1)/2
	post_tauZ <- matrix(0, nrow=nSamp, ncol=1, dimnames=list(NULL, c("tau")))
	for(i in 1:nSamp){
		ty <- Y[cbind(sample(seq_len(nPost), n, replace=TRUE), seq_len(n))]
		post_tauZ[i,] <- tsTau(x, ty, tauOnly=TRUE)
	}
	
	post_tauZ2 <- post_tauZ[!is.na(post_tauZ)]
	
	tau_mu <- mean(post_tauZ2)
	S <- post_tauZ2*D
	S_mu <- mean(S)
	S_sd <- sqrt(D*(2*n+5))/3
	Z <- S_mu/S_sd
	tailS <- c(FALSE,TRUE)[(sign(Z)==-1)+1L]
	pvalue <- pnorm(Z, lower.tail=tailS)*2
	
	return(c(tau=tau_mu, Z=Z, pvalue=pvalue))
}