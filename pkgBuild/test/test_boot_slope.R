

# ================================
# = Setup Options for Simulation =
# ================================
n <- 20
x <- 1:n
slope <- 0.5
n_post <- 50
ar_coeff <- 0.85


# ==================================================
# = Simulate a Time Series and Draw from Posterior =
# ==================================================
y_true <- c(slope*x+arima.sim(model=list(ar=c(ar_coeff)), n=n))
y_post <- sapply(y_true, function(x)rnorm(n_post, mean=x, sd=1))
y_mu <- apply(y_post, 2, mean)


# ======================================
# = Calculate Trends on Posterior Mean =
# ======================================
get_trend_lm <- function(x, y){
	mod <- lm(y~x)
	b <- coef(mod)[2]
	p <- summary((mod))$coef[2,4]
	data.frame(method='lm', b=b, tau=NA, bp=p, taup=NA, tauZ=NA)
}
get_trend_Kendall <- function(x, y){
	mod <- Kendall::Kendall(x, y)
	tau <- mod$tau[1]
	p <- mod$sl
	
	# calculate long-way from equations as check
	# n <- length(x)
	S <- mod$S[1]
	# D <- mod$D[1]# n*(n-1)/2 # checks out
	varS <- mod$varS # (sqrt(D*(2*n+5))/3)^2 # checks out
	Z <- S/sqrt(varS) #(3*S)/sqrt(D*(2*n+5))
	
	# p <- pnorm(Z, lower.tail=FALSE)*2 # checks out
	# # p <- pnorm(S, mean=0, sd=sqrt(varS), lower.tail=FALSE)*2 # this also works
	
	data.frame(method='Kendall', b=NA, tau=tau, bp=NA, taup=p, tauZ=Z)
}
get_trend_MannKendall <- function(y){
	mod <- Kendall::MannKendall(y)
	tau <- mod$tau[1]
	p <- mod$sl
	
	S <- mod$S[1]
	varS <- mod$varS # (sqrt(D*(2*n+5))/3)^2 # checks out
	Z <- S/sqrt(varS) #(3*S)/sqrt(D*(2*n+5))
	
	data.frame(method='MannKendall', b=NA, tau=tau, bp=NA, taup=p, tauZ=Z)
}

get_trend_zyp <- function(y, x){
	mod <- zyp::zyp.trend.vector(y, x, method="yuepilon")
	b <- mod['trend']
	tau <- mod['tau']
	taup <- mod['sig']
	
	data.frame(method='zyp', b=b, tau=tau, bp=mod['trendp'], taup=taup, tauZ=NA)
}

get_trend_white <- function(x, y, tauOnly=FALSE, checkTies=FALSE, ...){
	mod_stat <- forecast::auto.arima(y, xreg=x, d=0, max.p=2, max.q=2, seasonal=FALSE, allowdrift=FALSE, approximation=TRUE)
	b <- unname(coef(mod_stat)['x'])
	y_stat <- residuals(mod_stat) + b*x[!is.na(y)]
	
	if(!checkTies & tauOnly){
		tau <- cor(x, ty, method="kendall", ...)
	}else{
		mod <- Kendall(x, y_stat)
		tau <- mod[['tau']][1]
	}
	
	if(!tauOnly){
		bp_se <- sqrt(vcov(mod_stat)['x','x'])
		bp_tail <- c(FALSE,TRUE)[(sign(b)==-1)+1L]
		bp <- pnorm(b, 0, sd=bp_se, lower.tail=bp_tail)*2
	}
	
	if(checkTies){
		if(!tauOnly){
			taup <- mod[['sl']][1]
			varS <- mod$varS
			S <- mod[['S']][1]
			Z <- S/sqrt(varS)
			data.frame(method='white', b=b, tau=tau, bp=bp, taup=taup, tauZ=Z)
		}else{
			varS <- mod$varS
			return(matrix(c(tau, varS), ncol=2, nrow=1, dimnames=list(NULL, c("tau", "varS"))))
		}
	}else{
		if(!tauOnly){
			S <- post_tauZ2[,tau*D]
			S_mu <- mean(S)
			S_sd <- sqrt(D*(2*n+5))/3
			Z <- S_mu/S_sd
			data.frame(method='white', b=b, tau=tau, bp=bp, taup=taup, tauZ=Z)
		}else{
			return(matrix(c(tau), ncol=1, nrow=1, dimnames=list(NULL, c("tau"))))
		}
	}
	
}

get_trend <- function(x, y, method=c("lm","Kendall", "MannKendall","zyp","white")){
	method <- match.arg(method)
	switch(method,
		lm = get_trend_lm(x, y),
		Kendall = get_trend_Kendall(x, y),
		MannKendall = get_trend_MannKendall(y),
		zyp = get_trend_zyp(y, x),
		white = get_trend_white(x, y)
	)
}
methods <- c("lm","Kendall", "MannKendall","zyp","white")

data.table::rbindlist(lapply(methods, get_trend, x=x, y=y_mu), fill=TRUE)



# ======================================
# = Posterior/ Bootstrap Trend and Tau =
# ======================================
post_trend <- function(x, Y, nSamp=100){
	
	nPost <- nrow(Y)
	n <- length(x)
	D <- n*(n-1)/2
	post_tauZ <- matrix(0, nrow=nSamp, ncol=1, dimnames=list(NULL, c("tau")))
	for(i in 1:nSamp){
		ty <- Y[cbind(sample(seq_len(nPost), n, replace=TRUE), seq_len(n))]
		post_tauZ[i,] <- get_trend_white(x, ty, TRUE)
	}
	
	post_tauZ2 <- as.data.table(post_tauZ)[complete.cases(post_tauZ)]
	
	tau_mu <- post_tauZ2[,mean(tau)]
	S <- post_tauZ2[,tau*D]
	S_mu <- mean(S)
	S_sd <- sqrt(D*(2*n+5))/3
	Z <- S_mu/S_sd
	tailS <- c(FALSE,TRUE)[(sign(Z)==-1)+1L]
	pvalue <- pnorm(Z, lower.tail=tailS)*2
	
	return(c(tau=tau, Z=Z, pvalue=pvalue))
}
post_trend(x, y_post)



# ====================
# = Plot Time Series =
# ====================
plot(x, y_mu, ylim=range(y_post), type='n')
apply(y_post, 1, function(z)lines(x, z))


