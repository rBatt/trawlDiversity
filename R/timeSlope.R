# The fruit of: https://github.com/rBatt/trawl/commit/a04d8e166d7093538c698368d48542a0e60bd16c
timeSlope<- function(x, ...){
	# Old function
	# if(sum(!is.na(x))<3){
	# 	NA
	# }else{
	# 	as.numeric(lm(x~I(0:(length(x)-1)))$coef[2])
	# }
	nona <- !is.na(x)
	if(sum(nona)<3){
		NA
	}else{
		lx <- length(x)
		if(lx>=4E3){
			tvec <- (1:lx)[nona]
			x2 <- x[nona]
			cor(x2, tvec)*(sd(x2)/sd(tvec))

		}else{
			tvec <- cbind(1,1:lx)[nona,] # the 1 is a dummy vector for the intercept
			ttvec <- t(tvec)
			(solve(ttvec%*%tvec)%*%ttvec%*%x[nona])[2]
		}
	}
}