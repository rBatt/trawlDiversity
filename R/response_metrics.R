#' Response Metrics for Quadratic Model
#' 
#' For quadratic models predicting responses across an environmental gradients, extract the response and metrics (optimum value, e.g.)
#' 
#' @param alphas a data.table with alpha values; rows are posterior iterations, a1,a2,a3,a4,a5 are columns (see details)
#' @param X a data.table with columns for \code{bt} (bottom temperature) and \code{depth}; bt associated with a2 and a3, depth with a4 and a5
#' @param n_samp integer indicating (random) subsample of posterior to use; if NULL (default), no subsampling is performed
#' @param n_grid number of depth and temperature values over which to estimate response; if NULL (default), response is calculated only for observed combinations in \code{X}
#' @param b0 an intercept parameter value (like a1 in \code{alphas})
#' @param b1 a linear coefficient value (like a2 or a4 in \code{alphas})
#' @param b2 a quadratic coefficient value (like a3 or a5 in \code{alphas})
#' 
#' @details
#' The columns in \code{alphas} are parameters estimated from the MSOM model. The MSOM model has a form of \code{y=a1*1 + a2*bt+a3*bt^2 + a4*depth+a5*depth^2}. Thus, the columns in alpha are interpreted as follows:
#' \tabular{ll}{
#' \code{a1} \tab an intercept term, corresponds to \code{b0} \cr
#' \code{a2} \tab a linear coefficient associated with bottom temeprature, corresponds to \code{b1} \cr
#' \code{a3} \tab a quadratic coefficient associated with bottom temperature, corresponds to \code{b2} \cr
#' \code{a4} \tab a linear coefficient associated with depth, corresponds to \code{b1} \cr
#' \code{a5} \tab a quadratic coefficient assoicated with depth, corresponds to \code{b2} \cr
#' }
#' where \code{b0}, \code{b1}, \code{b2} are arguments to the response metric functions (\code{psi.opt}, \code{psi.tol}, and \code{psi.max}). Note, then, that the response metric functions do not consider models with 2 predictors, but only 1, hence the different notation.
#' 
#' @export get_response psi.opt psi.tol psi.max
#' @name response_metrics
NULL

#' \code{get_response}: Get model predicted response across a range of depth and temperature values
#' @rdname response_metrics
get_response <- function(alphas, X, n_samp=NULL, n_grid=NULL){
	if(!is.null(n_samp)){
		alphas <- alphas[sample(1:nrow(alphas),n_samp)]
		n <- n_samp
	}else{
		n <- nrow(alphas)
	}
	
	if(!is.null(n_grid)){
		d_vals <- X[,seq(min(depth), max(depth), length.out=n_grid)]
		b_vals <- X[,seq(min(bt, na.rm=TRUE), max(bt, na.rm=TRUE), length.out=n_grid)]
		vals <- expand.grid(b=b_vals, d=d_vals)
	}else{
		vals <- as.matrix(X[,list(b=bt,d=depth)])
	}
	
	pred_resp <- function(b, d, alphas){
		a <- matrix(alphas, nrow=5)
		x <- as.matrix(cbind(1, b, b^2, d, d^2))
		x%*%a
	}

	pr <- pred_resp(vals[,"b"], vals[,"d"], t(alphas[,list(a1,a2,a3,a4,a5)]))
	
	if(!is.null(n_grid)){
		dn <- list(btemp=as.character(b_vals),depth=as.character(d_vals),iter=NULL)
		pr2 <- array(pr, dim=c(n_grid, n_grid, n), dimnames=dn)
		pr3 <- apply(pr2, c(1,2,4), mean)
		return(pr3)
	}

	return(data.table(X,resp=rowMeans(pr)))	
}

#' \code{psi.opt}: Optimum environmental value
#' @rdname response_metrics
psi.opt <- function(b1,b2){-b1/(2*b2)}

#' \code{psi.tol}: Tolerance across environmental values
#' @rdname response_metrics
psi.tol <- function(b2){1/sqrt(-2*b2)}

#' \code{psi.max}: Response achieved at environmental optimum
#' @rdname response_metrics
psi.max <- function(b0,b1,b2){1/(1+exp((b1^2)/(4*b2)-b0))}