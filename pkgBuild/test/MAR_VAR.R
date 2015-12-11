#' ---
#' title: "Multivariate/ Vector Autoregressive Model with Covariates"
#' author: "Stan users mailing list"
#' date: "December 11, 2015"
#' output: 
#'   pdf_document:
#'     fig_caption: true
#' ---
#'   
#' Multivariate autogressive (MAR) and vector autoregressive (VAR) are the same thing, ecologists call them MAR. Here we fit a MAR(1) model and include covariates. The model is fit using `stats::ar()`, `vars::VAR()`, and a model written in Stan. A VAR model was written in Stan by Rob Trangucci, Ryan Batt added covariates to it and put this document together.  Find related discussions on the Stan users mailing list [here](https://groups.google.com/forum/#!topic/stan-users/PmGMnqXRD28) and [here](https://groups.google.com/forum/#!topic/stan-users/8RerHVzxjUQ).  
#' 
#' The model presented here is framed after [Ives et al. (2003)](https://paperpile.com/shared/XfSD1V); in particular, Eq. 27 (quoting from the paper):   
#'  $$X_t = A + BX_{t-1} + CU_t + E_t$$  
#' "where $U_t$ is a $q \times 1$ vector containing the values of $q$ covariates at time $t$, and $C$ is a $p \times q$ matrix whose elements $c_ij$ give the strength of effect of covariates $j$ on species $i$. The covariates $U_t$ can be any factors that affect the system ..."  
#' 
#' The purpose of this document is to fit a MAR model in Stan. The data are simulated with an ecological process in mind (nutrients, phytoplankton [algae], and zooplankton [herbivores]), but that's just a loose framing and shouldn't be read-into too much.  


#+ load_packages, echo=TRUE
# =================
# = Load Packages =
# =================
library(vars)
library(rstan)


#+ sim_setup, echo=TRUE
# =========
# = Setup =
# =========
set.seed(1337)
nT <- 2E2
nX <- 2 # don't change
nU <- 1 # don't change
e_sd <- 1

A <- matrix(c(0.1, 10), nrow=nX, ncol=1)
B <- matrix(c(0.75, -0.8, 0.2,  0.95), nrow=nX, ncol=nX)
C <- matrix(c(0,1), nrow=nX, ncol=nU)

X <- matrix(NA, nrow=nT, ncol=2, dimnames=list(NULL, c("zoop","phyto")))
U <- matrix(rlnorm(nT*nU), nrow=nT, ncol=nU, dimnames=list(NULL, c("nuts")))
E <- matrix(rnorm(nT*nX) * e_sd, nrow=nT, ncol=nX)



#+ simulate, echo=TRUE
# ============
# = Simulate =
# ============
X[1,] <- c(3,10)

for(i in 2:nT){
	X[i,] <- A + B%*%matrix(X[i-1,],ncol=1) + C%*%matrix(U[i,],ncol=1) + E[i,]
}


#+ time_series_fig, fig.width=5, fig.height=6, fig.cap="Time series of zooplankton (X[,1], black), phytoplankton (X[,2], blue), and nutrients (U, red)", fig.align="center", echo=TRUE
# =================
# = Plot Simulate =
# =================
par(mfrow=c(2,1), mar=c(2,2,0.1,2), cex=1, ps=9, mgp=c(1,0.1,0), tcl=-0.1)
plot(X[,1], 
	type="l", 
	ylim=c(min(X,0), max(X)), 
	ylab="X (zoops [black] phytos [blue])", 
	xlab="time step"
)

lines(X[,2], col="blue")

par(new=TRUE)
plot(U, 
	type="l", 
	col="red", 
	xlab="", ylab="", 
	xaxt="n", yaxt="n", 
	lwd=0.5, 
	lty="dashed"
)
axis(side=4, col="red")
mtext("U (nuts [red])", side=4, line=1)

plot(X[,2:1], type="o", lwd=0.25, pch=20)
#' Nutrients come into the system as a log-normal random variate; this representation loosely represents the distribution of nutrient pulses into a lake that would be observed in nature. This description and the simulation doesn't emphasize the distinction between nutrient loading and nutrient concentration. When there is a spike in the nutrients (red dashed line), phytoplankton jump up (`C[2,1]`). Then in the next time step, the zooplankton benefit from a greater presence of phytoplankton in the previous time step (`B[1,2]` in the simulation, or `B[1,1,2]` in the Stan model output). On the other hand, phytoplankton abundance is intensely suppressed by zooplankton abundance in the previous time step (`B[2,1]` in simulation, or `B[1,2,1]` in Stan output). Thus, a nutrient pulse quickly increases phytoplankton, which are subsequently suppressed by zooplankton. Further, the diagonal elements of `B` are less than 1 and greater than 0, indicating a steady decline in the abundance of both phytoplankton and zooplankton, but the decline is balanced a steady (constant) input from `A`. The `A` parameter is quite small for zooplankton ($A_1 = `r A[1]`$) and much larger for phytoplankton ($A_2 = `r A[2]`$).  
#' 
#' The top panel shows the time series of these state variables, where the spikes in nutrients and plankton are all plainly visible, as are the shifted timings and differences in the relative sizes of the spikes.  
#' 
#' The bottom panel of the figure shows the phase plot of the two state variables ($X$). The effects of large nutrient pulses are seen as the rapid excursions to the right (more phytoplankton), which are followed by the system moving up (more zooplankton) and to the left (fewer phytoplankton), and eventually spiralling counterclockwise back in toward the equilibrium point (which the system may never perfectly reach due to shocks [$E$] and nutrient pulses [$U$] bombarding the system at each time step).  


#+ stan_model, echo=TRUE
# ==============
# = Stan Model =
# ==============
cat(readLines('VAR_p_cov.stan', warn=FALSE), sep="\n")


#+ fit_models, echo=TRUE
# ==============
# = Fit Models =
# ==============
# MAR
ar(X) # can't handle covariates :(

# VAR
VAR(X, exogen=U)

# Stan VAR/ MAR
library(rstan)
stan_data <- list(T = nT, M = nX, P = 1, Y = X, U = U, nU = nU)

model_w_data <- stan(file = 'VAR_p_cov.stan', data=stan_data, iter=200, chains=4, cores=4)

print(model_w_data, pars=c("A","B","C"))

traceplot(model_w_data, inc_warmup=FALSE, pars=c("A","B","C"))


#+ summary, echo=TRUE
# ===========
# = Summary =
# ===========
#' #Summary  
#' All three models estimated the coefficients and intercepts ~correctly (`ar` didn't print out the intercept terms, oddly, but it at least got the coefficients pretty close). Both VAR and the Stan model estimated all parameters quite accurately.