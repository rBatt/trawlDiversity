
# ==============================
# = Simulating Spreading Range =
# ==============================
# ---- setup ----
nspp <- 300
nsite <- 300
nr <- nsite
nc <- nspp
n_iter <- 300

# ---- transfer matrix ----
transf <- diag(nrow=nsite)
transf[row(transf)-col(transf)==1] <- 1
transf[row(transf)-col(transf)==-1] <- 1

# ---- start with each species in different site ----
starting <- matrix(0, nrow=nsite, ncol=nspp)
start_occ <- sample(1:nsite, nspp) # the sites that start with a species
starting[as.matrix(data.frame(row=start_occ, col=1:nspp))] <- 1

# ---- dynamics ----
occupancy <- vector('list', n_iter)
occupancy[[1]] <- starting
for(i in 2:n_iter){
	occupancy[[i]] <- pmin(transf%*%occupancy[[i-1]],1)
}

# =================================
# = Calculate Range and Diversity =
# =================================
# ---- range size ----
rangeSize <- lapply(occupancy, colSums)
rangeSize_mu <- sapply(rangeSize, mean)

# ---- alpha diversity ----
alphaDiv <- lapply(occupancy, rowSums)
alphaDiv_mu <- sapply(alphaDiv, mean)

# ---- Euclidean beta ----
euc <- function(x){
	sum(scale(x,scale=FALSE)^2)/(nrow(x)-1)
}
betaDiv <- sapply(occupancy, euc) #

# ---- Jaccard beta ----
jaccard <- function(x){
	n <- nrow(x)
	d <- ade4::dist.binary(x, method=1) # requires ade4 package
	sstot <- sum(d^2)/n
	sstot/(n-1)
}
betaDiv_jaccard <- sapply(occupancy, jaccard)

# ---- Hellinger beta ----
hellinger <- function(x){
	n <- nrow(x)
	Y <- vegan::decostand(x, "hellinger") # requires vegan package
	s <- scale(Y, center=TRUE, scale=FALSE)^2
	sum(s) /(n-1)
}
betaDiv_hellinger <- sapply(occupancy, hellinger)


# ========
# = Plot =
# ========
# pdf("~/Desktop/range_alpha_beta.pdf", width=3.5, height=6)
dev.new(width=3.5, height=6)
par(mfrow=c(4,1), mar=c(1.75,1.75,0.5,0.5), tcl=-0.15, mgp=c(0.85,0.2,0), ps=8, cex=1)

plot(1:n_iter, rangeSize_mu, type='l')

plot(1:n_iter, alphaDiv_mu, type='l')

plot(1:n_iter, betaDiv, type='l')
par(new=TRUE)
plot(1:n_iter, betaDiv_jaccard, type='l', xaxt='n', yaxt='n', xlab='', ylab='', col='blue')
axis(side=4)
par(new=TRUE)
plot(1:n_iter, betaDiv_hellinger, type='l', xaxt='n', yaxt='n', xlab='', ylab='', col='red')
legend(
	"topright", 
	legend=c("Euclidean", "Jaccard", "Hellinger"), 
	col=c("black","blue","red"), 
	lty=1, y.intersp = 0.75, bty='n'
)

plot(rangeSize_mu, betaDiv, type='l')
par(new=TRUE)
plot(rangeSize_mu, betaDiv_jaccard, type='l', xaxt='n', yaxt='n', xlab='', ylab='', col='blue')
axis(side=4)
par(new=TRUE)
plot(rangeSize_mu, betaDiv_hellinger, type='l', xaxt='n', yaxt='n', xlab='', ylab='', col='red')
# dev.off()






