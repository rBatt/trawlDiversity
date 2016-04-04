

set.seed(1)
n <- 20
mu.x1 <- 2
sd.x1 <- 2
mu.x2 <- 3
sd.x2 <- 5

a0 <- 0.5
a1 <- 0.1
a2 <- -0.8
a3 <- 1.2
a4 <- -0.1

ox1 <- rnorm(n=n, mean=mu.x1, sd=sd.x1)
ox1.2 <- ox1^2
ox2 <- rnorm(n=n, mean=mu.x2, sd=sd.x2)
ox2.2 <- ox2^2

ox1.mu <- mean(ox1)
ox1.sd <- sd(ox1)
ox2.mu <- mean(ox2)
ox2.sd <- sd(ox2)
mu.vec <- c(ox1.mu, ox1.mu^2, ox2.mu, ox2.mu^2)
sd.vec <- c(ox1.sd, ox1.sd^2, ox2.sd, ox2.sd^2)


# ox1.mu <- mean(ox1)
# ox1.2.mu <- mean(ox1.2)
# ox1.sd <- sd(ox1)
# ox1.2.sd <- sd(ox1.2)
# ox2.mu <- mean(ox2)
# ox2.2.mu <- mean(ox2.2)
# ox2.sd <- sd(ox2)
# ox2.2.sd <- sd(ox2.2)
# mu.vec <- c(ox1.mu, ox1.2.mu, ox2.mu, ox2.2.mu)
# sd.vec <- c(ox1.sd, ox1.2.sd, ox2.sd, ox2.2.sd)

sx1 <- scale(ox1)
sx1.2 <- sx1^2
sx2 <- scale(ox2)
sx2.2 <- sx2^2

eps <- rnorm(n, mean=0, sd=0.05)
y <- a0 + a1*ox1 + a2*ox1.2 + a3*ox2 + a4*ox2.2 + eps

omod <- lm(y~ox1+ox1.2+ox2+ox2.2)
smod <- lm(y~sx1+sx1.2+sx2+sx2.2)

summary(omod)
summary(smod)

o_ahat <- coef(omod)
s_ahat <- coef(smod)

cor(fitted(omod), fitted(smod))

# ---- Banner ----
rescale.coefs <- function(beta,mu,sigma) {
    beta2 <- beta ## inherit names etc.
    beta2[-1] <- sigma[1]*beta[-1]/sigma[-1]
    beta2[1]  <- sigma[1]*beta[1]+mu[1]-sum(beta2[-1]*mu[-1])
    beta2
}


o_ahat
beta <- rescale.coefs(s_ahat, c(0,mu.vec), c(1,sd.vec))
y2 <- beta[1] + beta[2]*ox1 + beta[3]*ox1.2 + beta[4]*ox2 + beta[5]*ox2.2 +eps


plot(y, y2)
abline(a=0, b=1)

