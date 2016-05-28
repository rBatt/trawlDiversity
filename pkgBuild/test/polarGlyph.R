

# ---- start function ----
polarGlyph <- function(x, y, ctr_pt, rad_range=c(0,2*pi), g.cex=1, col="black", lwd=1, lcol="black"){

	xo <- order(x)
	x <- x[xo]
	y <- y[xo]

	x <- x/max(x)
	x <- x - min(x)

	y <- y/max(y)
	y <- y - min(y)



	r <- y/max(y)#/10 #sqrt(x^2 + y^2)
	# phi <- x/max(x)*2*pi #atan2(y, x)
	rad_seq <- do.call('seq', c(as.list(rad_range),length.out=list(length(x))))
	phi <- x/max(x)*diff(range(rad_seq)) + rad_seq[1]

	yr <- diff(par()$usr[3:4])
	xr <- diff(par()$usr[1:2])
	cxy <- mean(par("cin"))/par("pin")*c(xr,yr)*g.cex #par("cxy") # par("cin")/par("pin")
	cx <- cxy[1]
	cy <- cxy[2] #[2]

	pol_x <- r * cos(phi) * cx
	pol_y <- r * sin(phi) * cy

	len <- length(pol_x)
	xc <- ctr_pt[1]
	yc <- ctr_pt[2]

	# pseq <- seq(0,2*pi, length.out=100)
	pseq <- do.call('seq', c(as.list(rad_range),length.out=list(100)))

	# plot(pol_x, pol_y, ylim=c(-1,1)+ctr_pt[2], xlim=c(-1,1)+ctr_pt[1],type="n")
	polygon(c(xc+pol_x, rep(xc,len)), c(yc+pol_y, rep(yc, len)), col=col, border=NA)
	lines(cos(pseq)*cx+xc, sin(pseq)*cy+yc, type='l', lwd=lwd, col=lcol)
}


# ---- example inputs ----
x <- sample(1:75, 50) #rnorm(50)
y <- cumsum(rnorm(50))[rank(x)] + x*0.05
ctr_pt <- c(x[30], y[30])
rad_range <- c(0.5*pi,pi*1.5) #c(1.5*pi,pi*2.5) #c(0, pi)
g.cex <- 0.75
col <- adjustcolor("blue", 0.5)

dev.new()
plot(x[order(x)],y[order(x)],type='l')
polarGlyph(x, y, ctr_pt, rad_range=c(0.5*pi,pi*1.5), g.cex=1, col=adjustcolor("red", 0.5), lwd=1, lcol="red")
polarGlyph(x, y, ctr_pt, rad_range=c(1.5*pi,pi*2.5), g.cex=1, col=adjustcolor("blue", 0.5), lwd=1, lcol='blue')


# points(xc, yc)
# lines(cos(phi)*cx/5+xc, sin(phi)*cy/5+yc, type='l', col='green')

# points(pol_x[1],pol_y[1])


