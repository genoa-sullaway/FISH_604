############################
#
# FISH 604
# Module 2: Basics (part 1)
# Franz Mueter
# 
############################

#### Measures of location and spread:
# 1. Arithmetic, geometric, and harmonic means (slide 6)

# Random draws from normal distribution with mean = 0, sd = 1
n <- 100 # sample size
x <- rnorm(n, mean=10)  # Draw random numbers, mean = 10, sd = 1
hist(x)

AM <- mean(x)
abline(v=AM, col=4, lwd=3)

GM <- prod(x)^(1/n)
abline(v=GM, col=3, lwd=3)

HM <- 1/mean(1/x)
abline(v=HM, col=2, lwd=3)

# Moderately skewed distribution
n <- 1000
x <- rnorm(n, 1, 0.2)^3   
hist(x)

AM <- mean(x)
abline(v=AM, col=4, lwd=3)

GM <- prod(x)^(1/n)
abline(v=GM, col=3, lwd=3)

HM <- 1/mean(1/x)
abline(v=HM, col=2, lwd=3)

# Note that geometric mean is equal to antilog of the arithmetic mean of log(x):
GM
exp(mean(log(x)))

# Note that HM < GM < AM (always!)


### Median and mean:
x <- rnorm(n, 1, 0.2)^3   # right-skewed
hist(x)
abline(v=mean(x), col=4, lwd=3)
abline(v=median(x), col=2, lwd=3)
# Median < Mean for right-skewed data

x <- rnorm(n, 1, 0.2)^0.3   # left-skewed
hist(x)
abline(v=mean(x), col=4, lwd=3)
abline(v=median(x), col=2, lwd=3)
# Mean < Median for left-skewed data


###################################################################
# Example figures:
# 1. Dot plot
x <- rlnorm(25)
dotchart(x, paste("Obs.", 1:25), cex=1.2, pch=16, col=4)
title("Distribution of log-normal variable, n = 25", cex.main=1.6)

# Histograms
par(mfrow=c(2,2))
hist(rnorm(100), nclass=12, cex=1.4)
hist(rlnorm(100), nclass=12, cex=1.4)
x1 <- rnorm(250, 25, 4)
x2 <- rnorm(180, 38, 6)
hist(c(x1,x2), nclass=20, cex=1.4)
hist(c(rnorm(100), 6.1), nclass=12, cex=1.4)

# Example standard error plot (simulated data)
# (see function 'error.bar', appended below)
par(mfrow=c(1,1))
error.bar(1982:2008, rnorm(27, 10), lower = rep(0.8, 27), 
  upper=rep(2, 27), pch=16, cex=1.5, col=4, xlab="")

# Example boxplot (simulated random normal data):
fac <- factor(rep(1:10, each = 40))
x <- rnorm(400)
boxplot(x~fac, col=3)
abline(h=0, lty=2)


# See 'Mod 2 (1) variance.R' for examples of estimating variance 
# (delta method, bootstrap)




############################################################################
# Function to plot points with standard errors:
# (See also function 'errbar' in 'Hmisc' package and 'geom_errorbar" in ggplot2)

error.bar <- function(x, y = NULL, lower, upper, incr = T, bar.ends = T, gap = T, add = F, horizontal = F, ...,
xlab = deparse(substitute(x)), xlim, ylim)
{
draw.null.warn <- function(draw, gap)
{
if(!any(draw)) {
warning("Not enough room for a gap.")
draw <- !draw
gap <- 0
}
invisible(list(draw = draw, gap = gap))
}
if(missing(x))
stop("no data for x or y")
if(missing(y)) {
if(missing(xlab))
xlab <- "Index"
y <- x
x <- time(x)
}
n <- length(x)
if(length(y) != n)
stop("length of y must equal the length of x")
center <- if(horizontal) x else y
if(missing(lower))
stop("you must provide lower")
if(length(lower) > 1 && length(lower) != n)
stop("length of lower must be 1 or equal to the length of x")
#if incr=T lower is assumed >=0
if(incr) lower <- center - abs(lower) else lower <- rep(lower, length = n)
if(any(lower >= center))
warning(paste("There are values of 'lower' which are greater or equal to ", if(
horizontal) "x" else "y"))
if(missing(upper))
upper <- 2 * center - lower
else {
if(length(upper) > 1 && length(upper) != n)
stop("length of upper must be 1 or equal to the length of x")
if(incr)
upper <- center + upper
else upper <- rep(upper, length = n)
}
if(any(upper <= center))
warning(paste("There are values of 'upper' which are smaller or\nequal to ", if(
horizontal) "x" else "y"))
if(!add)
if(horizontal) {
if(missing(ylim))
plot(x, y, xlim = if(missing(xlim)) range(c(lower, upper), na.rm = T)
 else xlim, xlab = xlab, ...)
else plot(x, y, xlim = if(missing(xlim)) range(c(lower, upper), na.rm = T)
 else xlim, ylim = ylim, xlab = xlab, ...)
}
else {
if(missing(xlim))
plot(x, y, ylim = if(missing(ylim)) range(c(lower, upper), na.rm = T)
 else ylim, xlab = xlab, ...)
else plot(x, y, ylim = if(missing(ylim)) range(c(lower, upper), na.rm = T)
 else ylim, xlim = xlim, xlab = xlab, ...)
}
if(horizontal) {
if(gap)
gap <- 0.75 * par("cxy")[1]
draw <- x - lower > gap
z <- draw.null.warn(draw, gap)
draw <- z$draw
gap <- z$gap
segments(lower[draw], y[draw], x[draw] - gap, y[draw])
draw <- upper - x > gap
z <- draw.null.warn(draw, gap)
draw <- z$draw
gap <- z$gap
segments(x[draw] + gap, y[draw], upper[draw], y[draw])
if(bar.ends) {
size.bar <- par("cxy")[2]
segments(lower, y - size.bar, lower, y + size.bar)
segments(upper, y - size.bar, upper, y + size.bar)
}
}
else {
if(gap)
gap <- 0.75 * par("cxy")[2]
draw <- upper - y > gap
z <- draw.null.warn(draw, gap)
draw <- z$draw
gap <- z$gap
segments(x[draw], y[draw] + gap, x[draw], upper[draw])
draw <- y - lower > gap
z <- draw.null.warn(draw, gap)
draw <- z$draw
gap <- z$gap
segments(x[draw], y[draw] - gap, x[draw], lower[draw])
if(bar.ends) {
size.bar <- par("cxy")[1]/4
segments(x - size.bar, upper, x + size.bar, upper)
segments(x - size.bar, lower, x + size.bar, lower)
}
}
}


