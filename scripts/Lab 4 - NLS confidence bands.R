# Confidence band for NLS model

# One way to construct confidence bands for non-linear models is through the 
# use of the delta method, which is implemented in the package 'propagate'

library(propagate)

# Based on example in Lab 4 - OLS WLS GLS script

# Import data
dat <- read.csv("Lab 4 OLS-WLS-GLS example.csv")

# fit model: y = a * x^b
fit.OLS <- nls(y ~ a*x^b, data=dat, start = list(a=2.5, b=0.1))
summary(fit.OLS)

# Plot data and fitted model:
(cf <- coef(fit.OLS))
plot(y~x, data=dat, pch=16)
curve(cf[1]*x^cf[2], add=T, lwd=2, col=4)

# Simulate confidence bands:
new.x <- data.frame(x=seq(min(dat$x), max(dat$x), length=20))
pred <- predictNLS(fit.OLS, newdata=new.x)
conf <- pred$summary

# Add confidence bands to plot
lines(conf$"Sim.2.5%" ~ new.x$x, lwd=1, col=4)
lines(conf$"Sim.97.5%" ~ new.x$x, lwd=1, col=4)


