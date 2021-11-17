########################################

# Module 7: Generalized Additive Models
# Franz Mueter

# Last modified: 9-19-2021

########################################


################### LOWESS smooth - choose 'span':
GAK1 <- read.csv("data/GAK1_mod.csv")
str(GAK1)
GAK1$year  # Note that year is in decimal format!

# Select two years of surface data (2002 & 2003, depth 0) for analysis
dat <- GAK1[GAK1$year > 2002 & GAK1$year < 2004 & GAK1$depth == 0,]
dat

attach(dat)  # For simplicity, but generally NOT recommended unless 
             # you know exactly what you are doing!
# This creates a separate copy of 'dat' on your 'search path'. See:
search()
# and allows you to refer to variables in 'dat' by name without 
# specifying the data frame each time

plot(year, temperature)

# LOWESS fit with default level of smoothing (span = 0.75):
fit1 <- loess(temperature ~ year)

# Add fitted values to plot:
x <- seq(2002, 2004, length=100)  # Sequence of x-values for plotting
y1 <- predict(fit1, data.frame(year=x))  # Predicted y-values based on loess fit
lines(x,y1)
# Fit is too smooth and does not capture seasonal peaks! 

# Select a narrower (smaller) bandwidth:
fit2 <- loess(temperature ~ year,  span=0.2) 
y2 <- predict(fit2, data.frame(year=x))
lines(x,y2, col=2)
# This may perhaps not be smooth enough? (overfitting!)

# Try an intermediate bandwidth!
fit3 <- loess(temperature ~ year, span=0.3)
y3 <- predict(fit3, data.frame(year=x))
lines(x,y3, col=4,lwd=2)

# Summary output looks quite different from a 'typical' model:
summary(fit3)
names(fit3)

# Alternative way to specify degree of smoothing via 
# equivalent number of parameters 
# (=trace of smoother matrix):
fit4 <- loess(temperature ~ year, enp.target=12)
y <- predict(fit4, data.frame(year=x))
lines(x,y, col=3,lwd=3)


###############   #### Spline smoothing 
## The following approach provides an automatic choice of the "optimum"
## amount of smoothing (best smoothing parameter) via cross-validation

# The 'smooth.spline' fits a cubic smoothing spline to x and y data,
# the function only takes a single predictor variable

fit.spl <- smooth.spline(year, temperature)    # does NOT take a formula
fit.spl  # Note 'spar' and equivalent df!!
plot(year, temperature)
x <- seq(2002, 2004, length=100)
p <- predict(fit.spl, data.frame(year=x))
lines(p$x[,1], p$y[,1], col=4,lwd=2)
# The default fitting algorithm uses a "Generalized Cross-validation" (GCV) 
# criterion to determine the optimum amount of smoothing. This criterion 
# minimizes the sum of squared residuals adjusted for model complexity 
# (i.e. the "equivalent degrees of freedom" in the model output)

# Note that the generalized cross-validation in this case results in a model 
# that seems to grossly overfit the data (with a negative 'spar' parameter 
# and using 47 degrees of freedom for 47 data points (that is, the data are
# interpolated). This results in some patterns in the fitted values that 
# are not supported by the data. 

# The "leave-one-out" cross-validation option (cv=T) results in a more 
# reasonable fit:
fit.spl2 <- smooth.spline(year, temperature, cv=T)
fit.spl2 # More reasonable values for spar and Df
p2 <- predict(fit.spl2, data.frame(year=x))
lines(p2$x[,1], p2$y[,1], col=3, lwd=2)


# For illustration, for practice (in case you have to), and to understand how
# it works, let's minimize the 'leave-one-out' cross-validation score "by hand" 
# (Something like this is done automatically by 'smooth.spline' if no 'spar' 
# value is specified and if cv=T)

# The algorithm proceeds by leaving out one observation at a time, fitting
# a cubic smoothing spline to all the other observation, predict the value
# for the observation that was left out, and computes the squared residuals. 
# Once every observation has been predicted based on a spline fit to the
# remaining observations, the mean squared prediction error (MSPE) is
# computed as: (observed - predicted)^2

CV.MSPE <- vector("numeric")    # set up vector to hold output
trial.spar <- seq(-1,1,by=0.005) # Vector of trial values for smoothing parameter
trial.spar
trials <-length(trial.spar)
(n <- length(temperature))  # Number of observations

#### WARNING: for loop is slow!!
for(i in 1:trials) {
  CV <- vector("numeric")
  for(k in 1:n) {
    # fit cubic spline to all observations except k 
    fit <- smooth.spline(year[-k], temperature[-k], spar=trial.spar[i])
    # compute predicted value for left-out observation k
    p <- predict(fit, x=year[k])
    # compute squared residual: (observed - predicted)^2
    CV[k] <- (temperature[k] - p$y)^2
    }
    CV.MSPE[i] <- sum(CV)/n  # Mean squared prediction error
  }
plot(trial.spar, CV.MSPE, type="l")
# Notice the "flat" MSPE criterion over a wide range of trial values 
# (from -1 to 0). This may be the cause for the failure of GCV to 
# find the true minimum as the estimated 'spar' value was about -0.8

# Let's 'zoom in' on the region of interest:
plot(trial.spar[250:340], CV.MSPE[250:340], type="l")

# Value of 'spar' that minimizes cross-validation sum of squares:
trial.spar[which.min(CV.MSPE)]

# Compare to the leave-1-out result from 'smooth.spline':
fit.spl2   # Uses different (and much faster) algorithm!


################### Regression splines: 

# Example from Faraway ("Extending the linear model with R"):
# Download package of datasets if needed
if(!("faraway" %in% installed.packages()))
{
  install.packages("faraway")
}
library(faraway)
data(exa)
?exa  # See how data were simulated
plot(y ~ x, data=exa)   # Simulated data
lines(exa$x, exa$m, lwd=2, col=4)     # True means

# Illustrate piecewise linear fit:
spl.lin <- function(x, c) ifelse(x > c, x-c, 0) 
knots <- 0:9/10
knots

# design matrix for regression with knots at these points:
dm <- outer(exa$x, knots, spl.lin)
dim(dm)         # Design matrix of 10 "explanatory variables"
View(dm)        # Examine the design variables
matplot(exa$x, dm, type="l", col=1, xlab="x", ylab="splines", cex.lab=2)

# Compute and show regression fit (i.e. regress y on new linear variables):
fit.lin <- lm(exa$y ~ dm)  # Multiple linear regression on 10 variables
fit.lin
plot(y ~ x, data=exa, pch=16, cex=0.4, xlab = "x", ylab = "y", cex.lab=2)
lines(exa$x, exa$m, lwd=2)     # "True" means
lines(exa$x, fitted(fit.lin), col=2, lwd=2)   # piecewise linear fit

# Illustrate smoothing using cubic B-spline:
library(splines)
x1 <- seq(0,1,length=1000) # Sequence for plotting
spl.cub <- bs(x1, df=12)     # set up cubic B-spline of order 12 (12 knots)
matplot(x1,spl.cub, type="l", col=1, xlab="x", ylab="splines", cex.lab=1.5)

# Set up design matrix of 12 cubic B-splines for our example data:
dm <- bs(exa$x, 12)
dim(dm)  
head(dm)

# Fit linear regression on 12 "explanatory variables":
fit.cbs <- lm(exa$y ~ dm)
fit.cbs
plot(y ~ x, data=exa, pch=16, cex=0.4, xlab = "x", ylab = "y", cex.lab=2)
lines(exa$x, exa$m, lwd=2)      # True means (from simulated data set)
lines(exa$x, fitted(fit.cbs), col=2, lwd=2)   # piecewise linear fit


#################### Fitting a smooth surface in 2 dimensions:
# We use the temperature data at GAK 1 to fit a loess model in two 
# dimensions, that is we fit a smooth, 2-dimensional surface to the
# data for each combination of depth and (decimal) year:

# Select 2003/2004 data, all depths:
j <- GAK1$year > 2002 & GAK1$year < 2004
dat2 <- GAK1[j,]
dat2$depth <- -dat2$depth # reverse sign for plotting

# Fit LOESS model:
fit.lo <- loess(temperature ~ year + depth, data = dat2, span=0.2)

# Set up a grid of values for computing predicted values:
# (For 'image' values in x and y need to be in increasing order)
x <- seq(2002,2004,length=100)
y <- seq(-250, 0, length=100) 

# This creates a grid of all combinations of sampling dates 
# and depths in the x and y vectors:
(grd <- expand.grid(year=x, depth=y))

# Predicted values from model over entire grid:
z <- predict(fit.lo, grd)
library(fields)  # required for the oceanographic color scheme:
image(x, y, z, col = tim.colors(100), 
   xlab="Time", ylab="Depth (m)",
   cex.axis=1.2, cex.lab=1.4)
# show when/where measurements were taken:
points(dat2$year, dat2$depth, pch=16, cex=0.3) 

########## Legend (separate plot)
win.graph(3,1)
par(mar=c(2,.5,0,.5))
xx <- seq(min(z, na.rm=T), max(z,na.rm=T), length=80)
i <- xx[2]-xx[1]
plot(c(xx[1]-i,xx[80]+i),c(0,1),type="n",axes=0,xlab="",ylab="")
rect(xx-i, 0, xx+i, 1, col=tim.colors(80), border=NA)
axis(1, line=0, cex.axis=1)
dev.off()
##########

##### Diagnostics:
# Residuals over time:
plot(dat2$year, resid(fit.lo)); abline(h=0, lty=2)
# Residuals by depth:
plot(dat2$depth, resid(fit.lo)); abline(h=0, lty=2)
# Note higher variance at the surface!
plot(fitted(fit.lo),resid(fit.lo), 				# Evidence of heteroscedasticity!
	xlab="Fitted values", ylab="Residuals")
abline(h=0, col=4); grid()

# Residuals by depth and over time:
r <- resid(fit.lo)
j <- r > 0
image(x, y, z, col = tim.colors(100), 
      xlab="Time", ylab="Depth (m)",
      cex.axis=1.2, cex.lab=1.4)
symbols(dat2$year[j], dat2$depth[j], circles = r[j], 
	fg = 4, bg = 4, inches=0.05, add=T) 
symbols(dat2$year[!j], dat2$depth[!j], circles = abs(r[!j]),
	fg = 2, bg = 2, inches=0.05, add=T) 
# There are some patterns suggesting that temperatures are 
# underestimated at the surface during summer (blue, positive 
# residuals at the highest surface layer temperatures) 


################
# log-transform because of heteroscedasticity:
fit.log <- loess(log(temperature) ~ year + depth, data = dat2, span=0.2)

# Predicted values:
z <- predict(fit.log, grd)
# Plot fitted values on back-transformed scale
# (note bias correction when back-transforming from log-scale)
image(x,y, exp(z+fit.log$s^2/2),col = tim.colors(100), 
   xlab="Time", ylab="Depth (m)",
   cex.axis=1.2, cex.lab=1.4)

# Diagnostics:
plot(dat2$year, resid(fit.log)); abline(h=0, lty=2)
plot(dat2$depth, resid(fit.log)); abline(h=0, lty=2)
plot(fitted(fit.log),resid(fit.log), 				# Evidence of heteroscedasticity!
	xlab="Fitted values", ylab="Residuals")
abline(h=0,lty=2); grid()
# No obvious issues in the residuals

# Residuals by depth and over time:
r <- resid(fit.log)
j <- r > 0
image(x,y, exp(z+fit.log$s^2/2),col = tim.colors(100), 
      xlab="Time", ylab="Depth (m)",
      cex.axis=1.2, cex.lab=1.4)
symbols(dat2$year[j], dat2$depth[j], circles = r[j], 
	fg = 4, bg = 4, inches=0.05, add=T) 
symbols(dat2$year[!j], dat2$depth[!j], circles = abs(r[!j]),
	fg = 2, bg = 2, inches=0.05, add=T) 
# Still some patterns in residuals suggesting non-independence
# Apparent auto-correlation in space (with depth) and over time
# as evident in 'clusters' of positive/negative residuals

detach("dat")  # Detach data frame we attached in the beginning
