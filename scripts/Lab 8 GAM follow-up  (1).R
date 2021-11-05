##################################################

# Module 6: Generalized Additive Model
#             with negative binomial distribution!

# This script file follows up on Lab 8 to model
# the time trend in Chum salmon catches in the 
# Nenana River (data courtesy of Parker Bradley)

# by Franz Mueter
# Last modified: 10-7-2021

##################################################


# Import and examine data:
Nenana <- read.csv("Nenana fish.csv")
head(Nenana)
# Convert dates to Julian day for plotting and modeling 
Nenana$Julian <- strptime(Nenana$Date,"%m/%d/%Y")$yday
# Select Chum salmon and elimiate late season (nu chum caught!)
Chum <- Nenana[Nenana$Julian <= 190, c(1:5,8)]
# Choose a single sampling site on the right side of the river (R):
plot(Chum ~ Julian, data=Chum, subset = Side=="R", col=2)
Chum.R <- Chum[Chum$Side == "R",] 

# In lab 9, we determined that a negative binomial model provided the best fit:
library(MASS)
fit <- glm.nb(Chum ~ Julian + offset(log(Volume)), data = Chum.R)
summary(fit)

# Plot fitted model and confidence band
plot(Chum/Volume * mean(Volume) ~ Julian, data=Chum, col=2)
# Compute predicted values based on the average volume filtered:
dat <- data.frame(Julian = 130:190, Volume = mean(Chum.R$Vol))
Fitted <- predict(fit, newdata=dat, se=T)
mu <- exp(Fitted$fit)  # predicted values on response scale
lower <-  exp(Fitted$fit - 1.96 * Fitted$se.fit)
upper <-  exp(Fitted$fit + 1.96 * Fitted$se.fit)
x <- dat$Julian
lines(x, mu)
lines(x, upper, lty=2)
lines(x, lower, lty=2)

# Examine residuals over time:
plot(Chum.R$Julian, resid(fit))
abline(h=0, lty=2)
# There is a clear pattern in the residuals that is not
# captured by the simple exponential decline. 

# A logical next step would be to consider covariates
# that may explain the periods of higher or lower
# abundances. 

# If there is no obvious explanation (and a suitable covariate),
# we have two options: 
#  1. Model the mean catches as a function of Julian day using 
#     a more flexible approach that allows for "peaks" and "valleys"!
#  2. Fit the exponential decline as before and account for the 
#     pattern of variability around the fitted line by including
#     autocorrelation in the residuals

# Here we explore the first option, because our goal is to explore the
#  pattern in mean counts over time!

# You just learned about flexible, smooth functions that might do
# the job in the 'Generalized Additive Modeling' module and we 
# will apply a GAM with an appropriate error distribution here:

library(mgcv)
# (Unlike 'glm') 'gam' provides a negative binomial family for estimating
#  overdispersed counts. We can implemented this as follows:
fit.gam <- gam(Chum ~ s(Julian) + offset(log(Volume)), data = Chum.R, family = negbin(theta = 1.6))
# Note that this fixes theta at 1.6, close to the value that we estimated before (see 'fit')
summary(fit.gam)
# Fitted model on log-scale:
plot(fit.gam, shade=T, resid=T, pch=16)
plot(Chum.R$Julian, resid(fit.gam), type="b")

# Note that there is still some evidence of a systematic pattern early 
#  in the season, hence we may want to increase the degrees of freedom 
#  of the smooth fit to allow for an even more flexible fit:
fit.gam <- gam(Chum ~ s(Julian, k=15) + offset(log(Volume)), data = Chum.R, family = negbin(theta = 1.6))
summary(fit.gam)
plot(fit.gam, shade=T, resid=T, pch=16)

# Residual diagnostics:
par(mfrow=c(2,2))
gam.check(fit.gam)
par(mfrow=c(1,1))


# In this case, our goal is to capture as much of the higher-frequency variability 
# (on the order of several days) as possible, henc it is OK to increase the degrees
# of freedom to allow for more flexibility in the smooth fit. 

# However, in most situations, specifically if the non-linear function reflects a
# functional relationship that captures some biological process (e.g. physiological, 
# density-dependent, predator-prey, age-dependent, etc), I strongly encourage you 
# to consider limiting the degree of smoothing to no more than 3 or 4 equivalent 
# degrees of freedom, which can capture a variety of biologically reasonable 
# relationships such as linear, dome-shaped, asymptotic, sigmoidal, etc 

# To limit the smooth to 'p' degrees of freedom, set 'k' in the smoother to 'p+1'
# e.g.:   gam(y ~ s(x, k=4)) will use no more than 3 d.f. for the smooth function


# Plot fitted model and confidence band on the back-transformed (original) scale:
plot(Chum/Volume * mean(Volume) ~ Julian, data=Chum, col=2, type="n",
       xlab = "Day of Year", ylab = "Number of fish per set", cex.lab=1.6)
# Compute predicted values based on the average volume filtered:
dat <- data.frame(Julian = 130:190, Volume = mean(Chum.R$Vol))
Fitted <- predict(fit.gam, newdata=dat, se=T)
mu <- exp(Fitted$fit)  # predicted values on response scale
lower <-  exp(Fitted$fit - 1.96 * Fitted$se.fit)
upper <-  exp(Fitted$fit + 1.96 * Fitted$se.fit)
x <- dat$Julian
polygon(c(x,rev(x)), c(upper, rev(lower)), col="lightgrey", lty=0)
lines(x, mu); lines(x, upper, lty=2); lines(x, lower, lty=2)
points(Chum/Volume * mean(Volume) ~ Julian, data=Chum, col=4)


## Note that we fixed theta above, which is probably adequate as we already
#  had an estimate from the 'glm.nb' fit. However, we can also estimate 
#  theta, and there are two approaches to doing so:
#  1. We can use the 'negbin'  family with two values that are the endpoints
#    of an interval over which to search for theta:
fit.gam2 <- gam(Chum ~ s(Julian, k=15) + offset(log(Volume)), data = Chum.R, 
                family = negbin(theta = c(1,10)), optimizer = "perf")
summary(fit.gam2)
plot(fit.gam2, shade=T, resid=T, pch=16)

# The second option is to use the 'nb' family to estimate theta:
fit.gam3 <- gam(Chum ~ s(Julian, k=15) + offset(log(Volume)), data = Chum.R, 
                family = nb())
summary(fit.gam3)
plot(fit.gam3, shade=T, resid=T, pch=16)

# The results differ considerably and the AIC criterion suggests that the
# smaller theta value (fit.gam2) provides a better fit with about the
# same number of parameters (df):
AIC(fit.gam, fit.gam2, fit.gam3)

# Residual diagnostics for AIC-best model:
par(mfrow=c(2,2))
gam.check(fit.gam2)
par(mfrow=c(1,1))
# No major concerns revealed by the diagnostic plots
# and there is no longer any obvious temporal 
# autocorrelation based on Durbin-Watson test:
plot(Chum.R$Julian, resid(fit.gam2), type="b")
library(lmtest)
dwtest(resid(fit.gam2)~1)

# Plot final fitted model and confidence band on the back-transformed (original) scale:
plot(Chum/Volume * mean(Volume) ~ Julian, data=Chum, col=2, type="n",
     xlab = "Day of Year", ylab = "Number of fish per set", cex.lab=1.6)
dat <- data.frame(Julian = 130:190, Volume = mean(Chum.R$Vol))
Fitted <- predict(fit.gam, newdata=dat, se=T)
mu <- exp(Fitted$fit)  # predicted values on response scale
lower <-  exp(Fitted$fit - 1.96 * Fitted$se.fit)
upper <-  exp(Fitted$fit + 1.96 * Fitted$se.fit)
x <- dat$Julian
polygon(c(x,rev(x)), c(upper, rev(lower)), col="lightgrey", lty=0)
lines(x, mu); lines(x, upper, lty=2); lines(x, lower, lty=2)
points(Chum/Volume * mean(Volume) ~ Julian, data=Chum, col=4)

# Conclusions:
# There is an overall significant decrease in the number of chum
# salmon caught from mid-May through early August. Peak counts
# occur early in the season (suggesting that the actual peak
# may have occurred even earlier), with two minor peaks around
# Julian Day 145 (late May) and 168 (mid-June). It is unclear
# what causes these pulses, but that could be explored through
# including appropriate covariates (e.g. precipitation or 
# discharge at upstream locations at some appropriate lag)


