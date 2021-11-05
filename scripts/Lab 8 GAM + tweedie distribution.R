##################################################

# Time trend in Chum salmon catches in the 
# Nenana River (data courtesy of Parker Bradley)

# by Franz Mueter
# Last modified: 10-07-2021

##################################################

# Import and examine data:
Nenana <- read.csv("data/Nenana fish.csv")
head(Nenana)
# Convert dates to Julian day for plotting and modeling 
Nenana$Julian <- strptime(Nenana$Date,"%m/%d/%Y")$yday
# Select Chum salmon and elimiate late season (nu chum caught!)
Chum <- Nenana[Nenana$Julian <= 190, c(1:5,8)]
# Choose a single sampling site on the right side of the river (R):
plot(Chum ~ Julian, data=Chum, subset = Side=="R", col=2)
Chum.R <- Chum[Chum$Side == "R",] 

# I previously determined that a negative binomial model with an 
# offset to account for volume filtered provided a better fit
# than Poisson model or zero-inflated linear models:
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
plot(Chum.R$Julian, resid(fit), type="b"); abline(h=0, col=4)

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

# Apply a GAM with an appropriate error distribution here:

library(mgcv)
# Using 'gam' with a negative binomial family for estimating overdispersed counts:
fit.gam <- gam(Chum ~ s(Julian) + offset(log(Volume)), data = Chum.R, family = negbin(theta = 1.6))
# Note that this fixes theta at 1.6, close to the value that we estimated before (see 'fit')
summary(fit.gam)
plot(fit.gam, shade=T, resid=T, pch=16)
plot(Chum.R$Julian, resid(fit.gam), type="b"); abline(h=0, col=4)

# Note that there is still some evidence of a systematic pattern early 
#  in the season, hence we increase the degrees of freedom of the smooth
#  fit to allow for an even more flexible fit:
fit.gam <- gam(Chum ~ s(Julian, k=15) + offset(log(Volume)), data = Chum.R, family = negbin(theta = 1.6))
summary(fit.gam)
plot(fit.gam, shade=T, resid=T, pch=16)

# In this case, our goal is to capture as much of the higher-frequency variability 
# (on the order of several days) as possible, hence it is OK to increase the degrees
# of freedom to allow for more flexibility in the smooth fit. 

# Plot fitted model and confidence band on the back-transformed (original) scale:
plot(Chum/Volume * mean(Volume) ~ Julian, data=Chum.R, col=2, type="n",
     xlab = "Day of Year", ylab = "Number of fish per set", 
     cex.lab=1.6, ylim=c(0,60))
# Compute predicted values based on the average volume filtered:
dat <- data.frame(Julian = 130:190, Volume = mean(Chum.R$Vol))
Fitted <- predict(fit.gam, newdata=dat, se=T)
mu <- exp(Fitted$fit)  # predicted values on response scale
lower <-  exp(Fitted$fit - 1.96 * Fitted$se.fit)
upper <-  exp(Fitted$fit + 1.96 * Fitted$se.fit)
x <- dat$Julian
polygon(c(x,rev(x)), c(upper, rev(lower)), col="lightgrey", lty=0)
lines(x, mu); lines(x, upper, lty=2); lines(x, lower, lty=2)
points(Chum/Volume * mean(Volume) ~ Julian, data=Chum.R, col=4)

# To facilitate plotting the output for different models, it pays to
# to write a quick function to do so:
plot.nenana <- function(Model, Data=Chum.R) {
  Model <- fit.gam2
  # 'Model' should be a model with an appropriate 'predict' method
  plot(Chum/Volume * mean(Volume) ~ Julian, data=Data, col=2, type="n",
       xlab = "Day of Year", ylab = "Number of fish per set", 
       cex.lab=1.6)
  dat <- data.frame(Julian = 130:190, Volume = mean(Data$Vol))
  Fitted <- predict(Model, newdata=dat, se=T)
  mu <- exp(Fitted$fit)  # predicted values on response scale
  lower <-  exp(Fitted$fit - 1.96 * Fitted$se.fit)
  upper <-  exp(Fitted$fit + 1.96 * Fitted$se.fit)
  x <- dat$Julian
  polygon(c(x,rev(x)), c(upper, rev(lower)), col="lightgrey", lty=0)
  lines(x, mu); lines(x, upper, lty=2); lines(x, lower, lty=2)
  points(Chum/Volume * mean(Volume) ~ Julian, data=Data, col=4)
}  



## We can estimate theta, and there are two approaches to doing so:
# 1. We can use the 'negbin'  family with two values that are the endpoints
#   of an interval over which to search for theta:
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
# larger theta value (fit.gam2) provides a better fit with about the same
# number of equivalent parameters:
AIC(fit.gam, fit.gam2, fit.gam3)


##################################################################
# The Tweedie distribution

# A distribution that generalizes both the Poisson and negative 
# binomial, which has gained prominence in modeling count data in 
# recent years is the so-called 'Tweedie' distribution. It estimates
# an additional parameter called p, which ranges between 1 
# (equivalent to the Poisson distribution) and 2 (equivalent 
# to the negative binomial)

# Compare results with a GAM based on the tweedie distribution:
# (the default for the tweedie distribution also uses a log-link)
fit.gam4 <- gam(Chum ~ s(Julian, k=15) + offset(log(Volume)), data = Chum.R, 
                family = tw())
summary(fit.gam4)
# Note that the power parameter (p) for the Tweedie distribution 
# was estimated to be 1.214, which is closer to the Poisson than
# to the negative binomial, perhaps suggesting that the negative
# binomial over-estimates the increase in variance with the mean
# (as it assumes an increase in variance that is proportional to
#  the square of the mean counts)

AIC(fit.gam3, fit.gam4)
# This suggests that the negative binomial model provides a better
# fit to the data. The GAM smooth has about the same number of 
# parameters (edf=9.70 vs. edf=9.96 for the NB model)

# Plot fitted model and confidence band on the back-transformed 
# (original) scale:
plot.nenana(fit.gam4)
# Compare to:
plot.nenana(fit.gam3)

# The tweedie parameter is generally difficult to estimate and may
# be highly uncertain, hence it is a good idea to check the estimate.
# We can use a handy function in the 'tweedie' package to construct
# a likelihood profile for the power parameter p:
library(tweedie)
# The function requires a linear model formula, hence I 
# approximated the fit using a high-order polynomial
# (hihger orders did not converge)
out <- tweedie.profile(Chum ~ poly(Julian, 8), data = Chum.R, 
                       p.vec = seq(1.025, 1.975, 0.025), do.plot=FALSE)
print(out$p.max)
plot(out, type="l")
# The maximum occurs at a value close to that estimated by 'mgcv'
# (1.29 vs. 1.21)
# The results from profiling also include an approximate 
# confidence interval:
out$ci

# The profile implies that the model with p=1.29 fits considerably
# better than higher values of p, including values close to 2, 
# which should be equivalent to a negative binomial model. This
# seems to contradict the results above suggesting the negative
# binomial model has a much higher likelihood (by over 7 units)
# and a lower AIC:
logLik(fit.gam3)
logLik(fit.gam4)
AIC(fit.gam3, fit.gam4)

####  Need to check into this some more!!!

# Fit model with p fixed at externally estimated value
# (using the 'Tweedie' family:
fit.gam5 <- gam(Chum ~ s(Julian, k=15) + offset(log(Volume)), data = Chum.R, 
                family = Tweedie(out$p.max))
summary(fit.gam5)
AIC(fit.gam3, fit.gam4, fit.gam5)
# The result is very similar, but has a slightly higher AIC despite
# having one less paramter (because we fixed p)

# Compare fits:
plot.nenana(fit.gam5) 
plot.nenana(fit.gam4)
plot.nenana(fit.gam3)


############
# Explore tweedie distribution using different amounts of "extra" variance 
# (as quantified using the 'phi' parameter - change phi and see how
# the variance and the number of zeros change)
N <- 2000
#(x <- rTweedie(rep(4,N), phi=3))
(x <- rTweedie(rep(4,N), phi=6))
hist(x)    
sum(x==0)   # Number of zeros



