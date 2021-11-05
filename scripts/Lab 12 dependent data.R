############################################################### 

# Lab 12 - Dependent data
#          Temporal and spatial autocorrelation

# Last modified: October 30, 2021
# Last run with R v. 4.0.3

###############################################################

library(nlme)
library(mgcv)
library(here)
### 1. Import global temperature data
Global.temp <- read.csv("data/Global temperatures.csv")

# Data from GISS Surface Temperature Analysis (GISTEMP)
# Link: http://data.giss.nasa.gov/gistemp/tabledata/GLB.Ts.txt
# Note that the data are different from those presented in 
# the lecture!

# Our goal is to fit a trend to the data and test the trend for 
# significance, while accounting for any autocorrelation in the 
# time series

## 1. Graphically examine temperature over time:
plot(Anomaly ~ Year, Global.temp, type="b",
     ylab = "Global mean temperature anomaly", 
     cex.axis=1.2, cex.lab=1.2, col=4)

## 2. Fit a linear regression, examine output, and add fitted line 
summary(fit.lm <- lm(Anomaly ~ Year, Global.temp))
abline(fit.lm,col=2)
# The linear model suggests that global average temperatures have been
#  increasing by 0.0947 +/- 0.0074 degrees per decade
#              (= slope +/- 2 * std err)


## 3. Examine anomalies over time (after removing linear trend)
oldpar <- par()  # save previous parameter settings
par(mar=c(2,6,2,1))    # Change margin width
plot(Global.temp$Year, resid(fit.lm), type="b", xlab="",
     ylab = "Global mean temperature anomaly\n(linear trend removed)", 
     cex.axis=1.3, cex.lab=1.2, col=4)
abline(h=0,lty=4,col=2)
par(oldpar)    # revert to previous parameter settings

# It is quite obvious that there is a remaining trend over time.
# That means that our simple regression statistics in step 2 are 
# not valid and underestimate variability in the parameters


## 4. Deal with autocorrelation
# We have (at least) two options: 
# A. We can try to approximate the trend using (for example)
#    a non-linear smoother to capture decadal or longer-scale 
#    variability in the data, which may eliminate any auto-
#    correlation in the data. In that case we lose the ability
#    to say something about the avarge annual increase

# OR:

# B. we can quantify autocorrelation and use it to either filter 
#    the data series prior to regression or incorporate the 
#    autocorrelation in the linear regression to account for its 
#    effect and get valid standard errors and confidence intervals

# For now, we start with the second option (B), that is we first identify 
# and quantify any autocorrelation in the residuals:

# Extract and save residuals:
r <- resid(fit.lm)
n <- length(r)   # number of years

# plot residuals for times t=2 through t=n (years 1881-2010) and plot them against 
# residuals one year prior (times t=1 through t=n-1  or years 1880-2009):
plot(r[1:(n-1)], r[2:n], asp=1, col=4) 
# Note that I fixed the aspect ratio to scale both axes the same
abline(0,1, lty=2)   # Add 1:1 line
# We can obtain an estimate of lag-1 autocorrelation by computing the
# simple Pearson's product moment correlation between the time series 
# and the time series lagged by one year:
cor(r[-1],r[-n])  # Note alternative way of selecting the shifted series!
# This is a valid estimate of first-order autocorrelation

# R has a convenient function that produces lagged plots for any lag:
lag.plot(r, do.lines=F, col=4, pch=16, diag.col=2)
# The default is a lag 1 plot. We can examine the series for correlation
# at other lags using the 'set.lags'' argument:
lag.plot(r, set.lags=2, do.lines=F, col=4, pch=16, diag.col=2)

# An estimate of the second-order autocorrelation can also be obtained by
# simply correlating the lagged time series, here achieved by removing the 
# first two observations and the last two observatons from r, respectively:
cor(r[-(1:2)],r[-((n-1):n)]) 

# lag plot and correlation at lag 3:
lag.plot(r, set.lags=3, do.lines=F, col=4, pch=16, diag.col=2)
cor(r[-(1:3)],r[-((n-2):n)])
# Note that we lose k observations when we shift the series to create
# a paired series at lag k (n-k overlapping years)!

# This methods allows an examination of one lag at a time. To quantify
# autocorrelation at multiple lags simultaneously we can use the 'acf'
# function, which  is designed to compute and plot autocorrelations
# at multiple lags up to a specified maximum:
acf(r, lag.max=25)
# The plot shows vertical lines indicating the magnitude of autocorrelations
# at each lag up to the specified maximum of 25 years

# The 'pacf' function first removes any autocorrelation up to lag k, then
# computes autocorrelation at lag k. This helps identify the order of the
# auto-regressive process!
pacf(r, lag.max=25)
# Note that there is some evidence for lag 4 autocorrelation. However, the 
# 95% confidence band should be treated as approximate and does NOT adjust for 
# multiple tests, therefore some estimated coefficients may be expected to
# exceed the lower or upper critical values by chance

# We can save and extract the estimated autocorrelation coefficients:
coefs <- acf(r)
coefs
# Note that the lag-1 coefficient (0.624) is very similar to the coefficient
# we estimated earlier. They are different because a slightly different
# definition is used by 'acf''

# An alternative test for lag-1 autocorrelation that we saw before is the
# Durbin-Watson test, which takes a model or a formula as input:
library(lmtest)
dwtest(lm(Anomaly ~ Year, Global.temp))
# The test strongly supports our conclusion that there is significant
# autocorrelation in the residuals.

#######################################
##  Filtering (this is OPTIONAL as filtering is not used much in ecology)
# The filtering assumes that any trends in the data are due to autocorrelation
# We can produce a filtered time series by removing the
# autocorrelation from a time series as follows:
X <- Global.temp$Year  # Save as 'X' and 'Y' for simplicity
Y <- Global.temp$Anomaly
n <- length(X)
r <- resid(lm(Y~X))  # save residuals from linear model (same as before)
ar1 <-  cor(r[-1],r[-n])   # Estimate of AR(1) parameter

# With our estimate of the AR(1) parameter (assuming that is the correct
# structure), we can compute a 'filtered' residual time series using the 
# anomalies Y that we computed above by removing the "predicted anomaly" 
# at time t [= phi * Y(t-1)], where phi is the AR(1) parameter, from the 
# observed values: 
V <- Y[-1] - ar1*Y[-n]
# In the absence of a trend, and if the time series has a first-order
# autoregressive structure, this should be a 'white noise' series (= 
# uncorrelated random noise with mean zero and some variance)

# We can then use the filtered series to test for a trend over time. First, 
# we need to remove one year from our "explanatory variable" (year) as we
# 'lost' one year of data when filtering:
U <- X[-n]  

# Let's plot the filtered variable:
plot(U,V, typ="b",col=4)

# We now regress the filtered time series on Year to test for and 
# estimate the linear trend:
summary(fit <- lm(V ~ U))
abline(fit)
acf(resid(fit))  
# Residuals from filtered series still show some evidence of autocorrelation
# at lag 4. We could construct a more complex (4th order) filter but we'll 
# revisit 4th order correlation below in a different context.

# Although we can now draw valid inferences, the estimated slope for the
# filtered time series is on a different, not easily interpreted scale!
#######################################

# Because of the drawbacks of filtering, and because it is not used much
# in ecology, I generally prefer an alternative approach that directly 
# incorporates the autocorrelation as an AR(1) term into a generalized 
# least-squares regression with correlated errors. 

# This approach (as we saw previously) is implemented in the 'nlme' library:
library(nlme)
# In this approach the autocorrelation coefficient (or "auto-regressive" 
# coefficient) is simultaneously estimated with the regression parameters 
# as just another parameter in the model (when fitting the regression 
# via least-squeres or maximum likelihood)

# We can fit the same linear trend model and specify the assumed AR(1) 
# correlation structure as follows, while providing our previous estimate 
# of the AR(1) parameter as a starting value (if needed):
X <- Global.temp$Year
Y <- Global.temp$Anomaly
fit.gls <- gls(Y ~ X, correl = corAR1(value = 0.59, form = ~ X), method="ML")
# Note that we should fit the model via ML to allow for mode comparisons 
summary(fit.gls)

# The estimate of phi as estimated internally is slightly higher than our 
# earlier estimates (estimated 'outside' the model).
# Note also that we can obtain an approximate confidence interval for 'phi':
intervals(fit.gls)
# A 95% confidence interval for the estimated AR(1) parameter does not
# include 0, suggesting signficant positive autocorrelation at lag 1

# We can also formally compare models with and without autocorrelation
# using a likelihood ratio test:
fit.gls0 <- update(fit.gls, correl = NULL)  # refit model without autocorrelation
lrtest(fit.gls0, fit.gls)  
# The likelihood ratio test also clearly rejects the hypothesis that the
# AR(1) coefficient is zero (as assumed by the simpler model 'fit.gls0')


################################################################################
################## Exercise 1:
# There was some evidence of higher-order (lag 4) autocorrelation in the residuals 

# Therefore, fit a 'gls' model that incorporates fourth-order autocorrelation. 
#  An auto-regressive term of order 4 implies that the residual temperature
#  at time t depends on the residuals at time period t-1, t-2, t-3, AND t-4
#  through 4 autoregressive paramters:
#       e(t) = phi1*e(t-1) + phi2*e(t-2) + phi3*e(t-3) + phi4*e(t-4) + nu
# where the 'nu' are independent, normally distributed errors

# We can incorporate this kind of correlation structure in the call to 'gls'
# using 'correl = corARMA(p=4)'. Look at the help file for 'corARMA', refit 
# the linear regression of Y on X inlcuding a fourth order correlation 
# structure, and save the model output (be sure to use 'method = "ML"):
#---------

#---------

# Examine the model output and formally compare models with 1st order and 4th
# order autocorrelation in the residuals (e.g. likelihood ratio test).
#---------

#---------

# Compare the slopes and their standard errors (and/or confidence intervals) 
# between the three models  with no, 1st order,and 4th order autocorrelation
#  (fit.lm, fit.gls, fit.gls4)
#---------

#---------

# What can you conclude? 
#---------

#---------

##################################################################################

# Plot fits:
plot(Anomaly ~ Year, Global.temp, type="b", xlab="",
     ylab = "Global mean temperature anomaly", 
     cex.axis=1.3, cex.lab=1.5, col=4)
# 95% pointwise confidence interval from simple linear regression (iid): 
se <- predict(fit.lm, se=T)$se.fit
lines(X, fitted(fit.lm), lwd=2)
lines(X, fitted(fit.lm) + 2*se, lty=2, col=2, lwd=2)
lines(X, fitted(fit.lm) - 2*se, lty=2, col=2, lwd=2)

# 95% pointwise confidence interval using adjustment factor to inflate standard errors:
# AR(1) parameter estimate:
r <- resid(lm(Y~X))  # residuals from linear model
ar1 <-  cor(r[-1],r[-n])   # Estimate of AR(1) parameter
# Adjustment factor:
adj <- sqrt((1+ar1) / (1-ar1))
adj  
# Plot adjusted confidence bands:
lines(X, fitted(fit.lm) + 2*se*adj, lty=2, col=3, lwd=2)
lines(X, fitted(fit.lm) - 2*se*adj, lty=2, col=3, lwd=2)

# Note that we can fit the same model (fit.gls) using the 
# 'arima' function for time series analysis, but it has a
# very different syntax (see help file)
with(Global.temp, arima(Anomaly, order=c(1,0,0), xreg = Year))
# Compare regression coefficients to those obtained with 'fit.gls'

# 'ggplot' version of the simple linear regression results 
# with 95% confidence intervals:
library(ggplot2)
ggplot(Global.temp, aes(Year, Anomaly)) + geom_path() + 
        geom_smooth(method="lm")

## As an alternative to correcting for or accounting for autocorrelation, we 
# can try to model any pattern in the residuals by modifiying the model 
# appropriately to make the observed pattern part of the estimated 
# "mean structure"" of the data. (option A above!)

# For example, we could incorporate additonal covariates (besides 'Year') in our
# model of global mean temperatures that may explain the observed trend around the
# linear increase (which is presumably due to changes in CO2 concentrations).

# The trend may be related to decadal scale variability in the atmosphere/ocean
# system. In the absence of a suitable covariate, we could simply capture the
# mean trend through a smooth function - we can use a gam model to do so and to
# help us chose the appropriate amount of smoothing:

fit.gam <- gam(Y~s(X))
summary(fit.gam)
plot(fit.gam, se=T, resid=T, pch=1, rug=F, shift=mean(Y),shade=T)
par(mfrow=c(2,2))
gam.check(fit.gam)  # No obvious issues
par(mfrow=c(1,1))

# We can compare the resulting model to our earlier models that fit a linear
# trend with residual autocorrelation. The models could be compared using 
# the AIC or likelihood ratio tests.

# VERY IMPORTANTLY, when comparing a gls (or 'lme' fit) with fits from other 
# models that were fit via maximum likelihood (such as 'gam' or 'lm' models)
# YOU HAVE TO FIT THE 'gls' model via 'ML' (as we did above): 
AIC(fit.gls, fit.gam)

# Exercise 2 (Questions):
# What does the result suggest about the trend in long-term global 
# temperatures? What are the advantages and disadvantages of the 
# linear trend line vs. the smooth trend in helping us characterize 
# long-term climate variability?
#----------------
# Your answer:

#----------------


# We are not quite done yet as we have not yet fully examined the GAM fit. 
# In particular, we have not yet evaluated the new residuals for 
# independence, which you get to do in the next exercise

##################################################################################
########## Exercise 3:
# Plot residuals from the GAM model over time and graphically and 
# statistically check for any autocorrelation in the residuals using 
# lag plots and more formal tests based on the (partial) autocorrelation
# function or other tests.
#---------------

#---------------


# What do you conclude?
#---------------

#---------------

################################################################################

# If you find significant autocorrelation, it is relatively easy to include
# autocorrelation in a GAM modeling framework using the 'gamm' function,
# which combines the ability to use smooth transformations of the predictor
# variables from the 'mgcv' package with the mixed-effects modeling 
# capabilities from the 'nlme' package (which allows for different 
# correlation structures in the residuals).

# Using the 'gamm' function, we can fit the same gam model, while also allowing 
# for auto-correlated errors:
fit.gamm <- gamm(Y ~ s(X), correlation = corAR1())
# The resulting fit has two components, one relating to the 'GAM' part, the
# other relationg to the 'lme' (random effects) part. We did not include a
# random effects component, hence we are primarily interested in the output
# relating to the 'gam' fit: 
summary(fit.gamm$gam)  # Examine the 'gam' part of the output (a list)

# However, the 'lme' part contains the estimate of the auto-regressive 
# parameter (and a bunch of other stuff that we won't worry about!):
summary(fit.gamm$lme)
# A formal model comparison of this model with 'fit.gam' is not straightforward
# but we can re-fit the same model using gamm without the autocorrelation:
fit.gamm0 <- gamm(Y ~ s(X))
# We need to use the 'lme' part to compare models as that has the overall
# likelihood of the fitted model:
lrtest(fit.gamm0$lme, fit.gamm$lme)
AIC(fit.gamm0$lme, fit.gamm$lme)

# Based on the evidence from these comparisons, the model with temporal 
# autocorrelation provides a more adequate fit to the data

# Using a different and maybe more intuitive comparison, we can assess whether  
# the auto-regressive coefficient (phi) is significantly different from zero by 
# examining its 95% confidence interval:
intervals(fit.gamm$lme)
# Based on the confidence interval for phi is there evidence that there is any 
# remaining autocorrelation that we should account for in the model?

# We can plot the 'gamm' fit using the 'gam' part of the model:
plot(fit.gamm$gam, se=T, resid=T, pch=1, rug=F, shift=mean(Y),shade=T)
# Note that the fit is virtually identical but the confidence interval
# is slightly wider to reflect the additional uncertainty that results 
# form having non-independent errors (and hence estimating an 
# additional parameter: phi)!

# Compare:
par(mfrow=c(1,2))
plot(fit.gam, se=T, resid=T, pch=1, rug=F, shift=mean(Y),shade=T)
plot(fit.gamm$gam, se=T, resid=T, pch=1, rug=F, shift=mean(Y),shade=T)
par(mfrow=c(1,1))



##################################################################################
##################################################################################

### A 'simple' spatial example:  (OPTIONAL)

# This portion of the lab is entirely optional but recommended for those
# of you who deal with spatial data. Be warned that this is a semester-long
# course on spatial statistics condensed into a single example but it should
# give you a flavor for how to deal with spatial correlations.


# We have a random sample of pollock abundances from the eastern Bering Sea
# for two recent years (2009, 2010) and wish to estimate mean densities by 
# year and compare densities between years while accounting for spatial
# correlation in the samples

library(maps)
library(mapdata)
library(mapproj)
library(lattice)
require(ggplot2)
require(nlme)

# Import data
EBS <- read.csv("data/EBS shelf.csv")
str(EBS)
EBS$Year <- factor(EBS$Year)

ggplot(EBS, aes(Longitude, Latitude)) + geom_point() + facet_wrap(~Year)

# We first convert latitude & longitude to x/y coordinates that reflect 
# actual distances using the equal-distance Albers projection:
# (I chose the parameters based on a 'rule of thumb' in the help file)
x <- mapproject(EBS$Longitude, EBS$Latitude, "albers", param=c(55.9, 60.8))
EBS$x <- x$x
EBS$y <- x$y

# Spatial distribution of sampling stations based on x,y coordinates,
# which ensure that distances reflect actual distances:
ggplot(EBS, aes(x, y)) + geom_point() + facet_wrap(~Year)

# Estimate means of log-transformed CPUE by year:
EBS$logP <- log(EBS$pollock+1)  # Create new variable for simplicity
f.lm <- lm(logP ~ Year, data=EBS)
summary(f.lm)
# Note that the intercept corresponds to the estimated mean log-density in 2009
# and the 'Year2010' paramter corresponds to the difference in mean log-density
# between 2009 and 2010

# Is there a significant difference in pollock density between 2009 and 2010?

# Examine residuals for autocorrelation. First, we can do a quick graphical 
# exploration by plotting the residuals in space:
library(sp)

# We should examine autocorrelation for one year at a time because the spatial 
# autocorrelation, if present, should be evident within years only (because the 
# fish may have a quite different spatial distributions between years):
EBS$r <- resid(f.lm)   # Extract residuals
coordinates(EBS) <- c("x","y")
# Note that this creates a more complex object that includes the data frame as
# a 'slot'. However, you can still use it like a data frame for most purposes.
j <- EBS$Year == "2010"  # Extract 2010 data only
bubble(EBS[j,], zcol="r")
# Note strong spatial patterns of all positive and all negative correlations!

# A formal test for spatial autocorrelation is a little less straightforward
# than for a time series because we now operate in two dimensions:

# The idea is to plot the squared difference in the values of the variable of 
# interest (here: log(CPUE)) against the distance separating the samples!

# This is called a 'Variogram' (for historical and mathematical reasons the 
# 'Semi Variogram', or half the Variogram, is used). 

# For each pair of observations x1,x2, the corresponding semi-variogram 
# is (x1-x2)^2/2, which is plotted against the distance between x1 and x2

# This is sort of the inverse of plotting correlations against distance.
# That is, whereas correlations are highest at short distances, the variogram
# (which is smaller the more similar observations are to each other) will be 
# smaller at short distances if there is autocorrelation (i.e. observations 
# that are close to each other in space are, on average, more similar to each 
# other than observations that are very far apart).

# We should examine autocorrelation for one year at a time because the spatial 
# autocorrelation, if present, should be evident within years only (because the 
# fish may have a quite different spatial distribution between years):
j <- EBS$Year == "2010"   # Select year
r <- resid(f.lm)[j]   # Extract residuals
# Compute pairwise distances among locations based on distances in 'x'
d <- dist(coordinates(EBS)[j,])
d  # Large (88x88) matrix of pairwise distances (values below diagonal only!)
d <- as.vector(d)   # Need to convert to a vector for 'Variogram' function
# The vector d, contains only one set of pairwise distances (values from below
# the diagonal) and has length: n * (n-1) / 2, where n is the number of observations
# (see help files for 'dist' and 'Variogram')

# The 'Variogram' function from the 'nlme' package computes the semi-variogram,
#  given a variable (residuals in this case) and the corresponding pairwise 
#  distances 
SemiVar <- Variogram(r, d)
head(SemiVar, 10)  
# Each row of the result corresponds to a pair of observations and the column
# 'variog' contains the semivariogram values and the column 'dist' the 
# corresponding distances between samples 

# There is a default plotting option that simply plots the values against
# distances with a scatterplot smoother, which indicates an increase in
# the semivariogram with distance, as would be expected:
plot(SemiVar, xlim=c(0,0.1))
# Note that the variance increases with distances, that is: the further 
# two points are apart the less similar to each other they tend to be.
# (The inverse of this is that the closer two observations are together
# the more similar - i.e. positively correlated - they are)

# Because of the large number of pairwise distances, it is typically more
# useful to plot binned values of the semi-variogram (which is the default
# for many spatial packages such as 'geoR')

# We create bins for the distances in increments of 0.005 and plot the 
# distribution of the semi-variogram values in each bin:
bins <- cut(SemiVar$dist, seq(0,0.1, by=0.005))  
plot(bins, SemiVar$variog)
# This clearly shows that differences in observations that are the smallest 
# distance apart, tend to have low variance, i.e. are similar to each other!

# The modeling functions such as 'gls' and 'lme' then fit a model that 
# quantifies how the semi-variogram increases with distance, for example 
# using an exponential, Gaussian, or spherical model.

# We can estimate these (and other) models directly using the 'gls' function,
# which simultaneously finds the maximum likelihood estimates of the parameters 
# of a model relating the response to the explanatory variable (in our simple 
# case a model with two separate means for 2009 and 2010) and the parameters 
# of the variogram model, for example an exponential model

# We include a "nugget" effect in the model in case there is measurement 
# error at very short distances (that is, repeated measurements taken at 
# the same exact point do not have identical values but have 
# some non-zero variance due to small-scale variability.

# We can include an appropriate correlation structure using one of the 
# spatial variogram functions, such as 'corExp' for an exponential model
# The 'form' specifies that the distances of points are specified by
# the coordinates x and y within a given year (~ x+y | Year)
f.exp1 <-gls(logP ~ Year, data=EBS, correl=corExp(form= ~ x+y | Year, nugget=T))
summary(f.exp1)
# The output includes the year effect, as well as estimates of the 'range' 
# and nugget parameters. The former corresponds to the distance over which 
# autocorrelation extends (i.e. roughly where the exponential function starts
# to flattten out) and the latter to the variance at zero distance.
# Note that the estimated nugget effect is very small and we may not need it
# Hence let's fit a model without the nugget effect and compare:
f.exp2 <-gls(logP ~ Year, data=EBS, correl=corExp(form= ~ x+y | Year, nugget=F))
summary(f.exp2)
AIC(f.exp1, f.exp2)  # Suggests that model without nugget effect is adequate!

# We can easily visualize the estimated variogram model:
plot(Variogram(f.exp2), ylim=c(0,1.1))
# Note that the y-axis ranges from 0 to 1 with 0 being no difference in observed 
# values and 1 representing the maximum difference (corresponding to the "sill") 

# Unfortunately, these functions are not well documented and it is not quite 
# clear how the semi-variogram is standardized to range from 0 to 1 or how 
# exactly the 'range' parameter is defined without some digging

# Other spatial packages such as 'geoR' are much more flexible and allow 
# estimation of the range, nugget, and "sill" (i.e. the maximum value of 
# the semi-variogram)

################
# Question: Based on these models, which account for spatial autocorrelation 
# in CPUE, what can you conclude about differences in pollock CPUE between 
# 2009 and 2010? Compare to the result from the linear model! (f.lm)

# Answer: To test the null hypothesis that there is no significant difference in
# CPUE between 2009 and 2010, we can simply use the Wald's t-test in the 'summary'
# output for the 'Year2010' coefficient. The estimated difference (on the log-scale)
# is 0.0284: 
round(coef(summary(f.exp2)),4)
# The standard error is large, hence the t-test cannot reject the null hypothesis
# that the coefficient (difference between years) is equal to zero. Note that our
# conclusions change after taking spatial autocorrelation into account. 
# Compare to (see also Exercise 4 below):
round(coef(summary(f.lm)),4)
################

# Although the exponential model is probably adequate, we can fit other 
# correlation structures and compare, for example, a Gaussian and a 
# spherical model as follows:
f.Gaus <-gls(logP ~ Year, data=EBS, correl=corGaus(form= ~ x+y | Year, nugget=F))
summary(f.Gaus)
plot(Variogram(f.Gaus), ylim=c(0,1.1))

f.Spher <-gls(logP ~ Year, data=EBS, correl=corSpher(form= ~ x+y | Year, nugget=F))
summary(f.Spher)
plot(Variogram(f.Spher), ylim=c(0,1.1))
# Note the relatively poor fit for both the Gaussian and spherical model 
# based on visual examination alone

# This tends to happen quite frequently and it is important to check 
# your fits!! It is not uncommon for these models to fail to converge. 
# In that case, don't give up! You can provide starting values for the 
# parameters of the autocorrelation function

# For example, if the spherical model fails to converge, we could use the 
# estimate of the range parameter from the exponential model above as a 
# starting value. This can be specified by the 'value' argument of the 
# corSpher() function, which takes a single value for the range parameter 
# if nugget = F, or a vector of two values if nugget = T:
# (see help file for 'corSpher')
f.Spher <-gls(logP ~ Year, data=EBS, correl=corSpher(value = 0.0175, 
                                     form= ~ x+y | Year, nugget=F))
summary(f.Spher)
plot(Variogram(f.Spher), ylim=c(0,1.1))

# This spherical model does not fit particularly well based on a viusal 
# assessment alone, hence we may want to inlcude a nugget for this model 
# (with starting value for range from previous fit):
f.Spher2 <-gls(logP ~ Year, data=EBS, 
              correl=corSpher(value=0.033, form= ~ x+y | Year, nugget=T))
summary(f.Spher2)
plot(Variogram(f.Spher2), ylim=c(0,1.1))   # Much better!

# We can compare these (non-nested) models using AIC, which clearly favors
# either the exponential model without a nugget or the spherical model with
# nugget effect:
AIC(f.exp2, f.Gaus, f.Spher2)

# For parsimony we may want to chose the exponential model, but the choice has 
# little effect on parameter estimates and their uncertainty. As long as you 
# include a "reasonable" spatial correlation structure, you
# should be OK to use the model for drawing inferences.


################################################################################
############## Exercise 4:
# Compare the estimated differences in mean CPUE between 2009 and 2010 
# from these models and the associated standard errors with those from 
# the simple linear model!
#---------------

#---------------

###############################################################################


# We can easily include additional covariates in the model that may help 
# explain the observed spatial patterns in density, such as a dome-shaped 
# temperature effect, reflecting an optimum temperature for walleye pollock:

# For model comparisons, you need to fit all models using 'ML':
f.exp3 <-gls(logP ~ Year + poly(Temperature,2), data=EBS, 
             correl=corExp(form= ~ x+y | Year, nugget=F), method="ML")
summary(f.exp3)
# Check the variogram:
plot(Variogram(f.exp3), ylim=c(0,1.1))  # looks reasonable
# re-fit 'Year' only model using ML for model comparison:
(f.exp4 <- update(f.exp2, method="ML"))
summary(f.exp4)

# Compare models via AIC:
AIC(f.exp3, f.exp4)
# Or use a likelihood ratio test:
anova(f.exp3, f.exp4)
# This suggests that the model with a temperature effect fits 
# substantially better than one without! 

# Hence, we can conclude that the local CPUE of pollock appears to be 
# substantially affected by temperature (assuming the same relationship 
# between temperature and CPUE in both 2009 and 2010), but accounting 
# for differences in temperature does not affect our conclusion that there 
# is no significant difference in CPUE between 2009 and 2010 
# (see 't-table' from summmary output). Note also that the estimate of the 
# range (which quantifies the extent of spatial autocorrelation) is smaller 
# after including the temperature effect (compare estimates of 'range' 
# from the summary outputs). This is to be expected because including 
# temperature (which is itself spatially autocorrelated) in the model  
# removes some of the autocorrelation in the residuals.

# Let's examine the temperature effect by setting up a data frame that 
# includes a range of temperature values for 2010 (we only need one year 
# because the temperature effect is assumed to be the same for both years):
new.df <- data.frame(Year = factor("2010"), Temperature = seq(-1.7, 6.3, by=0.1))
# Note that you need to make the 'year' variable  a factor, to match the 
# original data frame (EBS), othrewise 'predict' will not work. 

# Compute and save predicted values and plot result:
new.df$pred <- predict(f.exp3, newdata=new.df)
plot(pred ~ Temperature, data=new.df, type="l", ylab="log(pollock CPUE)", col=2)
# The temperature effect is quite substantial (note range of y-values) and we
# can compute a 'pseudo' R2 value (analogues to other linear models) as 1 minus
# the portion of the variability in log(CPUE) that remains unexplained. The 
# latter is estimated by dividing the sum of squared residuals from the model 
# by the sum of squared residuals from the 'null' model. The null model is a 
# model that includes an intercept (overall mean) only, so we can simply remove 
# the mean and use the sum of the squared anomalies for the denominator (same 
# as fitting an intercept only model (logP ~ 1) and taking the sum of squared 
# residuals):
1 - sum(resid(f.exp3)^2) / sum((EBS$logP - mean(EBS$logP))^2)
# Thus, temperature accounts for nearly a quarter of the variability in CPUE,
# which is quite a bit for highly variable CPUE data!

# Note that this does not account for patterns in the data due to spatial 
# autocorelation, which "explains" additional variability in CPUE. Unfortunately,
# the actual spatial pattern that is 'induced' by autocorrelation is not easy to 
# extract from the 'gls' output ('geoR' has functionality to do so) because the
# focus of 'gls' is on the main effects in the model (i.e. Year and Temperature)
# rather than the spatial pattern due to autocorrelation.

# We could compare the magnitude of the temperature effect to the estimated 
# differences between 2009 and 2010 by superimposing the estimates for 2009:
new.df <- data.frame(Year = factor("2009"), Temperature = seq(-1.7, 6.3, by=0.1))
new.df$pred <- predict(f.exp3, newdata=new.df)
lines(pred ~ Temperature, data=new.df, type="l", ylab="log(pollock CPUE)", col=4)
# Note that the estimated CPUE is virtually identical between the two years
# hence the two lines are almost indistinguishable. Of course the Year 
# effect was not significant in the model, so this is to be expected.

# This concludes the spatial modeling lab!! Thanks for hanging in there....
