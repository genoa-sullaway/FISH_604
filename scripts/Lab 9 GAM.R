################################################

# Lab 9   Generalized Additive Modeling
#          Skate CPUE data
# 
# Author: Franz Mueter

# Last run using R 4.1.0 
# Last modified: October 19, 2021

################################################
library(here)
# I. Skate data, Gulf of Alaska

### 1. Import data:
CPUE <- read.csv("data/cpue.csv")
head(CPUE)
str(CPUE)

# Examine spatial distribution of data:
# (you will need to be online for 'leaflet' as it 
#  used OpenStreetMap)
library(leaflet)
library(plyr)
CPUE$year.col <- mapvalues(CPUE$YEAR, unique(CPUE$YEAR), 
                           c("slategray", "tomato", "violet", "yellow", "steelblue"))
leaflet(CPUE) %>% addTiles() %>% 
  addCircles(~LONGITUDE, ~LATITUDE, color=~year.col) 
# Note this creates 

# Scatterplot by depth with LOESS smooth:
scatter.smooth(CPUE$DEPTH, CPUE$SKATES)

## 1a.
hist(CPUE$SKATES)
hist(CPUE$SKATES^0.25)

#---------------------------------------------
## Exercise 1b.  What is troubling about the distribution 
#       of the fourth-root transformed CPUE data if 
#       we want to statistically model CPUE?

# Answer:

#---------------------------------------------

### 2. Analysis of CPUE-where-present

## 2a. Create subscript to only use non-zero catches
sub <- CPUE$SKATES > 0

## 2b. Attach data frame and add scatterplot smoother 
#      to fourth-root transformed values:

with(CPUE, scatter.smooth(DEPTH[sub], SKATES[sub]^0.25)) # spans are different window widths for loess smooth

attach(CPUE)   # Using 'attach' again because the functions 
               # below don't have a 'data' argument!!!

## 2c. Compare different spans:
fit <- loess.smooth(DEPTH[sub], SKATES[sub]^0.25, span=0.2) #red
lines(fit, col=2)
fit <- loess.smooth(DEPTH[sub], SKATES[sub]^0.25, span=0.5) #green
lines(fit, col=3, lwd=2)
fit <- loess.smooth(DEPTH[sub], SKATES[sub]^0.25, span=1) #blue
lines(fit, col=4, lwd=2)

fit <- loess.smooth(DEPTH[sub], SKATES[sub]^0.25, span=2) #blue
lines(fit, col=4, lwd=2)

#--------------------------------------------------------------
## Exercise 2c - Questions:
#   i.  Does the default setting for span (=2/3) produce a 
#       "reasonable" fit? 
#   ii.  Choose a "best" fit (span) subjectively (try other values)
#        and justify your choice

# Answers:

#--------------------------------------------------------------

### 2d. Determine the "best" smoothing spline fit by 
#       leave-one out cross-validation (CV) as 
#       implemented in 'smooth.spline':
plot(DEPTH[sub], SKATES[sub]^0.25)
fit2 <- smooth.spline(DEPTH[sub], SKATES[sub]^0.25)
lines(fit, col=2)
summary(fit)

### 2e. Compare your subjective "best fit" to the one 
#       that minimizes the CV score!
#--------------------------------------------------------------
# Your response:


#--------------------------------------------------------------

### 3a. Fit LOESS smooth in four dimensions:
skate.lo <- loess(SKATES^0.25 ~ DEPTH * AD * JULIAN * NET.WIDTH,
 	data=CPUE, subset=sub)
summary(skate.lo)

### 3b. Diagnostic plots:
r <- resid(skate.lo)
# Residuals against fitted values:
scatter.smooth(fitted(skate.lo), r, lpars=list(col=2)); abline(h=0, lty=2)
# Residuals by depth:
scatter.smooth(DEPTH[sub], r, lpars=list(col=2)); abline(h=0, lty=2)
# (note that we select a subset of the 'DEPTH' vector. This has to
#  be the same subset used when fitting the model above!)
hist(r)
qqnorm(r)
abline(0,1)
#----------------------------------------------------------
# Exercise: Examine additional diagnostic plots for
# AD
# JULIAN
# NET.WIDTH

# and other suitable diagnostics

# Solution (Code): 

# Summarize your conclusions from residual diagnostics
# Answer:

#----------------------------------------------------------


### 3c. Re-fit model with smaller span (more degrees of freedom)
skate.lo <- loess(I(SKATES^0.25) ~ DEPTH * AD * JULIAN * NET.WIDTH,
 	data=CPUE, subset=sub, span=0.3)

scatter.smooth(DEPTH[sub], r, lpars=list(col=2)); abline(h=0, lty=2)

scatter.smooth(AD[sub], r, lpars=list(col=2)); abline(h=0, lty=2)

scatter.smooth(NET.WIDTH[sub], r, lpars=list(col=2)); abline(h=0, lty=2)

scatter.smooth(JULIAN[sub], r, lpars=list(col=2)); abline(h=0, lty=2)

CPUE[sub,][fitted(skate.lo) < 2 | fitted(skate.lo) > 8,] 
plot(fitted(skate.lo), SKATES[sub]^0.25)
abline(0,1)

#---------------------------------------------------------------
# Exercise: Examine residual plots by re-running your code above 
# with the new fit! Comment on the residual diagnostics.
# (Since we named the model the same, just re-run lines 101-108)

# Answer:

#---------------------------------------------------------------


# Identify outliers and compare to full range of data values:
# Plot of residuals against fitted values:
r <- resid(skate.lo)
scatter.smooth(fitted(skate.lo), r, lpars=list(col=2)); abline(h=0, lty=2)
# The plot does not show any large outliers in the residuals, but two
# fairly extreme fitted values (< 2 and > 8):
CPUE[sub,][fitted(skate.lo) < 2 | fitted(skate.lo) > 8,] 

# We can compare these values to the range of observed values
# which we can summarize as follows
apply(CPUE[sub,c("DEPTH", "AD", "JULIAN","NET.WIDTH", "SKATES")], 2, range)
# This shows minimum and maximum depth, minimum and maximum AD, etc
# Note that the unusual values are both near the southeastern corner 
# of the survey area (see lat, long, and AD)

### 3d. Observed against fitted values and R^2 value:
plot(fitted(skate.lo), SKATES[sub]^0.25)
abline(0,1, col=2)

#--------------------------------------------------------
# Exercise: Compute "Pseudo-R2" and interpret the output

# Solution:



#--------------------------------------------------------


### 3e. Plot fitted surface by depth & AD

# Set up grid for predictions:
(x <- seq(min(DEPTH), max(DEPTH), length=40))
(y <- seq(min(AD), max(AD), length=40))

# We predict values over this grid, for a sampling date of
# June 30 and for mean net width. First, we need to set up
# a data frame with the predictor variables:
pred <- expand.grid(DEPTH=x, AD=y)
pred$JULIAN <- 180  # Arbitrarily set to day 180 (June 30)
pred$NET.WIDTH <- mean(NET.WIDTH[sub])  # Set to mean net width
head(pred) # Depth / AD combinations for predictions
plot(DEPTH ~ AD, data=pred)  # Grid of values for predictions (40 x 40)

# Compute predicted values on our 40 x 40 grid:
p <- predict(skate.lo, newdata = pred, se=T)
names(p)
p$fit   # Note that the result is a 40x40 matrix
dim(p$fit)

# Plot the resulting fit:
image(x, y, p$fit, xlab= "Depth", ylab = "AD", col=topo.colors(60))
contour(x, y, p$fit, add=T)
points(DEPTH[sub], AD[sub], pch=16, cex=0.5)

# Plot standard errors:
image(x, y, p$se.fit, xlab= "Depth", ylab = "AD")
contour(x, y, p$se.fit, add=T)

# Note large uncertainty in the corners, where 
# the model is extrapolating

# We may want to eliminate the regions where the standard errors
# exceed some threshold, say larger than 2:
p$fit[p$se.fit > 2] <- NA
# This replaces any fitted values in 'p$fit' whose standard
# error exceeds 2 with 'NA', which suppresses plotting:
image(x, y, p$fit, xlab= "Depth", ylab = "AD", col=topo.colors(60))
contour(x, y, p$fit, add=T)
points(DEPTH[sub], AD[sub], pch=16, cex=0.5)

# A different view of the fitted surface:
persp(x, y, p$fit)

detach(CPUE)	# Detach CPUE data frame


####################################################################
# II. Pollock example:

### 1. Import data
pollock <- read.csv("data/pollock.csv")
pollock
str(pollock)

### 2. Fit GAM model:
library(mgcv)
pollock.gam <- gam(log.surv ~ s(Bloom) + s(Wind), data = pollock)
summary(pollock.gam)

#------------------------------------------------------------------
# Exercise: What do the (approximate) F-tests suggest about the 
# significance of the terms? What do the estimated degrees of freedom
# suggest about linearity or non-linearity of the terms? Could the 
# terms be approximated by low-order polynomials?

# Answer:
# bloom and wind speed are related to survivasl
# EDF of 1 means you can replace that terms with a linear forms
# if we want could replace the second term (edf =2) with a quadratic term, then you would just conduct the model in a linear framework

#------------------------------------------------------------------

# Plot model fits:
par(mfrow=c(1,2))
plot(pollock.gam, shade=T, seWithMean=T, resid=T, pch=1)
# Note the argument 'seWithMean=T', which accounts for uncertainty
# about the overall mean. Generally preferred as a better confidence
# interval (i.e. it more accurately reflects true uncertainty, as
# opposed to the default seWithMean=F, which ignores the uncertainty
# in the intercept and other linear parameters)

#------------------------------------------------------------------
# Exercise: Based on the spread of the smooth terms along the y-axis, 
# does one of the variables have a stronger effect on log-survival?

# Answer: 

#------------------------------------------------------------------

### 3. Diagnostic plot:
par(mfrow=c(2,2))
gam.check(pollock.gam)
par(mfrow=c(1,1))
r <- resid(pollock.gam)
scatter.smooth(fitted(pollock.gam), r)
scatter.smooth(pollock$Bloom, r)
scatter.smooth(pollock$Wind, r)

# Check for serial correlation (temporal autocorrelation):
plot(pollock$Year, r, type="b")
# Durbin-Watson test to check for serial correlation:
library(lmtest)
dwtest(pollock.gam)
# There is significant autocorrelation that we should take
# into account (but we won't, right now!)


### 4. Fit semi-parametric model that fixes degrees of freedom 
#      for bloom date (linear):
pollock.sem <- gam(log.surv ~ Bloom + s(Wind), data = pollock)
summary(pollock.sem)


### 5. Bootstrap confidence intervals for predicted values

## 5a. Re-fit the model to fix the degrees of freedom for the wind 
# term at df=2 (k=3):
fit <- gam(log.surv ~ Bloom + s(Wind, k=3, fx=T), data = pollock)
summary(fit)

#you can use AIC because likelihood approach
#GCV is used internally by GAM to determine the amount of smoothing (the edf, and/or k values)

## 5b. Parametric bootstrap:
#if you residuals area  little weird, and/or you have small sample sizes consider a bootstrap for the confidence intervals. 

# Conditional on a model, trying to simulation other random realizations that oculd arise 
# from the same model. one way to do that is take all residuals, and resample with
# replacement, then use residuals to add points to the line, calculate CI. only log(survival) y values
# change, x values do not change. 
# this creates a dataset that mimics our data, creating another realization of
# what could have happened in nature. 

# Bootstrap samples for the parametric bootstrap are 
# obtained by "scrambling" (re-sampling) the residuals 
# randomly and adding them to the fitted values from 
# the model to get a new 'bootstrapped' response variable
# that forms the bootstrap dataset along with the 
# original predictor variables. We then fit the same 
# model to the bootstrap dataset and save the results.
# We do this many times to assess uncertainty in the
# resulting outputs (parameters or predicted values).

# Here we compute confidence bands for 40 predicted
# values at the mean bloom date (150) and at 40 levels
# of wind mixing that span the range of observed wind
# mixing (285-590):
new <- data.frame(Bloom = rep(150, 40), 
                  Wind =  seq(285,590,length=40))
# Here create a new data frame to add predicted values to plot fitted lines CIs etc.
# The 'new' data frame is used by the 'predict' function
# below to compute predicted values at the desired levels
# of wind mixing and bloom date.

# We need both the fitted values and the residuals from
# our original model for the parametric bootstrap:
.fitted <- fitted(fit)
.resid <- resid(fit)
# Make a copy of the data for bootstrapping, only including
# the variables we need:
boot.dat <- pollock[,c("log.surv", "Bloom", "Wind")]

R <- 1999 # Number of bootstrap replicates
# Create matric for output:
out <- matrix(NA, nrow=R, ncol=40)

for(i in 1:R) {
  # Create bootstrapped value for log-survival:
  boot.dat$log.surv <- .fitted + sample(.resid, replace=T)
  # Fit the model to this new bootstrapped data set:
  fit.boot <- gam(log.surv ~ Bloom + s(Wind, k=3, fx=T), data = boot.dat)    
  # Save predicted values:
  out[i,] <- predict(fit.boot, new)
}

# Let's look at some of the bootstrap fits:
x <- new$Wind  # Wind values at which to evaluate fit and CIs
y.hat <- predict(fit, new) # Predicted values
plot(x, y.hat, ylim = c(0.5,3.5), type="l", lwd=3, col=4)
for(i in 1:15) lines(x, out[i,], col=2)

## 5c.  Plot (percentile-based) confidence intervals
# Compute lower (2.5th) and upper (97.5th) percentiles for a 
# percentile-based confidence interval at each predicted value:
ci <- apply(out, 2, quantile, probs=c(0.025, 0.975))
ci
lower <- ci[1,]
upper <- ci[2,]

plot(x, y.hat, ylim = range(c(lower, upper)), type="l", lwd=2)
lines(x, lower, lty=2)
lines(x, upper, lty=2)

# Predicted values and confidence band from GAM fit:
p <- predict(fit, newdata = new, se=T)
lines(x, p$fit + 2*p$se.fit, lty=2, col=2)
lines(x, p$fit - 2*p$se.fit, lty=2, col=2)


# ----------------------------------------------------------------------
# Exercise: Do the confidence intervals returned by gam() 
# agree approximately with those obtained by bootstrapping?

# Answer

# ----------------------------------------------------------------------


## 5d. Accounting for uncertainty in estimating edf

# In the bootstrap approach above, we ignored uncertainty in the
# appropriate amount of smoothing by fixing k = 3. To properly account
# for this uncertainly, we can redo the bootstrap and estimate the 
# best value for the equivalent degrees of freedom (edf) at each step

# We re-fit the model without fixing df. To prevent over-fitting, we 
# assume that the degrees of freedom are no more than 4 (k=5), which
# allows for a lot of flexibility in what should be a relatively simple
# relationship:
fit <- gam(log.surv ~ Bloom + s(Wind, k=5), data = pollock)
summary(fit)  # Note edf is close to but not exactly 2!

# The only difference between this version of the bootstrap and the
# above approach is that the appropriate d.f. of the smooth 
# term is now estimated each time when fitting the model to a 
# new bootstrap data set (rather than fixed at df=2 (k=3).

# Set up a new matrix for output from modified bootstrap:
R <- 1999
out2 <- matrix(NA, nrow=R, ncol=40)

for(i in 1:R) {
  # Create bootstrapped value for log-survival:
  boot.dat$log.surv <- .fitted + sample(.resid, replace=T)
  # Fit the model to this new bootstrapped data set:
  fit.boot <- gam(log.surv ~ Bloom + s(Wind, k=5), data = boot.dat)    
  # Save predicted values:
  out2[i,] <- predict(fit.boot, new)
}

# -------------------------------------------------------------------

# Exercise:
# Overlay the resulting bootstrap confidence bands on the previous plot 
# (follow example above for first bootstrap results)! Does the bootstrap 
# confidence interval differ much from the previous intervals?

# Solution:


# Answer:


# ----------------------------------------------------------------


##########################################################
# III. Fit GAM to skate data:

# Recall the skate data from section (I.) above.
# If you did not save your results from that section,
# re-run the code in section I first 

### 1. Fit Generalized Additive Model
skate.gam <- gam(SKATES^0.25 ~ s(DEPTH) + s(AD) + s(JULIAN) 
 	+ s(NET.WIDTH), data=CPUE, subset = sub)
summary(skate.gam)

## 1a. Plot results with partial residuals
par(mfrow=c(2,2))
plot(skate.gam, resid=T, seWithMean=T)  # Plot fits with partial residuals
plot(skate.gam, seWithMean=T)  # Plot fits without partial residuals

## 1b.Diagnostic plots:
gam.check(skate.gam)

#-------------------------------------------------------------------
## Exercise: Compare the estimated R^2 value for this model (see 
#      'summary' output above to the one we estimated from 
#      the loess fit (see I.3.d).  
#      Which model fits "better" based on the R2 values? Why? 
#      (compare equivalent number of parameters or edf)

# Solution & answer:
summary(skate.lo)
summary(skate.gam)

skate.lo #51.5 parameters 
skate.gam #24 edf

# higher parameters = better model -- how is this related here?! 
#-------------------------------------------------------------------

### 2. Model "interaction" between depth and AD:
skate.gam2 <- gam(SKATES^0.25 ~ s(DEPTH, AD, bs='tp', k=40) + s(JULIAN) 
 	+ s(NET.WIDTH), data=CPUE, subset = sub)
summary(skate.gam2)

anova(skate.gam, skate.gam2, test="Chisq")
AIC(skate.gam, skate.gam2)

skate.gam2 <- gam(SKATES^0.25 ~ s(DEPTH, AD, bs='tp', k=50) + s(JULIAN) 
                  + s(NET.WIDTH), data=CPUE, subset = sub)
summary(skate.gam2)

anova(skate.gam, skate.gam2, test="Chisq")
AIC(skate.gam, skate.gam2)

# We can visualize the 'interaction' as follows:
vis.gam(skate.gam2, c("DEPTH", "AD"), plot.type="contour", color="topo")
points(AD ~ DEPTH, data = CPUE[sub,], cex= 0.5, pch=16)

vis.gam(skate.gam2, c("DEPTH", "AD"), plot.type="contour", color="topo",
        too.far=0.05)
points(AD ~ DEPTH, data = CPUE[sub,], cex= 0.5, pch=16)

#-------------------------------------------------------------------
# Exercise: Which model would you choose as the "best" model 
# based on these comparisons?

# Answer: 


#-------------------------------------------------------------------

### 3. Re-fit without 'NET.WIDTH'
skate.gam3 <- update(skate.gam2, ~.-s(NET.WIDTH))

AIC(skate.gam,skate.gam2,skate.gam3)
# -------------------------------------------------------------------
# Exercise: Compare this model to skate.gam2 using anova() and AIC(). 
# Which is the overall best model among the three models we fit so far?

# Solution & Answer: 2 and 3 are the same but because 3 is simpler would opt for that, and it has fewer edf.  

# -------------------------------------------------------------------


### 4. Add YEAR as categorical variable to best model:
skate.gam4 <- update(skate.gam3, ~ . + factor(YEAR)-1)
summary(skate.gam4)
AIC(skate.gam2, skate.gam3, skate.gam4)

plot(skate.gam3, select=2, seWithMean=T)
plot(skate.gam4, select=2, seWithMean=T)

# Adding year substantially improves the model, suggesting 
# large differences in skate CPUE among years. Therefore
# pooling years is clearly not appropriate and we should
# probably have included year from the get-go!
# -------------------------------------------------------------------
summary(skate.gam4)
plot(skate.gam4, all.terms=T, select=3, seWithMean=T, rug=F)

tapply(CPUE$JULIAN[sub], CPUE$YEAR[sub], mean)

AIC(skate.gam3, skate.gam4)
unique(CPUE$YEAR)
plot(skate.gam4)

########################### OPTIONAL:

# This is a follow up on the bootstrap exercise in section II and 
# illustrates the use of the 'boot' package, an efficient way for 
# conducting bootstrap analysis. The basic steps are the same as in
# section II.5 above, but uses a slightly more complex syntax that
# takes some getting used to. The main reference for the 'boot' 
# package is:
#   Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and 
#   Their Application. Cambridge University Press.

# To run this example, make sure to run section II, 1 - 5a above, then 
# proceed with the following steps:

## 5b. Bootstrap function 'boot'
# (the subscript i in the function is used in the function 'boot' below):
# the function resamples the residuals, creates a new bootstrap dataset,
# fits the model to the bootstrapped dataset, and returns predicted values

# First, we append fitted values from the model to the original 
# date frame, as we will need them for bootstrapping:
dat <- cbind(pollock, fit = fitted(fit))

# Writing a bootstrap function to work with 'boot':
bf <- function(resids, i, x) {
  # create a bootstrap data set by adding re-sampled residuals 
  # to the fitted values:
  dat$log.surv <- dat$fit + resids[i] 
  # Create a vector of x-values at which to compute predicted values
  new <- data.frame(Bloom = rep(150, 40), Wind =  x)
  # Fit the model to this new bootstrapped data set:
  fit.boot <- gam(log.surv ~ Bloom + s(Wind, k=3, fx=T), data = dat)    
  # Return predicted values:
  predict(fit.boot, new)   
}

## 5c. Run the bootstrap
library(boot)
# The 'boot' function uses the above function ('bf') to do the bootstrapping,
# i.e. it runs bf with 1999 different sets of resampled residuals from 'resid(fit)':) 
x <- seq(285,590,length=40)          
pollock.boot <- boot(resid(fit), bf, x=x, R = 1999)

# Confidence intervals for individual predictions 
# (predicted survival at the first, second, 10th and 40th wind values)
# The 'boot.ci' function can return four different types of confidence intervals,
# we extract only the percentile-based interval:
boot.ci(pollock.boot, type = "perc", index = 1)   # Confidence intervals for fitted value at x[1] = 285
boot.ci(pollock.boot, type = "perc", index = 2)   # CI at x[2] = 292.8
boot.ci(pollock.boot, type = "perc", index = 10)
boot.ci(pollock.boot, type = "perc", index = 40)

# Plot the best model fit to the data (t0) as well as a few individual bootstrap
# fits over the full range of x-values (wind):
plot(x, pollock.boot$t0, ylim = c(0.5,3.5), type="l", lwd=3, col=4)
for(i in 1:15) lines(x, pollock.boot$t[i,], col=2)

## 5d.  Plot (percentile-based) confidence intervals
# Compute lower (2.5th) and upper (97.5th) percentiles for a 
# percentile-based confidence interval at each predicted value:
ci <- apply(pollock.boot$t, 2, quantile, probs=c(0.025, 0.975))
ci
lower <- ci[1,]
upper <- ci[2,]

plot(x, pollock.boot$t0, ylim = range(c(lower, upper)), type="l", lwd=2)
lines(x, lower, lty=2)
lines(x, upper, lty=2)

# Predicted values and confidence band from GAM fit:
p <- predict(fit, newdata = data.frame(Bloom = rep(150, 40), 
                                       Wind =  x), se=T)
lines(x, p$fit + 2*p$se.fit, lty=2, col=2)
lines(x, p$fit - 2*p$se.fit, lty=2, col=2)

## 5e. Account for uncertainty in estimating edf
# Redo bootstrap while accounting for uncertainty in the "best" value for
# the degrees of freedom (edf), corresponding to the level of smoothing

# We re-fit the model without fixing df:
fit <- gam(log.surv ~ Bloom + s(Wind), data = pollock)
summary(fit)  # Note edf is close to but not exactly 2!

# The only difference between the new version of the function 'bf' below 
# and the earlier version is that the appropriate d.f. of the smooth 
# term is now estimated each time when fitting the model to a bootstrapped
# data set, rather than fixed at df=2 as in the above version:
bf <- function(resids, i, x) {
  dat$log.surv <- dat$fit + resids[i]
  new <- data.frame(Bloom = rep(150, 40), Wind = x)
  fit.boot <- gam(log.surv ~ Bloom + s(Wind), data = dat)
  predict(fit.boot, new)
}

### 5f. Run bootstrap with new bootstrap function:
pollock.boot2 <- boot(resid(fit), bf, x=x, R = 1999)

# Explore results as for 'pollock.boot' above