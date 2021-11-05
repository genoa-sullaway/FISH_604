########################################

# Module 10: Model selection & averaging
# Franz Mueter
# Last modified: 10-22-2019

########################################
library(here)
# Import and check data:
Branch <- read.csv("data/Branch sockeye.csv")
names(Branch)
str(Branch)
head(Branch)
# Create ordered factor for regime variable (for nicer plotting)
Branch$regime <- ordered(Branch$regime, levels = c("before", "after"))

# Attach new library (download first if needed) for computing
# small-sample Akaike Information Criterion and other libraries:
library(AICcmodavg)
library(MASS)
library(lattice)
library(ggplot2)
library(GGally)
library(mgcv)
library(lmtest)
library(nlme)
library(visreg)

##### Goal: Build a model to predict recruitment of sockeye salmon

###########
### Step 1: Examine response variable: Recruitment (R) 
dotplot(Branch$R)         
histogram(Branch$R)
# large outlier and/or long right tail suggests that a transformation
# may be needed
# Determine best Box-Cox transformation based on simple model of
# recruitment (R) as a linear function of spawner abundance (S):
boxcox(lm(R ~ S, data=Branch[-c(1:5),]))

# 95% confidence interval for "best" transformation includes 
#  a log-transformation, hence we assume that recruitment is
#  log-normally distributed, which is consistent with the idea
#  that salmon recruitment can be described by the Ricker model
#  (see Word doc)
#  Our response variable is therefore log(recruits-per-spawner):
Branch$logRS <- log(Branch$R / Branch$S)


###########
### Step 2: Examine predictor variables (check for multi-collinearity, examine for trends):
# Here's a different option for a pairwise scatterplot using ggplot2/GGoby:
# a bit slower than the others!)
ggpairs(Branch[,c(9,2,4:7)], lower=list(continuous='smooth'))
# Table of pairwise correlations:
cor(Branch[,c(9, 2,4:7)])
# No evidence of substantial collinearity! The largest correlation between
# explanatory variables is that between the PDO at two lags (0.37). 

win.graph(4,7)   # use 'quartz' on Mac OS
par(mfrow=c(4,1), omi=c(0.3,0,0.1,0), mar = c(0,4.5,0,1))
plot(SST.t2 ~ Year, data=Branch, type="b", xaxt="n", 
     ylab = "SST, lag 2", cex.lab = 1.6)
plot(SST.t3 ~ Year, data=Branch, type="b", xaxt="n", 
     ylab = "SST, lag 3", cex.lab = 1.6)
plot(PDO.t2 ~ Year, data=Branch, type="b", xaxt="n", 
     ylab = "PDO, lag 2", cex.lab = 1.6)
plot(PDO.t3 ~ Year, data=Branch, type="b", xaxt="n", 
     ylab = "PDO, lag 3", cex.lab = 1.6)
axis(1, cex.axis=1.5)
# There are considerable fluctuations in temperature and the
# PDO and a period of positive PDO values after the mid-60s.
# The lagged time series are of course identical, just shifted
# against each other by a year

dev.off()  # Close windows graphics device

###########
### Step 3: Conceptual model 
## Check for non-linearity:

# Examine relationship with SST:
# Linear model:
fit.lm <- lm(logRS ~ S + SST.t2 + SST.t3, data=Branch)
# Check for potential non-linearity in SST relationship:
fit.gam <- gam(logRS ~ S + s(SST.t2, k=5) + s(SST.t3, k=5), data=Branch)
summary(fit.gam)
plot(fit.gam, seWithMean=T)
AIC(fit.lm, fit.gam)
# AIC values are almost identical and the estimated degrees of freedom are very 
# close to 1 (SST at lag 2) or equal to 1 (lag 3), as is also evident from 
# the model fits, hence there is little evidence of non-linearity in the 
# relationship between log-survival and SST and we can simply use a linear
# regression

# Examine relationship with PDO:
fit.lm <- lm(logRS ~ S + PDO.t2 + PDO.t3, data=Branch)
fit.gam <- gam(logRS ~ S + s(PDO.t2, k=5) + s(PDO.t3, k=5), data=Branch)
summary(fit.gam)
plot(fit.gam, seWithMean=T)
AIC(fit.lm, fit.gam)
# The GAM model estimated one degree of freedom for each smoother, equivalent 
# to a linear model, hence the AIC values are identical and there is no
# evidence for non-linearity

# Thus we conclude that assuming linear relationships between the 
# environmental predictors and salmon recruitment is adequate!


###########
### Step 4: Constuct, fit, and evaluate the "global" model:

Branch.fits <- vector("list")  # Create empty list to hold output from all models

#--------------------------------------------------------------------------------
# Exercise 1: fit and examine the "full" or "global" model:
# Model A ("Global" model): SST, PDO and regime shift
#   ln(S/R) = a + b*S + c1*SST2 + c2*SST3 + cc1*PDO2 + cc2*PDO3 + d*regime + e

##  Fit the model and examine basic diagnostics. 
##  What do the model results suggest about the effects of spawner abundance, 
##    SST, PDO, and the regime shift on sockeye salmon survival?
##  Are any of the model assumptions violated?


# Solution:
Branch.fits$A <- lm(logRS ~ S + SST.t2 + SST.t3 + PDO.t2 + PDO.t3 + regime, 
  data=Branch)
summary(Branch.fits$A)

# Fit suggests a significant effect off SST during early ocean life that is very 
# consistent (very similar coefficients) at both lags, i.e. when age 1. or age 2. 
# smolts leave the river!

# No evidence of a significant PDO effect in addition to the 'local' SST effect,
# but there is evidence that log-survival changed at the time of the regime shift!

## Model diagnostics
par(mfrow=c(2,2))
plot(Branch.fits$A)
# Diagnostics suggest adequate model fit:
#   No large outliers
#   Approximately equal variances
#   No years with large Cook's distance, but one point with high 'leverage'(1956)
#--------------------------------------------------------------------------------

###### Before proceeding, make sure you have the correct model fit!!
# (you should have a multiple R-squared value of 0.5343)

# Identify single point with high leverage ('hat'), as evident in the
# plot of reisudals agains leverage:
hatvalues(Branch.fits$A)   # Leverage
# Observation 5 (1960) has high leverage:
Branch[hatvalues(Branch.fits$A) > 0.5, ]
# 1960 is a year with very high spawner abundance and high leverage!
# We could consider a transformation of spawner abundance, but that 
# would change the nature of the relationship 
# (which would no longer be a Ricker model!)

# Examine fit without high leverage point:
summary(update(Branch.fits$A, subset = -5), cor=F)
# Coefficients and standard errors are very similar and conclcusions
#  do not change, hence retain the 1960 data point!

## Plot resiuals against covariates and over time to check for patterns:
r <- rstandard(Branch.fits$A)  # Standardized residuals
par(mfrow=c(2,2), mar=c(4,2,1,1))
plot(Branch$SST.t2, r, xlab="SST"); abline(h=0, lty=2)
plot(Branch$PDO.t2, r, xlab="PDO"); abline(h=0, lty=2)
plot(Branch$S, r, xlab="Spawners"); abline(h=0, lty=2)
plot(Branch$Year, r, xlab="Year", type="b"); abline(h=0, lty=2)
# Plot of residuals against spawners shows high leverage point (S > 1200)
# Plot of residuals over time shows some evidence of autocorrelation (runs 
#  of positive and negative residuals)

visreg(Branch.fits$A)

## Tests for significant autocorrelation:
# 1. Examine residuals for autocorrelation at various lags:
acf(r)
# Some evidence of negative autocorrelation at lag 5, which 
#  approximately corresponds to the average generation length.
#  This implies that an unusually successful year class is 
#  followed by a poor year class 5 years later (possibly due to
#  density dependence - i.e. 'crowding' - in the stream)

# 2. Durbin-Watson test (requires 'lmtest'):
dwtest(Branch.fits$A)
dwtest(Branch.fits$A, alternative = "two.sided")
# Durbin-Watson test only tests for first-order autocorrelation!!!
# (See help file)

# 3. Compare models with and without 5th-order autocorrelation
#   (requires 'nlme')
fit.gls <- gls(logRS ~ S + SST.t2 + SST.t3 + PDO.t2 + PDO.t3 + regime, 
  data=Branch, method="ML", correlation = corARMA(p=5))
summary(fit.gls)
intervals(fit.gls)  # Note that confidence interval for Phi5 does NOT include 0

# Compare models using AIC:
AIC(fit.gls, Branch.fits$A)  # AIC suggests that AR(5) model fits better!
# Compare models using small-sample AIC (requires 'AICcmodavg')
AICc(fit.gls)                         
AICc(Branch.fits$A)

# While there is evidence for significant autocorrelation at lag 5 (based on AIC
# and significant phi value), the small-sample AIC criterion suggests that the
# improved fit does not justify the additional 5 phi parameters (although we could
# also fit the model with only 5th order autocorrelation)!

# Hence, for now, we assume no autocorrelation in the residuals (iid), but 
# we will return to this issue later!


###########
### Step 5: Fit set of alternative ('reduced') models:
# Model B: SST and PDO
# 	ln(S/R) = a + b*S + c1*SST2 + c2*SST3 + d1*PDO2 + d2*PDO3 + e
Branch.fits$B <- lm(logRS ~ S + SST.t2 + SST.t3 + PDO.t2 + PDO.t3, data=Branch, x=T)

# Model C: SST and regime shift
#   ln(S/R) = a + b*S + c1*SST2 + c2*SST3 + d*regime + e
Branch.fits$C <- lm(logRS ~ S + SST.t2 + SST.t3 + regime, data=Branch, x=T)

# Model D: SST
#	ln(S/R) = a + b*S + c1*SST2 + c2*SST3 + e
Branch.fits$D <- lm(logRS ~ S + SST.t2 + SST.t3, data=Branch, x=T)

# Model E: PDO
#	ln(S/R) = a + b*S + c1*PDO2 + c2*PDO3 + e
Branch.fits$E <- lm(logRS ~ S + PDO.t2 + PDO.t3, data=Branch, x=T)

# Model F: regime shift
#	ln(S/R) = a + b*S + d*regime + e
Branch.fits$F <- lm(logRS ~ S + regime, data=Branch, x=T)

# Model G: stock - recruit relationship only!
#	ln(S/R) = a + b*S + e
Branch.fits$G <- lm(logRS ~ S, data=Branch, x=T)



###########
### Step 6: Model selection

## 6a. Hypothesis testing

##################################################################################
# Exercise 2: Compare each of the models D, E, and F with model G (Ricker model) 
# using an F-test ('anova') and a likelihood-ratio test ('lrtest' in package 
# 'lmtest'). 

# Solution:
# Model D vs. G:
anova(Branch.fits$G, Branch.fits$D)
lrtest(Branch.fits$G, Branch.fits$D)
# Model E vs. G:
anova(Branch.fits$G, Branch.fits$E)
lrtest(Branch.fits$G, Branch.fits$E)
# Model F vs. G:
anova(Branch.fits$G, Branch.fits$F)
lrtest(Branch.fits$G, Branch.fits$F)

# Use a likelihood ratio test to compare models D and F with model C!

# Solution:
lrtest(Branch.fits$C, Branch.fits$D)
lrtest(Branch.fits$C, Branch.fits$F)

# What can you conclude about the effects of SST, the PDO, and the regime shift? 
# What is an obvious limitation of these tests?

# Results suggests that including SST (D) or a regime shift (F) significantly improved 
# the Ricker fit (Likelihood ratio test: Lambda=10.4, p = 0.0055 and Lambda=9.04, 
# p = 0.0026, respectively). In contrast, inclusion of the two PDO terms (E) did not
# lead to a significant improvement (p = 0.282) over the simple Ricker model, suggesting
# that we may eliminate PDO as an explanatory variable. Comparing model C to models D 
# and F suggests that the inclusion of both SST and a regime shift significantly 
# improved the fit, compared to models including only one of the variables (Lambda=7.15, 
# p=0.0075 and Lambda=8.50, p=0.014, respectively). In this example, the comparison of
# two subsets of nested models points to a single "best" model (C). In many cases, 
# however, hypothesis tests can suggest two or more potential candidate models that are
# not nested. The likelihood ratio test or the F-test cannot be used to compare 
# non-nested models (for example B vs. C or D vs. E), which is one of the shortcomings
# of hypothesis tests in the context of model selection. 

#--------------------------------------------------------------------------------


## 6b. Stepwise approach (AIC based)
step.best <- step(Branch.fits$A)
summary(step.best)
# 'Best'' model includes S, SST at lag 2, and regime shift!
# However, we should include SST at both lags for consistency!


## 6c. Information-theoretic approach:

## Illustration of computing AIC and AICc directly:
# Normal likelihood for model C, which can be computed from the number
# of observations (n), the number of paramters (k) and the residual sum
# of squares (RSS):
k <- 6   # Number of paramters (5 coefficients plus variance!)
n <- 36
RSS <- sum(resid(Branch.fits$C)^2)

# AIC = -2*log(L) + 2*k
n + n*log(2*pi) + n * log(RSS/n) + 2*k
# AICc = -2*log(L) + 2*k*(n/(n-k-1))
n + n*log(2*pi) + n * log(RSS/n) + 2*k*(n/(n-k-1))

# compare to:
AIC(Branch.fits$C)     
AICc(Branch.fits$C)

##################################################################################
# Exercise 3: Compute both the AIC values and the small-sample version (AICc) for 
#   all seven candidate models. 

# Solution: 
out <- matrix(NA, 7, 7)
dimnames(out) <- list(names(Branch.fits), c("k", "R2", "logL", "AIC", "delAIC", "AICc", "delAICc"))

for(i in names(Branch.fits)) {
  mod <- Branch.fits[[i]]
  logL <- logLik(mod)
  out[i, "k"] <-  attr(logL, "df")
  out[i, "R2"] <- summary(mod)$r.sq
  out[i,"logL"] <- logL
  out[i,"AIC"] <- AIC(mod)
  out[i, "AICc"] <- AICc(mod)
}
out

# Subtract minimum AIC (best model) from all values:
out[,"delAIC"] <- out[,"AIC"] - min(out[,"AIC"])
out[,"delAICc"] <- out[,"AICc"] - min(out[,"AICc"])
round(out,2)

# Which model is chosen as the AIC-best model? As the AICc-best model? 
# How strong is the support for other models?

# Solution: 
#  Model C is chosen as the best model in both cases, but the order of
#  the other models differs between the two selection criteria!

#--------------------------------------------------------------------------------

###########
### Step 7: Model diagnostics

##################################################################################
# Exercise 4: Thoroughly evaluate the 'best' model following the same approach that 
# we used for the 'global' model above

#   Does our best model (model C) provide an adequate fit to the data? 
#   Are all modeling assumptions met?

# Solution:

# See code for global model above! 
par(mfrow=c(2,2))
plot(Branch.fits$C)
# The same high leverage point is apparent but its residual is small
# and it has minimal influence on the fit (Cook's Distance close to 0)

r <- rstandard(Branch.fits$C)  # Standardized residuals
par(mfrow=c(2,2), mar=c(4,2,1,1))
plot(Branch$SST.t2, r, xlab="SST"); abline(h=0, lty=2)
plot(Branch$PDO.t2, r, xlab="PDO"); abline(h=0, lty=2)
plot(Branch$S, r, xlab="Spawners"); abline(h=0, lty=2)
plot(Branch$Year, r, xlab="Spawners", type="b"); abline(h=0, lty=2)
# Plot of residuals against spawners shows high leverage point (S > 1200)
# Plot of residuals over time shows some evidence of autocorrelation (runs 
#  of positive and negative residuals)

# Examine model fit:
visreg(Branch.fits$C)


#--------------------------------------------------------------------------------


###########
### Step 8: Model weights (relative weight of evidence for different models)
###         and multimodel inference (Model averaging)

# Akaike weights
D <- out[,"delAIC"]    
W.aic <- exp(-0.5*D) / sum(exp(-0.5*D))
# Small-sample Akaike weights
D <- out[,"delAICc"]
W.aicc <- exp(-0.5*D) / sum(exp(-0.5*D))

W <- round(cbind(W.aic, W.aicc),3)

# Table 4 in Word document:
Tbl4 <- data.frame(Model=row.names(out), out[,c(1,5)], Waic = W.aic, 
                   delAICc = out[,"delAICc"], Waicc = W.aicc)
Tbl4

## Model averaging:

# See Homework 4


# The package also provides a function to tabulate AIC values:
aictab(Branch.fits, names(Branch.fits))          # small-sample AICc (default)
aictab(Branch.fits, names(Branch.fits), sec=F)   # AIC

# Note that the function sorts the models from "best" (smalles AIC/AICc) to "worst"


################################################################################### 
### Appendix 1: Compute bootstrap model selection probabilities for 
###             sockeye salmon models for Branch River data
################################################################################### 

if(F)   # Prevents bootstrap from running when sourcing file
{
  
B <- 10000     # Number of bootstrap samples!
# We chose a large number to more accurately estimate weights
#  for the poor models that are rarely chosen as 'best'

dat <- Branch[,c(9, 2, 4:8)]
dat$regime <- as.numeric(dat$regime)
dat <- as.matrix(dat)
n <- nrow(dat)
dimnames(dat)

# set up matrices to store output for Ricker-a parameter (intercept),
# and for AICc values:
a <- matrix(NA, B, 7)
.AICc <- matrix(NA, B, 7)

dimnames(a) <- dimnames(.AICc) <- list(NULL, LETTERS[1:7])

# Constant to add to likelihood term (does not depend on model)
.const <- n + n*log(2*pi)

for(i in 1:B) {
# draw bootstrap sample:
  j <- sample(1:n, replace = T)
	x <- dat[j,]  # Select re-sampled rows 
	y <- x[,1]    # select response
	x <- x[,-1]   # all predictors
# model A:
	fit <- lsfit(x, y)
	a[i, "A"] <- fit$coef[1]
	k <- length(fit$coef) + 1
	.AICc[i, "A"] <- .const + n * log(sum(fit$resid^2)/n) + 2*k*(n/(n-k-1))
# model B:
	fit <- lsfit(x[,-6], y)
	a[i, "B"] <- fit$coef[1]
	k <- length(fit$coef) + 1
	.AICc[i, "B"] <- .const + n * log(sum(fit$resid^2)/n) + 2*k*(n/(n-k-1))
# model C:
	fit <- lsfit(x[,-c(4,5)], y)
	a[i, "C"] <- fit$coef[1]
	k <- length(fit$coef) + 1
	.AICc[i, "C"] <- .const + n * log(sum(fit$resid^2)/n) + 2*k*(n/(n-k-1))
# model D:
	fit <- lsfit(x[,-c(4,5,6)], y)
	a[i, "D"] <- fit$coef[1]
	k <- length(fit$coef) + 1
	.AICc[i, "D"] <- .const + n * log(sum(fit$resid^2)/n) + 2*k*(n/(n-k-1))
# model E:
	fit <- lsfit(x[,-c(2,3,6)], y)
	a[i, "E"] <- fit$coef[1]
	k <- length(fit$coef) + 1
	.AICc[i, "E"] <- .const + n * log(sum(fit$resid^2)/n) + 2*k*(n/(n-k-1))
# model F:
	fit <- lsfit(x[,-(2:5)], y)
	a[i, "F"] <- fit$coef[1]
	k <- length(fit$coef) + 1
	.AICc[i, "F"] <- .const + n * log(sum(fit$resid^2)/n) + 2*k*(n/(n-k-1))
# model G:
	fit <- lsfit(x[,-(2:6)], y)
	a[i, "G"] <- fit$coef[1]
	k <- length(fit$coef) + 1
	.AICc[i, "G"] <- .const + n * log(sum(fit$resid^2)/n) + 2*k*(n/(n-k-1))
}

# identify AICc-best model for each bootstrap sample:
# 1. remove smallest AICc from each row:
delAIC <- t(apply(.AICc, 1, function(x) x - min(x)))

# Count number of times each model was chosen as "best" model:
Best <- apply(delAIC, 2, function(x) sum(x==0))
# Proportion of times each model was chosen as best:
Best/B

# Compare to Akaike weights!
}
