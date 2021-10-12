##############################################

# Lab 7:  Generalized Linear Modeling
# Franz Mueter

# Last modified: 10/3/2019
# 
#####    Poisson regression
#####    Example: Density of salamanders

##############################################

# import data
salamander<-read.csv("salamander.csv")
names(salamander)
attach(salamander)
par(mfrow = c(2, 2))
plot(PCTCOVER, SALAMAN)
plot(FORESTAGE, SALAMAN)
plot(PCTCOVER, log(SALAMAN + 1))
plot(FORESTAGE, log(SALAMAN + 1))

# Fit a LOWESS smooth function on the un-tansformed scale and on log-transformed scale:
scatter.smooth(PCTCOVER, SALAMAN)		
scatter.smooth(FORESTAGE, SALAMAN)		
scatter.smooth(PCTCOVER, log(SALAMAN+1))		
scatter.smooth(FORESTAGE, log(SALAMAN+1))		
detach(salamander)

# create dummy variables of low (< 70%) vs. high (>70%) cover 
salamander$cover <- salamander$PCTCOVER > 70	
salamander$cover				
fit <- glm(SALAMAN ~ (PCTCOVER * FORESTAGE + I(PCTCOVER^2) +  
 	I(FORESTAGE^2)) * cover, data = salamander, family = poisson) 
summary(fit)

# goodness-of-fit test:
chisq <- sum(residuals(fit, type="pearson")^2)
chisq
pchisq(chisq, df.residual(fit), lower.tail = FALSE) 
# Alternatively (using deviance residuals):
pchisq(deviance(fit), df.residual(fit), lower.tail = FALSE) 

# Check for overdispersion
r <- residuals(fit)
y.hat <- fitted(fit)
par(mfrow=c(1,1))
plot(y.hat,r)
abline(h=c(-2,0,2),col=c(3,2,3))
# The plot suggests overdispersion (many residuals in the tails)

#---------------------------------------------------------------------------------
# Exercise 1: Estimate overdispersion:



#---------------------------------------------------------------------------------




#---------------------------------------------------------------------------------
### Exercise 2: Re-fit model without forest age
# Fit the same model as above ('fit') without any of the FORESTAGE terms:

#      fit2 <-   [Your call to 'glm' here]

# Compare the two models using 'anova()' using a chisq-test (test = "Chisq" argument to anova).
# See ?anova.glm for details!

# Based on the results, does it seem justified to drop forest age from the model 
# (or conversely, is the reduction in deviance enough to justify the extra parameters
# of the more complex model that includes forest age)? 



#---------------------------------------------------------------------------------


# Quasi-likelihood approach 
fit3 <- update(fit2, family = quasipoisson)

#---------------------------------------------------------------------------------
### Exercise 3: Overdispersion parameter
# Look at the model results produced by summary. Look for the estimate of the 
# overdispersion parameter and compare to the estimate obtained above! Compare
# the estimated coefficients and standard errors produced by 'fit2' and 'fit3'
# (using 'summary()')! How and why do they differ?
 

#---------------------------------------------------------------------------------


# Diagnostics:
par(mfrow=c(2,2))
plot(fit3, which=1:4)

# Predict values with normal confidence bands
par(mfrow=c(1,1))
new <- data.frame(PCTCOVER = 0:100, cover = c(0:100) > 70)
y.hat<-predict(fit3,new,se=T, type="response")
plot(SALAMAN~PCTCOVER,data=salamander)
lines(y.hat$fit~new$PCTCOVER,lwd=3)
lines((y.hat$fit-1.96*y.hat$se.fit)~new$PCTCOVER,col=3)
lines((y.hat$fit+1.96*y.hat$se.fit)~new$PCTCOVER,col=3)

# Predict values and confidence bands on original scale, back-transformd
y.hat<-predict(fit3,new,se=T)
j <- new$PCTCOVER <= 70
plot(SALAMAN~PCTCOVER,data=salamander)
lines(new$PCTCOVER[j], exp(y.hat$fit[j]), lwd=3)
lines(new$PCTCOVER[!j], exp(y.hat$fit[!j]), lwd=3)
lines(new$PCTCOVER[j], exp(y.hat$fit-1.96*y.hat$se.fit)[j],col=3)
lines(new$PCTCOVER[!j], exp(y.hat$fit-1.96*y.hat$se.fit)[!j],col=3)
lines(new$PCTCOVER[j], exp(y.hat$fit+1.96*y.hat$se.fit)[j],col=3)
lines(new$PCTCOVER[!j], exp(y.hat$fit+1.96*y.hat$se.fit)[!j],col=3)
# Note that the confidence band for PCTCOVER over ~50% goes through
# the roof because there are no observations. Because of the two-piece fit,
# the fitted values and the confidence intervals are discontinuous

# Predicted values and confidence bands on log-scale
plot(new$PCTCOVER[j], y.hat$fit[j],	xlim=c(0,100), ylim=c(-11.5, 7.5), 
    type="l", lwd=3)
lines(new$PCTCOVER[!j], y.hat$fit[!j], lwd=3)
lines((y.hat$fit-1.96*y.hat$se.fit)[j]~new$PCTCOVER[j],col=3)
lines((y.hat$fit-1.96*y.hat$se.fit)[!j]~new$PCTCOVER[!j],col=3)
lines((y.hat$fit+1.96*y.hat$se.fit)[j]~new$PCTCOVER[j],col=3)
lines((y.hat$fit+1.96*y.hat$se.fit)[!j]~new$PCTCOVER[!j],col=3)


# Fit the above model to plots with canopy cover over 70% only,
# using the negative binomial distribution and a log link for
# comparison with the quasipoisson
library(MASS)
fit4 <- glm.nb(SALAMAN ~ PCTCOVER + I(PCTCOVER^2), data = salamander,
   link="log", subset = cover) 
# 'cover' (logical variable) is used to extract plots with PCTCOVER > 70
summary(fit4)
anova(fit4, test = "Chi")

# residual diagnostics
par(mfrow=c(2,2))
plot(fit4, which=1:4)


# Because the Poisson is a special case of the negative binomial model, 
# We can directly test if the negative binomial model fits better using
# a likelihood ratio test (compute "by hand" or use 'lrtest')

# First, fit the same model assuming a Poisson distribution: 
fit4b <- glm(SALAMAN ~ PCTCOVER + I(PCTCOVER^2), data = salamander,
   family = poisson, subset = cover) 

# Compute likelihood ratio test statistic:
L1 <- -logLik(fit4b)   # minus negative log-Likelihood, Model 1
L2 <- -logLik(fit4)    # minus negative log-Likelihood, Model 2
q <- attr(L1, "df")   # model degrees of freedom are stored as attribute of the log-Likelihood
p <- attr(L2, "df")
LRT <- 2*L1 - 2*L2   # Likelihood-ratio test statistic!

# Under the null hypothesis (no difference between models), the distrubtion of 
# LRT follows a Chi-square distribution with p-q degrees of freedom
# The probability that a Chi-square with p-q d.f. is larger than LRT is given by:
1-pchisq(LRT, p-q)    # where pchisq computes the cumulative distribution function

# Graphical illustration:
x <- seq(0, 30,by=0.1)
plot(x, dchisq(x, p-q), type="l", lwd=2); abline(h=0)
abline(v=LRT, col=2, lwd=2)
# The probability that a Chi-square with 1 d.f. is as large as LRT=23.32 corresponds
# to the area under the curve to the right of the red line, which is vanishingly small
# Hence we can reject the null hypotheses and conclude that the negative binomial
# model is clearly the better model!

# As a simple alternative we can use 'lrtest' to get the same test:
library(lmtest)
lrtest(fit4, fit4b)


################################### Extra:
# The same models can be fit using function 'gamlss' in 
# the package of the same name.
# This function can also fit generalized linear models
# (as well as generalized additive models)
# and you can chose from a very large number of distributions
# (see help for 'gamlss.family')
library(gamlss)
fit1 <- gamlss(SALAMAN ~ PCTCOVER + I(PCTCOVER^2), 
    data = salamander[salamander$cover,], family = NBI)
summary(fit1)
fit2 <- gamlss(SALAMAN ~ PCTCOVER + I(PCTCOVER^2), 
    data = salamander[salamander$cover,], family = PO)
summary(fit2)
lrtest(fit1, fit2)

