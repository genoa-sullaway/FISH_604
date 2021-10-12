##########################################################

# Module 7: Linear Models 

# ANCOVA example & pollock data (for diagnostics)

# Franz Mueter
# Last modified: 9-20-2021
library(here)
##########################################################


####### Analysis of Covariance (ANCOVA) example

# Read in sex, length, weight data for sablefish:
sable <- read.csv("data/sablefish2.csv", as.is=F)
sable
names(sable)

#### Graphical analysis
# Compare male and female weights (overall):
plot(Weight~ Sex, data=sable)

# Weight at length by sex
plot(Weight ~ Length, data=sable, subset=Sex=="M", col=4, pch=16, 
	ylim=range(Weight),cex.axis=1.5,cex.lab=1.5)
points(Weight ~ Length, data=sable, subset=Sex=="F", col=2, pch=16)

# allometric weight-length relationship (W = a*L^b) is linear
# on the log-log scale:
# log-weight vs log-length by sex
plot(log(Weight) ~ log(Length), data=sable, subset=Sex=="M", col=4, pch=16, 
	ylim=range(log(Weight)),cex.axis=1.5,cex.lab=1.5)
points(log(Weight) ~ log(Length), data=sable, subset=Sex=="F", col=2, pch=16)
# Add linear regression lines fit separately by sex:
abline(lm(log(Weight) ~ log(Length), data=sable, subset=Sex=="F"), col=2)
abline(lm(log(Weight) ~ log(Length), data=sable, subset=Sex=="M"), col=4)

# Or, using ggplot:
library(tidyverse)
# Weight vs Length
ggplot(sable) + 
  geom_point(aes(x=Length, y = Weight, color=Sex))

# log-weight vs log-length by sex with separate regression lines:
ggplot(sable, aes(x=log(Length), y = log(Weight), color=Sex)) + 
  geom_point() +
  geom_smooth(method="lm", se=F)

# log-weight vs log-length by sex with same slope, different intercept:
sable <- mutate(sable, logLength=log(Length), logWeight=log(Weight))
fit <- lm(logWeight ~ Sex + logLength, data=sable)
ggplot(sable, aes(x=log(Length), y = log(Weight), color=Sex)) + 
  geom_point() +
  geom_line(data=fortify(fit), aes(x=logLength, y=.fitted))


##### Fit and compare different models

# One-way ANOVA (Without covariate):
fit1 <- lm(log(Weight) ~ Sex, data=sable)
summary(fit1)   
# No difference in average weight between males and females

# ANCOVA model:
# Assume separate intercept, ONE slope (same for males and females)
fit2 <- lm(log(Weight) ~ Sex + log(Length), data=sable)
summary(fit2)

# Fit separate slope and intercept for males and females:
fit.F <- lm(log(Weight) ~ log(Length), data=sable, subset=Sex=="F")
fit.M <- lm(log(Weight) ~ log(Length), data=sable, subset=Sex=="M")
summary(fit.F)
summary(fit.M)

#### Model comparisons
# General modeling approach: Fit model with separate slope and intercept
# and compare to model with separate intercepts only (same slope).

# First, we fit a model with separate slopes & intercepts by sex. 
# This is the same as the two separate models fit by sex above, 
# EXCEPT that a single error variance is estimated here, whereas 
# separate error variances are implied above (fit.F, fit.M)
fit3 <- lm(log(Weight) ~ Sex / log(Length), data=sable)
summary(fit3)

# Equivalence of nesting and interaction
# We can fit the same exact model using a different parameterization,
# where the parameters have a different interpretation:
fit4 <- lm(log(Weight) ~ Sex * log(Length), data=sable)
summary(fit4)
# The coefficients now correspond to the intercept and slope ('log(Length)') 
# for female fish and to the differences in intrecept (SexM) and slope 
# (SexM:log(Length) between males & females

# The model corresponding to 'fit2' is nested within fit 3 and fit4.
# This is most obvious for fit4, which simply has one additional 
# coefficient (the interaction term) compared to fit2. 

# We can compare these nested models in two (equivalent) ways:

# 1. Wald test 
# The interaction term in fit4 is highly significant
# based on the t-statistic (p = 5.8 * 10^-5)
summary(fit4)

# 2. F-test 
# We can compare any two nested models using an F-test. Because the 
# difference between fit2 and fit4 involves a single parameter,
# the two tests are functionally identical in this case (in fact, the 
# F-statistic is simply the square of the t-statistic as you can 
# easily verify:
-4.134^2  # Squared t-statistic
anova(fit2, fit4, test="F")  # Note F-statistic
# The Null hypothesis for an F-test, in general terms, is that the 
# "smaller" model (fit2 in our case) is the "true" model and that the
# additional parameters in the more complex model (here: fit4) do not
# reduce the residual sum of squares significantly (review 'F-tests'
# if this is not familiar) 

# Based on these tests, we conclude that both the intercept and slope
# in the log(Weight) - log(Length) relationship differ significantly
# between male and female sablefish



#######################################################################
# Pollock example for diagnostics
# Read in data:
pollock <- read.csv("pollock.csv")
row.names(pollock) <- pollock$Year

WP.fit <- lm(log.surv ~ Bloom + poly(Wind,2), data=pollock)
plot(WP.fit, 1)
plot(WP.fit, 2, id.n = 5)
plot(WP.fit, 3)
plot(WP.fit, 4)
plot(WP.fit, 5)
plot(WP.fit, 6, col=4, pch=16)
influence.measures(WP.fit)
influence(WP.fit)



