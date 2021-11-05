##################################################

# Lab 8:    Generalized Linear Model
#           Zero-inflated models

# Franz Mueter
# Last modified: 10-08-2021

##################################################
library(tidyverse)
library(here)

#### Chum salmon example

# Nenana River chum salmon catches
# Data courtesy of Parker Bradley

# The data consist of the counts of chum salmon (and 2 other species)
# collected in two types of nets (Fyke net, plane trap) along two sides
# of the Nenana River (river left (L) and river right (R)), as well as
# in mid-channel (M) during the summer of 2011

# Our goal for this exercise is to estimate the trend in mean count 
# of chum salmon over the course of the season to identify any
# seasonal patterns of migration

# Import and examine data:
Nenana <- read.csv("data/Nenana fish.csv")
head(Nenana)

# Convert dates to Julian day for plotting and modeling 
Nenana$Julian <- strptime(Nenana$Date,"%m/%d/%Y")$yday

library(ggplot2)
# CPUE by Julian day (overall, by gear and by location in river):
(p1 <- ggplot(Nenana, aes(Julian, Chum)) + geom_point())
p1 + facet_wrap(~Gear)
p1 + facet_wrap(~Side)

# For this analysis, we restrict data to early season to eliminate 
# excessive zeros(only two fish caught after July 10 (day 190))
Chum <- Nenana[Nenana$Julian <= 190, c(1:5,8)]
head(Chum)
(p2 <- ggplot(Chum, aes(Julian, Chum)) + geom_point())
p2 + facet_wrap(~Gear)
p2 + facet_wrap(~Side)

p3 <- ggplot(Chum) + geom_histogram(aes(Chum), binwidth=5)
p3
p3 + facet_wrap(~Side)
# The histograms is very skewed towards small counts (right-skewed)

# For simplicity, We focus on a single sampling site 
# on the right side of the river (R). First, we'll save 
# the subset as a separate data frame for convenience
Chum.R <- Chum[Chum$Side == "R",]   
table(Chum.R$Gear)  # Note that there are only fyke net samples, which
                    # avoids the need to account for gear differences

ggplot(Chum.R, aes(Julian, Chum)) + geom_point() + 
  ggtitle("Chum salmon counts, Nenana River, by Julian day, river right")


#-------------------------------------------------------------------------
# Exercise 1:  Examine the number of zero counts in 'Chum.R' 
# Use 'table' or 'barplot' or directly compute the proportion of zeros in the data
# (both overall and by sampling site ('Side')



# Does there seem to be an excess of zeros?


# Use a binomial model with 'Julian' as the only predictor variable to 
# examine the trend in the proportion of positive catches over time!

# Use 'Chum > 0' as your response variable in the model formula!

# Thus your response will consist of zeros and ones depending on whether any
# chum were caught (1) or not (0):
with(Chum.R, Chum > 0)

fit <- glm(...)

# With your fitted model ('fit'), you should be able to plot the
# estimated trend in the proportion of positive catches over time:
plot(Chum.R$Julian, fitted(fit, type="response"), type="l")
# Connecting the fitted values with a line works just fine here 
# because the Julian days are ordered

#-------------------------------------------------------------------------


### We first model the trend in counts over time using a
# Poisson regression with log link 
# (For the moment, we won't worry about extra zeros!)
fit1 <- glm(Chum ~ Julian, data = Chum.R, family = poisson)
summary(fit1)

# Goodness of fit:
# The null and residual deviances are (sort of) equivalent to 
# the total sum of squares and residual sum of squares in a 
# linear model. There is no R^2 in GLM models, but we 
# can get an equivalent that's close, which is the proportion 
# of deviance explained, defined as:

#                  null deviance - residual deviance
#           100 X ------------------------------------
#                          null deviance

null.dev <- summary(fit1)$null.deviance
resid.dev <- deviance(fit1)

100 * (null.dev - resid.dev) / null.dev
# About 65% of the variability in the counts (on the log scale)
# is accounted for by the estimated trend

### Visualizing the fitted values
plot(Chum ~ Julian, data=Chum.R, col=2)

# We can compute predicted values on the log-scale for each 
# sampling day using the estimated coefficients.
(cf <- coef(fit1))
y.hat <- cf[1] + cf[2] * Chum.R$Julian
# predicted value = 15.7920782 - 0.0934361 * Julian

# Remember that the Poisson model is fit on the scale of 
# the log-transformed response (using the 'log-link').
# Therefore, to back-transform to the scale of the response
# variable we exponentiate the predicted values and then 
# add add them to the plot:
lines(Chum.R$Julian, exp(y.hat))
# The curvature in the fitted line is directly related to the 
# fact that we are using a log-link (this fit is the equivalent 
# of a simple linear regression in a linear model)

# Alternatively, we can simply use the 'predict'' function to  
# compute predicted values on the scale of the response:
predict(fit1, type="response") 
# The values are identical to the back-transformed y.hat
# values that we computed above!

### We generally want to acknowledge and ideally visualize 
#   uncertainty in the fitted regression line, hence let's
#   add 95% confidence bands:
# We first compute CIs on the log-scale:
Fitted <- predict(fit1, type="link", se=T)
names(Fitted)    
# The result is a list with three components:
#   'fit'  is the fitted values (on the scale of the link function)
#   'se.fit' are the standard errors
#   'residual scale'   is the scale parameter in the GLM, which
#                      is assumed to be 1 for the Poisson

# Exponentinating the fitted values gives us y.hat' (again)
exp(Fitted$fit)  
# We can compute 95% confidence intervals as usual, on the scale 
# of the linear predictor (log scale):
lwr <- Fitted$fit - 1.96 * Fitted$se.fit
upr <- Fitted$fit + 1.96 * Fitted$se.fit
# To obtain correct confidence intervals on the back-transformed
# scale, we also need to exponentiate the lower & upper bounds: 
lower <-  exp(lwr)
upper <-  exp(upr)

# Let's add them to the plot:
lines(Chum.R$Julian, upper, lty=2)
lines(Chum.R$Julian, lower, lty=2)
# Voila! A simple Poisson regression with a 95% confidence band

# Of course, we can do all that in a simple call to 'visreg':
visreg(fit1, scale="response")
# Verify that the confidence bands are the same as computed above:
lines(Chum.R$Julian, upper, lty=2)
lines(Chum.R$Julian, lower, lty=2)


## Accounting for fishing effort (as measured by volume filtered)

# In the analysis so far we assumed that counts are comparable, 
# i.e. that each count represents the same amount of effort. 
# However, each set filtered a different volume of water. 
# We can account for the volume filtered by using an 'offset':

#     log(catch/volume)  = log(catch) - log(volume) = a + b*Julian
#             log(catch) = a + b*Julian + log(volume)
# where 'log(volume)' (which has no regression coefficient) is the 'offset'

fit2 <- glm(Chum ~ Julian + offset(log(Volume)), family=poisson, data = Chum.R)
summary(fit2)
AIC(fit1, fit2) # The offset does improve the model (as expected!)

# Model diagnostics
plot(fit2)
# Note the clear heteroscedasticity evident in the first plot, which is likely related to 
# overdispersion, which we will examine below!
# There is one very influential point (the observation labeled "2"):
Chum.R["2",]   # This corresponds to the largest observed count!
max(Chum.R$Chum)

# The predicted values are a little mode tricky now as they depend on the
# amount of water filtered. Therefore, we can only visualize a smooth trend 
# in catch-per-unit-effort (CPUE) for a given amount of water filtered. 
plot(Volume ~ Julian, data=Chum.R)
abline(h=mean(Chum.R$Vol), col=4, lty=2)
# Note that volume filtered decreased over time, which may explain 
# part of the decline in observed catches!

# Instead of the raw counts, we can plot counts adjusted to mean volume
# to show the decline in standardized CPUE over time:
plot(Chum/Volume * mean(Volume) ~ Julian, data=Chum, col=2)

# To visualize predicted values, we first need to compute predicted values
# based on the mean volume filtered. We set up a data frame that contains
# a range of Julian days and the mean Volume:
dat <- data.frame(Julian = 130:190, Volume = mean(Chum.R$Vol))
# Predicted catches assuming that the volume filtered is constant (at the mean) 
Fitted <- predict(fit2, newdata=dat, se=T)
mu <- exp(Fitted$fit)  # predicted values on response scale
lower <-  exp(Fitted$fit - 1.96 * Fitted$se.fit)
upper <-  exp(Fitted$fit + 1.96 * Fitted$se.fit)
x <- dat$Julian
lines(x, mu)
lines(x, upper, lty=2)
lines(x, lower, lty=2)

# Evaluating Poisson fit:
summary(fit2)
# The counts are clearly overdispersed (Compare residual deviance to 
# residual d.f., which should be roughly the same if the Poisson 
# assumption is valid). 
# To account for the apparent overdispersion, we fit a negative binomial
# model. 'glm' does not currently accommodate the negative binomial, hence
# we use function 'glm.nb' from the package MASS, which does:
library(MASS)
fit3 <- glm.nb(Chum ~ Julian + offset(log(Volume)), data = Chum.R)
summary(fit3)
# Note that the residual deviance is now very similar to the d.f.
# and that there is still a highly significant decrease in CPUE 
# over time (no surprise there!)

# Model diagnostics:
plot(fit3,1:5)
# Note that the distribution of the standardized residuals 
# looks much better and there is no overly influential outlier!
# The largest standardized residual ("119"), which also has 
# the largest Cook's Distance now corresponds to a count of 8 
# that occurs later in the time series (when the mean count 
# is very low)

# While the difference in the models here is striking and the
# negative binomial results in an obvious improvement over the
# Poisson based on visual inspection of residuals and comparison
# of the deviances, we can also conveniently compare the negative 
# binomial model fit and the Poisson fit, which is nested in the
# negative binomial model, using a likelihood ratio test: 

# Unfortunately, our usual approach using the 'anova' function
# does not work here as it does not know how to handle the
# 'glm.nb' output:
anova(fit2,fit3, test="LRT")
# You can compute the LRT by hand as we've done before (using
# 'logLik' to extract the likelihood from each model), or we 
# can use a convenient function in the 'lmtest' package: 
library(lmtest) 
lrtest(fit2, fit3)  
# Based on the reduction in deviance achieved by adding a single 
# parameter (the 'theta' parameter of 'glm.nb', the negative binomial
# model provides a better fit, and the difference is highly significant.
# This is consistent with our earlier observation that the Poisson 
# counts are highly over-dispersed!

# A considerable difference in the fits is also readily apparent when
# we compare the AIC values:
AIC(fit2, fit3)   
# Here lower is better, and the substantial reduction in AIC 
# (corresponding to a substantial difference in log-likelihoods)
# provides strong support for the negative binomial model.


#-------------------------------------------------------------------------
# Exercise 2: Follow the template for fit2 above to plot counts
#   over time (adjusted to mean volume) and superimpose the predicted
#   values from the NB model!


# Compare the confidence bands of the two models!


#-------------------------------------------------------------------------


### Zero-inflated models
# We are now ready to address the potential for extra zeros in the 
# data, which can also contribute to the overdisperion problem

# Recall the large number of zeros:
hist(Chum.R$Chum,nclass=20)
# The zeros are more clearly evident in a barplot that shows 
# the number of sets that caught a given number of fish:
x <- barplot(tabulate(Chum.R$Chum+1))
axis(1, at=x, lab=0:max(Chum.R$Chum))

# Note that most of the zeros occur later in the season when 
# the mean counts are very low, hence many of them may be a 
# natural part of the NB distribution, which has more zeros 
# when the mean count is low (rather than reflecting "extra" zeros)

# Whether there are extra zeros or not is difficult to evaluate, 
# but we can fit a hurdle model or a zero-inflated model 
# and compare the fits:

# We will use the 'pscl' package, which contains functions for fitting 
# hurdle models and zero-inflated models
if(!("pscl" %in% installed.packages())) {
  install.packages("pscl")
}
library(pscl)

# We fit a hurdle model that uses a negative binomial model for
# the counts, a binomial model for modeling the probability of 
# zeros (which is the default, with logit link) and we use an 
# offset as before. The same variables ('Julian' and the offset)
# are included as predictor variables for both the binomial and
# NB components (by default):
fit4 <- hurdle(Chum ~ Julian + offset(log(Volume)), 
               dist="negbin", data = Chum.R)
summary(fit4)
# Note that the output contains two sets of coefficients, 
# one for the hurdle component and one for the truncated 
# count component (modeled via NB)! The NB component includes 
# an estimate of the theta parameter (on the log-scale). Note that
# both coefficients for Julian day are negative and significantly
# different from zero, implying that there is a signficant decrease
# over time in both the probability of catching chums (fewer 
# 'successes' as the season progresses) and in the number of
# chum salmon being caught.

# Let's visualize the 'hurdle' component first (that is, the 
# proportion of positive catches):
d <- 130:190
p.logit <- 17.97077 -0.10738 * d
plot(Chum > 0 ~ Julian, data = Chum.R, pch="|", cex=2)
lines(d, exp(p.logit) / (1 + exp(p.logit)), type="l", lwd=3)

# The output also contains the estimate (plus standard errors) 
# of the overdispersion parameter 'theta' of the negative 
# binomial:     Var(y) = mu^2 + mu^2 / theta
# where a larger theta implies a smaller variance
fit4$theta

#-------------------------------------------------------------------------
# Exercise 3: What is the amount of "extra variance" relative 
#   to the Poisson variance (mu^2) at mean counts of 5 and 10? 

#   As per the above formula, the variance depends on both the
#   mean and on the 'theta' parameter.


#-------------------------------------------------------------------------


plot(Chum/Volume * mean(Volume) ~ Julian, data=Chum, col=2)

# Compute predicted values based on the average volume filtered:
dat <- data.frame(Julian = 130:190, Volume = mean(Chum.R$Vol))
# Currently, the 'predict' function for hurdle models does not compute
# standard errors, hence we only plot the fitted values
# (We could compute standard errors, but it's a bit more work!)
Fitted <- predict(fit4, newdata=dat, type="response")
lines(dat$Julian, Fitted, col=4)

# Compare the fit of the negative binomial without extra zeros 
# to the hurdle model:
AIC(fit3, fit4)

### Which model fits better according to the AIC criterion?

# Note that we can't get simple diagnostic plots either, 
# but we can extract residuals' to look at some relevant 
# diagnostic plots ('Deviance' residuals are not available, 
# hence we use 'Pearson' residuals)
r <- resid(fit4, type="pearson")

# Residuals versus fitted values:
plot(fitted(fit4), r); abline(h=0, lty=2, col=4)
# No strong heteroscedasticity, but 2 fairly large outliers (> 3)

hist(r)  # Outliers are also evident in the histogram

#-------------------------------------------------------------------------
# Exercise 4: Read the help file for 'zeroinfl' and fit an 
#              equivalent model using a zero-inflated 
#              negative binomial distribution. 

# Compare the zero-inflated model to the hurdle model!

# Does the model suggest a signficiant amount of "zero-inflation"?
# Why or why not? 

?zeroinfl
#-------------------------------------------------------------------------


# So far, we have not worried about the time series nature of the data but it 
# is very likely that there may be temporal autocorrelation in the counts


#-------------------------------------------------------------------------
# Exercise 5: Plot the residuals from 'fit4' over time 
# (Julian day) and connect residuals with a line 
# (which makes it easier to spot time trends)


# Is there evidence of autocorrlation?


#-------------------------------------------------------------------------

# We will tackle autocorrelation in a future lab!

