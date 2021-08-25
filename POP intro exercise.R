#######################################################################
# 
# FISH 604: Modern Applied Stats for Fisheries
#
# Module 1 - Introduction
#
# Import & plot Pacific Ocean Perch biomass series for in-class exercise
#
# Written by Franz Mueter
# Last modified: 7-30-2021
#
#######################################################################

# Import and view data:
(POP <- read.csv("Pacific Ocean Perch.csv"))

# Read in the function 'error.bar' from script for plotting error bars
source("ErrorBarFunction.R")

# Using POP data frame, plot estimates with 95% confidence intervals:
par(mar=c(2,5,1,1))
with(POP, error.bar(Year, Biomass/1000, (Biomass-1.96*StdErr)/1000, (Biomass+1.96*StdErr)/1000, F,
          xlab = "", ylab = "Survey biomass (1,000t)", cex.lab=1.8, cex.axis = 1.4,
          pch=16, cex=1.5, col=4))

# There are a number of ways in which we could estimate the increase in POP 
# biomass between 1984 and 2013. Here are a few:

#######################################################################
#### 1. Estimate based on difference between last and first year:
#  (in 1000s of tons)
POP.sub <- subset(POP, Year %in% c(1984, 2013))
POP.sub
with(POP.sub, error.bar(Year, Biomass/1000, (Biomass-1.96*StdErr)/1000, (Biomass+1.96*StdErr)/1000, F,
                    xlab = "", ylab = "Survey biomass (1,000t)", cex.lab=1.8, cex.axis = 1.4,
                    pch=16, cex=1.5, col=4))
abline(h=POP.sub$Biomass/1000, lty=2, col=2)
# Observed difference:
(diff <- round(POP[POP$Year==2013, "Biomass"] - POP[POP$Year==1984, "Biomass"]))

### Accounting for uncertainty
## a. Analytical solution using variance of difference (see Module 2):
#   var(B_2013 - B_1984) = var(B_2013) + var(B_1984)
# This equation for the variance of a difference is true if the
# two observations (B_1984 and B_2013) are independent. Given they
# were taken 29 years apart, that is a fair assumption.
Var.diff <- sum(POP.sub[,"StdErr"]^2)
# Standard error of the mean (= square root of variance):
sqrt(Var.diff)
# 95% confidence interval (assuming normality):
ci <- round((diff + c(-1.96, 1.96) * sqrt(Var.diff))/1000)  # in thousands of tons
ci

## b. Bootstrapping (assuming normal distribution for mean biomass):
POP.sub
B.1984 <- rnorm(10000, mean=POP.sub[1,2], sd=POP.sub[1,3]) 
B.2013 <- rnorm(10000, mean=POP.sub[2,2], sd=POP.sub[2,3]) 
# Estimated difference and its uncertainty:
hist(diff.boot <- (B.2013 - B.1984)/1000, col=5, main="Difference in biomass")  # in 1000 t
# 95% confidence interval:
(ci.boot <- quantile(diff.boot,c(0.025, 0.975)))
abline(v=ci.boot, col=4, lwd=2)
# Compare to analytical ('asymptotic') confidence interval 
# (they should be very similar as they both assume normality): 
ci
(round(ci.boot))


#######################################################################
#### 2. Simple linear regression: Estimate change based on an assumed 
#       linear increase, assuming independence / equal variances:

#  Fit linear regression of biomass over time:
fit.lm <- lm(Biomass/1000 ~ Year, data=POP)
summary(fit.lm)
# Plot data again:
with(POP, error.bar(Year, Biomass/1000, (Biomass-1.96*StdErr)/1000, (Biomass+1.96*StdErr)/1000, F,
                    xlab = "", ylab = "Survey biomass (1,000t)", cex.lab=1.8, cex.axis = 1.4,
                    pch=16, cex=1.5, col=4))
# Add fitted regression line
abline(fit.lm, col=2)
# The slope of the fitted line is 26.7 kt / year, or an increase of:
#      29 years * 26.7 kt / year= 775 kt between 1984 and 2013
# This is the same as the difference in predicted values for 1984 and 2013:
fitted(fit.lm)[13] - fitted(fit.lm)[1]
# The linear model may or may not be a reasonable model. 
# What are the assumptions for this model (or any simple linear
# regression model)? How may these assumptions be violated?

# We will deal with quantifying uncertainty in this context later, but we 
# can easily extract the predicted standard errors for each year from 
# the model output:
predict(fit.lm, se=T)
# We can then compute the difference between the predicted 2013 biomass
# and the predicted 1984 biomass and its uncertainty as above, but using
# predicted values and their standard errors instead of the first and last
# year biomasses and their standard errors.


#######################################################################
# The following approaches deal with various violations of the basic
# regression assumptions. This is just a quick preview and we will 
# cover these and other approaches in more detail later.

#### 3. Linear model using a weighted regression approach to account 
#       for unequal variances in the annual estimates of biomass:
fit.w <- lm(Biomass/1000 ~ Year, POP, weights = 1/StdErr^2)
summary(fit.w)
abline(fit.w, col=3)
coef(fit.w)[2] * 29
# Slope: 22.8 kt increase per year or 29*22.8 = 661 kt increase 1984-2013
# This will account for unequal variances, but not lack of independence.


#######################################################################
#### 4. Linear model using a generalized least squares approach to 
#       account for temporal autocorrelation:
library(nlme)
fit.ar <- gls(Biomass/1000 ~ Year, POP, correl = corCAR1(form = ~ Year))
summary(fit.ar)
abline(fit.ar, col=4)
coef(fit.ar)[2] * 29
# Slope: 28.9 kt increase per year or 29*28.9 = 840 kt increase 1984-2013
# This will account for lack of independence (autocorrelation), but not for
#  unequal variances. 


#######################################################################
#### 5. Linear model using a GLS approach to account for both unequal
#       variances and non-independence (temporal autocorrelation)
fit.gls.ar <- gls(Biomass/1000 ~ Year, POP, weights = varFixed(~StdErr^2),
                  correl = corCAR1(form = ~ Year))
summary(fit.gls.ar)
abline(fit.gls.ar, col=1, lty=2)
coef(fit.gls.ar)[2] * 29
# Slope: 22.8 kt increase per year or 29*22.8 = 661 kt increase 1984-2013


#######################################################################
#### 6. Non-linear increase (perhaps more accurate representation 
#       of population dynamics but not based on any specific mechanism):

library(mgcv)
fit.gam <- gam(Biomass/1000 ~ s(Year), data=POP)
summary(fit.gam)
x <- 1984:2013
p <- predict(fit.gam, newdata= data.frame(Year = x))
lines(x,p, col=5, lwd=2)
# Estimated increase between 1984 and 2013:
p[30] - p[1]
# Similar to the linear model, we can extract predicted values 
# for each year and their standard errors from the model object:
predict(fit.gam, se=T)


##################################################################
# As an alternative, using the ggplot2 package, you can 
# more easily add some of the standard fits, including
# confidence bands if desired:
library(ggplot2)
ggplot(POP, aes(x=Year, y=Biomass)) +
  geom_errorbar(aes(ymin = Biomass-1.96*StdErr, ymax=Biomass+1.96*StdErr)) +
  geom_point(color="blue", size=3) +
  geom_smooth(method="lm", se=T, fill="blue", alpha=0.2) +
  geom_smooth(method="loess", se=T, color="red", fill="red", alpha=0.2) +
  xlab("") + ylab(" Biomass (t)") 

