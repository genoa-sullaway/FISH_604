############################################################### 
# Lab 11

# Non-linear mixed effects model of sablefish growth
# Data from Sherri Dressel (ADF&G, Juneau)

# Last modified: October 30, 2021
###############################################################

library(nlme)
library(ggplot2)

# We will use sablefish data (females only) from Chatham Strait 
# to explore the use of mixed-effects models for modeling 
# length-at-age using a 3-parameter non-linear model. 
library(here)
### Import data:
sable <- read.csv("data/sablefish3.csv")
str(sable)
sable$Year <- factor(sable$Year)

# The 'sablefish.csv' file contains ages and lengths of female
# sablefish over a number of years. We are interested in testing 
# for differences in length-at-age (or growth) among years:  

### Exploratory analysis:
# Length-at-age by year
p1 <- ggplot(sable, aes(Age, Length)) + geom_point()
# Across all years:
p1
# By year:
p1 + facet_wrap(~Year)


### Fit von Bertalanffy models across years with random year effect:
# We first need to define a functional relationship that corresponds to the
# model we are trying to fit!

# Ludwig von Bertalanffy growth model:
LvB <- function(t, k, L.inf, t0=0) {
	L.inf*(1-exp(-k*(t-t0)))
}

# Explore some reasonable parameter values for starting values:
# We need reasonable starting values for fitting a non-linear model!
plot(Length ~ Age, data=sable)
# Some trial values for paramters
curve(LvB(x, k=0.05, L.inf=900), 0, 60, ylim=c(0,900), col=2, lwd=2, add=T)
curve(LvB(x, k=0.2, L.inf=900), 0, 60, col=3, lwd=2, add=T)
curve(LvB(x, k=0.07, L.inf=860, t0=-10), 0, 60, col=4, lwd=2, add=T)
# Blue line provides reasonable fit!

### Fit LvB model across all years:
START <- c(k = 0.07, L.inf = 860, t0 = -10)
fit.all <- nls(Length ~ LvB(Age, k, L.inf, t0), data=sable, start=START)
summary(fit.all)
cf <- coef(fit.all)
cf
# Visualize the fitted model:
plot(Length ~ Age, data=sable)
curve(LvB(x, cf[1],cf[2], cf[3]), col=2, lwd=3, add=T)


### As a next step, we fit separate models by year to examine the
#   distribution of parameters across years:
#   We use parameter estimates from the overall fit as starting values
START <- coef(fit.all)		
by.year <- nlsList(Length ~ LvB(Age, k, L.inf, t0) | Year, data=sable, 
	start=START, control = list(maxiter=200)) 	
# You may get a warning or the model may not converge, in spite of increasing
# the number of iterations. This is not uncommon in non-linear models!!

# Exclude first four years which have poor data (particularly 1990!)
by.year <- nlsList(Length ~ LvB(Age, k, L.inf, t0) | Year, data=sable, 
	start=START, subset = as.numeric(Year) > 4) 	

# Examine parameter estimates and confidence intervals by year:
plot(intervals(by.year))   # Note that some parameters are very 
                           # poorly estimated in some years!

# Examine distribution of coefficients:
(cf <- coef(by.year))   # Matrix of coefficients by year
hist(cf$k)
qqnorm(cf$k); qqline(cf$k)

hist(cf$L.inf)
qqnorm(cf$L.inf); qqline(cf$L.inf)

hist(cf$t0)
qqnorm(cf$t0); qqline(cf$t0)

# Thes have some perhaps troublesome patterns (likely due to having successive years of data)
# but, importantly, no large outliers

cor(coef(by.year), use="pair")
# Note that the coefficients are highly correlated across years, which is typical
# and increases the uncertainty about the value of individual coefficients


### Fit mixed-effects model across years with different random effects
#   structures. We let k and / or L.inf be random across years:
#  (I was unable to fit models with all three as random effect)
fit1 <- nlme(Length ~ LvB(Age, k, L.inf, t0), data=sable, 
	fixed = list(k + L.inf + t0 ~ 1), random = k + L.inf ~ 1 | Year, 
	start = START, method = "REML")
summary(fit1)
plot(ranef(fit1))

# Note the different growth patterns in 1999-2003 (smaller k, larger Linf)

# After fitting the random effects model, we typically examine
# the random effects for approximate normality (even though we
# already checked the stock-specific coefficients)
qqnorm(ranef(fit1)$k); qqline(ranef(fit1)$k)
qqnorm(ranef(fit1)$L.inf); qqline(ranef(fit1)$L.inf)
# These look very reasonable!

#############################
# Exercise 1:
# Which of the two parameters L.inf and k has the largest RELATIVE 
# variability across years. Divide the standard deviation of the 
# parameters across years by the mean parameter estimate to assess
# relative variability

summary(fit1)
VarCorr(fit1) 
 
#k
0.007488061/0.0680

#Linf
33.175824851/882.1714


### It looks like k has the largest relative variability across years. 


#############################

# Only L.inf random:
fit2 <- nlme(Length ~ LvB(Age, k, L.inf, t0), data=sable, 
	fixed = k + L.inf + t0 ~ 1, random = L.inf ~ 1 | Year, 
	start = START, method = "REML")# Model with L.inf random:
summary(fit2)

# Only k random:
fit3 <- nlme(Length ~ LvB(Age, k, L.inf, t0), data=sable, 
	fixed = k + L.inf + t0 ~ 1, random = k ~ 1 | Year, 
	start = START, method = "REML")
summary(fit3)

# Only t0 random:
fit4 <- nlme(Length ~ LvB(Age, k, L.inf, t0), data=sable, 
	fixed = k + L.inf + t0 ~ 1, random = t0 ~ 1 | Year, 
	start = START, method = "REML")
summary(fit4)

#############################
# Exercise 2:
# Compare the overall model with no differences among years (fit.all) and 
# the three random effects models. What can you conclude about variability 
# in growth among years?


df<-AIC(fit.all,
    fit2,
    fit3,
    fit4)
 
# We can conclude there is variability in growth among years and 
# that this is important to include in the model.

#############################

# We can use the 'augPred' function to generate an object for visualizing 
# the data along with the fitted model by year:
plot(augPred(fit1, ~ Age), col=4, lwd=2, cex=0.15)

# Diagnostics:
plot(fit1)   # residual against fitted values
plot(fit1, resid(.) ~ fitted(.) | Year)		# Residuals against fitted by year
plot(fit1, resid(.) ~ Age | Year)		      # Residuals against age by year
plot(fit1, Year ~ resid(.))						# Boxplots of residuals by year
plot(fit1, factor(Age) ~ resid(.))           # Boxplots of residuals by age
# Notice that there appears to be some model mis-specification. For example, residuals
# at younger ages are mostly negative before becoming positive on average at ages 7-9
# In other words the LvB model overestimates length at younger ages and underestimates 
# age at some intermediate ages. To capture the true trend, a more complex growth model
# would be required, for example the 4-parameter Gompertz model 

# The bias seems to be consistent across years:
plot(fit1, factor(Age) ~ resid(.) | Year) 



