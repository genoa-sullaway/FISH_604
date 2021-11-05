############################################################### 

# Lab 11 - Linear mixed effects model of salmon survival rates
# Data from Randall Peterman (SFU)
# Pyper et al. (2001), Mueter et al. (2002)

# Last modified: October 30, 2021
# Last run with R 4.0.1

###############################################################

# We will use pink salmon stock-recruitment data with an
# environmental covariate to fit several different generalized, 
# multi-stock, mixed-effects Ricker models of recruitment

library(nlme)
library(tidyverse)
library(here)

### Import data
pinks <- read.csv("data/pinks.csv")
head(pinks)
str(pinks)
# Look at the structure of the 'pinks' data frame. It has data 
# on spawners (S), recruits (R), log(R/S) (SR), and local SST 
# anomalies for each of 29 pink salmon stocks from Alaska for 
# different numbers of brood years (17-37 years):
# (Numbers are in thousands)
table(pinks$stock)

# Convert 'stock' to ordered factor that maintains south-to-north order:
unique(pinks$stock)
pinks$stock <- factor(pinks$stock, levels=unique(pinks$stock))
levels(pinks$stock) # Check levels

# Let's first examine spawner-recruit relationships by stock:
p1 <- ggplot(pinks, aes(S,R)) + geom_point() 
p1 + facet_wrap(~stock)

# The absolute abundances are on very different scales, hence 
# we should use different x- and y-axes by panel:
p1 + facet_wrap(~stock, scales="free")

# To fit a Ricker stock recruitment model, log-survival (=log(R/S)) 
# is modeled as a linear function of spawner abundance: 
#                                 log(R/S) = a + b*S

p2 <- ggplot(pinks, aes(S,log(R/S))) + geom_point() + facet_wrap(~stock, scales="free")
p2
# Add linear fits by stock:
p2 + geom_smooth(method="lm")

# Note how noisy the data are, which is typical of spawner-recruit data! 

# While the function above fits simple linear regressions and adds the
# fitted lines, it does not save any output.

### Fit the generalized Ricker model that includes effects of spawner 
# abundance (S) and SST on log-survival. The model can be fit as a
# multiple linear regression

# We start by fitting a generalized Ricker model for each stock separately
# to compare the slopes and intercepts across stocks:
by.stock <- lmList(SR ~ S + SST | stock, data = pinks) # this means you are fitting the model to each stock independently

# Examine the parameter estimates (coefficients from multiple linear 
# regression) and their standard errors, t-statistic and p-values for
# testing the null hypothesis that the coefficients are equal to zero:
summary(by.stock) #first table is intercept, second table is slope
coef(by.stock) # use to calculate predicted values, if you do it by hand it is estimate + S + SST * whatever temperature  you are intersted in

intervals(by.stock, level = 0.95)		# Estimates with 95% CIs
plot(intervals(by.stock))

# Note that the density-dependent parameters (beta) are highly stock-specific
# because of the large differences in the absolute number of spawners, 
# whereas the intercept (productivity parameter a) and the SST parameter 
# are much more consistent across stocks, if variable. For our analysis, we 
# assume that both productivity and the SST effect can be modeled as consisting
# of an average productivity and an average SST effect across all stocks with 
# random stock-specific deviations around these average effects.

#---------------------------------------------------------------------------
# Exercise 1: Estimate the median number of recruits returning to Nome
#             at the long-term mean spawner abundance (299K) when
#             sea-surface temperatures (SST) are 1 degree BELOW average
#             and when they are 1 degree ABOVE average. That is: compute  
#             the predicted log-survival (log(R/S)) based on the stock-
#             specific model coefficients, back-transform to get median
#             survival (i.e. no need to apply bias correction), then 
#             solve for R. How much does average recruitment increase
#             when temperatures increases by 2 degrees?
# (Use the coefficients to calculate predicted values or use 'predict')

#by.stock <- lmList(SR ~ S + SST | stock, data = pinks) # this means you are fitting the model to each stock independently
 
# ----------- Your solution & answer:
#extract coefficients for Nome
coeff_df<-coef(by.stock) %>%
  rownames_to_column() %>%
  filter(rowname == "Nome")

#get mean sst for predicting
# sst<-pinks %>%
#   filter(stock == "Nome") %>%
#   summarise(avg_SST = mean(SST))  %>%
#   mutate(above_SST = 1+avg_SST, below_SST = 1-avg_SST )  

spawn <- 299
#log(R/S) = a = B*S + SST*y

predicted_below <- coeff_df$`(Intercept)` - spawn*coeff_df$S + coeff_df$SST*-1
exp(predicted_below)

predicted_above <- coeff_df$`(Intercept)` - spawn*coeff_df$S + coeff_df$SST*1
exp(predicted_above)

log(R/S)

#---------------------------------------------------------------------------

# When modeling the stock-specific effects as random effects, we assume 
# that they are approximately normally distributed. We can (and should) 
# check this assumption by examining the distribution of the coefficients 
# (productivity and SST effect) across stocks from the separate models:

a <- coef(by.stock)$"(Intercept)"  	
hist(a)
qqnorm(a); qqline(a)    # Approximately normally distributed

sst.cf <- coef(by.stock)$SST
hist(sst.cf)
qqnorm(sst.cf); qqline(sst.cf)   # Approximately normally distributed

# These are approximately normally distributed, hence it is reasonable 
# to fit a mixed-effects model to all stocks simultaneously:

# Multi-stock, mixed-effects model (see lecture):    

#           log(R.i/S.i) = a + a.i - beta.i*S.i + gamma*SST + g.i*SST

# where subscript i is for the ith stock, a, beta.i, and gamma are fixed 
# effects, and the a.i and g.i are random (stock-specific) effects.

### Fit three different models to first examine alternative structures  
#    for the random effects, then test for overall SST effect. 
#    Models are fit via maximum likelihood estimation (for model comparison)

# 1. Random intercept (productivity) and slope (SST effect):
#    Fit via 'ML' to allow for model comparisons
fit1 <- lme(SR ~ SST + stock:S, data = pinks, random =  ~ SST | stock, 
	method = "ML")
(sm <- summary(fit1))  # Somewhat verbose (as it includes var-cov matrix)

# Extract table of (fixed) coefficients:
coef(sm)  

# Random effects by stock:
ranef(sm)

# Extract variances/standard deviation of random effects & residuals
VarCorr(sm)

# Diagnostic plots:
plot(fit1)   # residual against fitted values
plot(fit1, resid(.) ~ fitted(.) | stock, scales=list(x = list(relation="free")))		# Residuals by stock
plot(fit1, resid(.) ~ SST | stock)		      # Residuals against SST

# The output includes the fixed effects coefficients (with their 
# standard errors returned by summary()) and the standard deviation 
# of the random effects (which measure the variability of the 
# stock-specific productivity and SST effect around the average 
# effect). The diagnostic plot is useful to show residuals by group 
# (against fitted values or against any of the covariates).


#---------------------------------------------------------------------------
# Exercise 2: Questions:
# 1. What is the average productivity across stocks and its uncertainty?'

# 2. What is the variability of the productivity parameter across stocks?

# 3. How does SST affect the response {log(R/S)} and how much does the 
#    estimated effect vary across stocks? 

#---------------------------------------------------------------------------

# Confidence intervals for all random effects parameters:
intervals(fit1)

ranef(fit1)   # (predicted) random effects
plot(ranef(fit1))
# Note that these are predicted values and NOT estimated coefficients 
# of the model. The only coefficients associated with these predicted 
# values are the variances (or standard deviations) of the random 
# productivity and SST effects and their covariance


# 2. Fit a random intercepts model only with overall fixed SST effect 
# (same SST effect across stocks, stock-specific productivity parameter):

fit2 <- lme(SR ~ SST + stock:S, data = pinks, random =  ~ 1 | stock, 
	method = "ML")

# Random intercept model without SST:
fit3 <- lme(SR ~ stock:S, data = pinks, random =  ~ 1 | stock, 
	method = "ML")

#---------------------------------------------------------------------------
# Exercise 3: 
# Compare the three models using the AIC or AICc criterion. What can you 
# conclude about the effects of SST on pink salmon survival?

# ----------- Your solution & answer:

# -----------

#---------------------------------------------------------------------------

# Note that we fit the models using maximum likelihood estimation (ML), 
# rather than restricted maximum likelihood estimation (REML). 
# We cannot compare models that were fit via REML unless both models
# have the same fixed effects structure because the likelihood is not
# comparable when fixed effects differ

# Let's update our model fits to use 'REML' and compare the results
# with those from 'ML' fits. The 'nlme' package provides some nice 
# functionality to compare random effects from different models:

fit2.REML <- update(fit2, method = "REML")
plot(compareFits(ranef(fit2), ranef(fit2.REML)))
# estimated productivity parameters based on 'ML' 
# and 'REML' are very close or identical!

fit1.REML <- update(fit1, method = "REML")
plot(compareFits(ranef(fit1), ranef(fit1.REML))) 
pairs(compareFits(ranef(fit1), ranef(fit1.REML)))
# The pairs function shows a scatterplot of the random intercepts and 
# slopes with lines connecting the coefficients from the tow models 
# for the same stocks.

# Again, the productivity parameter are almost identical between the ML 
# and REML fits, but the SST effects change slightly. The ML estimates
# are known to be biased, hence REML estimates are preferred.

# More importantly, we can compare the estimated SST effects based on 
# the stock-specific fixed effects fits with the predicted SST effects 
# from the mixed-effects model for each stock, which are the sum of the 
# overall fixed effect + the random SST effect for each stock:
cf.fix <- coef(by.stock)[,c("(Intercept)", "SST")]
cf.rand <- coef(fit1.REML)[,c("(Intercept)", "SST")]
plot(compareFits(cf.fix, cf.rand))
# Note how the SST effects "shrink" towards the mean, which is particularly
# noticeable for stocks with fairly extreme SST coefficients. For example,
# the SST coefficient for Humpy Creek (near Yakutat) was negative when the 
# generalized Ricker model is fit to Humpy Creek data alone. However, the
# Humpy Creek data are noisy and the time series relatively short. Therefore,
# the model attributes much of the noise in the data to residual variability
# instead of attributing it to a negative SST effect. In that sense, it 
# 'borrows' information from other stocks regarding the SST effect to 
# estimate an SST effect that is not "too different" from that of the other
# stocks, and attributes the remaining year-to-year variability in Humpy 
# Creek to 'other' unexplained (residual) variability. 

# We can also visualize the SST effects and their variability across stocks for 
# both the stock-specific fits and for the mixed-effects model by selecting
# and fixing spawner abundance for each stock and predicting log-survival
# at different temperatures:

# First, set up a data frame for predictions that includes the stock name
# and, for each stock, the mean spawner abundance mean(S) and a range of
# SST anomalies from -2 to 2:
sst.pred <- pinks %>% group_by(stock) %>% summarize(S=mean(S))
sst.pred <- expand_grid(sst.pred, SST=seq(-2, 2, by=0.1))
sst.pred
# Compute predicted values for each row in 'sst.pred':
sst.pred$by.stock <- predict(by.stock, newdata=sst.pred)
sst.pred$mixef <- predict(fit1.REML, newdata=sst.pred)

# Stock-specific fixed estimates of SST effects:
sst.pred %>% ggplot() + 
  theme_bw() +
  geom_path(aes(SST, by.stock, color=stock)) +
  ylab("Effect of log-survival")

sst.pred %>% ggplot() + 
  theme_bw() +
  geom_path(aes(SST, mixef, color=stock)) +
  ylab("Effect of log-survival")


#---------------------------------------------------------------------------
# Exercise 4: 
# As in exercise 1, estimate the median number of recruits returning to Nome
# at the long-term mean spawner abundance (299K) when sea-surface temperatures
# (SST) are 1 degree BELOW average and when they are 1 degree ABOVE average.

# However, this time, base your predicted values on the mixed-effects model
# with random intercept and slope (fit1.REML).
# That is: compute the predicted log-survival (log(R/S)) based on the sum
# of the overall fixed effects and of the stock-specific, random effects. 

# Because the 'predict' function does not work correctly with this model, 
# uou will need to directly calculate the predicted values Using the 
# coefficients from the model, which are extracted using:
fixef(fit1.REML)
ranef(fit1.REML)

# Note that you will only need the overall intercept, overall SST effect,
# and the stock-specific effect of spawner abundance for Nome:
(cf.fix <- fixef(fit1.REML)[c("(Intercept)", "SST", "stockNome:S")])
# and the corresponding random effects:
(cf.rand <- ranef(fit1.REML)["Nome",])
# (Also note that this is a data frame with 1 row, rather than a vector, which 
# matters if you want to extract, say, the SST effect for calculations):
cf.rand[2]   # Data frame with 1 row and 1 column
cf.rand[,2]  # Scalar

# The predicted log-survival is given by (see slide 35 in Mod 9):
#      alpha + a_i + beta_i * S + gamma * SST + g_i * SST
# where alpha is the overall intercept, a_i is the stock-specific random 
# intercept, beta_i is the fixed effect for stock i, gamma is the overall
# SST effect and g_i is the stock-specific SST effect.

# ----------- Your solution:

# -----------

# Compare the predicted SST effect to the stock-specific model predictions 
# that you computed in Exercise 1

# ----------- Your answer:

# -----------

#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
# Exercise 5: Load the 'lme4' package and re-fit the model using the
# 'lmer' function in this package. Note that 'lmer' uses a different
# syntax for specifying random effect directly as part of the 
# fixed-effects formula. Look at the help file and examples and fit
# the full model (equivalent to fit1.REML above) using 'lmer'.
# Compare model results between the two algorithms (parameter estimates,
# their uncertainty and / or the predicted values.

# ----------- Your solution & answer:

# -----------

#---------------------------------------------------------------------------
