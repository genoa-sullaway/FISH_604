######################################

# Lab 8:    Generalized Linear Model
#           Multinomial regression

# Franz Mueter
# Last modified: 10-07-2021

######################################
library(here)

#### For this lab, You need to install the 'VGAM' package!
if(!("VGAM" %in% installed.packages())) {
  install.packages("VGAM")
}

# Example: Multinomial logit model to model behaviors of 
# sea lion pups before and after branding
# Data from Carole Di Poi / Shannon Atkinson

library(lattice)

##### Import and examine data
brand <- read.csv("data/branding.csv")
brand
# Examine the data: Pups were observed for 14 days total, three
# times daily. Branding took place at the end of day 7.
# Each observations recorded the behavior of all pups that were
# present at the time (note that total N varies as some pups may
# have been in the water or undetected during the observation)
# Behaviors were classified into one of 5 types! The most common
# behavior was 'resting'

##### Prepare data for analysis
# First, extract the response variables (number of animals 
# displaying behavior j):
y <- as.matrix(brand[,4:8])
y

# Re-order treatment factor and 'time of day' variable for display purposes
# (chose order, instead of default alphabetical order)
brand$treatment <- ordered(brand$treatment, levels = c("pre-branding", "post-branding"))
brand$time <- ordered(brand$time, levels = c("8:00", "12:00", "18:00"))
brand$treatment  # Note this is now an ordered vector!
brand$time       # Note this is now an ordered vector!

# Compute observed proportions over time for the two treatment groups 
# and for each behavior:
# (First, compute row totals ('rowSums'), then divide each row by the 
# total to get the observed proportions:
obs <- sweep(y, 1, rowSums(y), "/") 
obs  # Proportions of animals exhibiting each behavior
rowSums(obs)  # Check! Row sums should be = 1 

# Stack data for plotting
# (equivalent to 'pivot_longer' in tidyr)
obs <- stack(as.data.frame(obs))
obs
names(obs) <- c("proportion", "activity")
head(obs)
# Add 'day', 'time', 'treatment' column. These are
# repeated as often as necessary to match rows in 'obs'
obs <- cbind(obs, brand[,1:3])
obs

##### Graphical exploration of data
xyplot(proportion ~ day | treatment, groups=activity, data=obs, cex=1.3, 
       scale=list(x="free"), auto.key=T)
# There are three daily observations for morning, noon, evening
# These are not easy to interpret - no apparent trends!

# Let's separate observations by time of day and pre- /post-branding:
xyplot(proportion ~ day | treatment + time, groups=activity, data=obs, cex=1.3, 
       type="b", lty=2, scale=list(x="free"), auto.key=T)


##### Fitting model to branding data
# Goal: model proportions of animals displaying each behavior over time to test 
#   for differences before and after branding. 
#   We need to account for possible differences in behavior by time of day and we 
#   also allow for the possibility that the proportion of animals exhibiting a 
#   given behavior changes (linearly) over time as he animals get older!


# Attach package for vectorized models to fit multinomial logit:
library(VGAM)

# First, for illustration only, we estimate overall average proportions:
#  that is, we fit an 'Intercept only' model to estimate the "mean response"
# or the average proportions of each behavior
fit <- vglm(y ~ 1, multinomial)   
summary(fit)  # Examine output
# The output requires some explanation. The estimated coefficients are the
# estimated log-odds (on the logit scale), relative to the last behavior 
# (column 5, which by default serves as the "reference behavior")
# thus the coefficients are:  log(p_j / p_ref)
# where p_j is the estimated proportion of animals displaying behavior j
# and p_ref is the proportion of animals displaying the reference behavior
# The reference behavior is used because proportions always add to 1, hence
# we need to estimate only 4 of the 5 proportions to uniquely determine
# the fifth (reference) proportion.

# We can extract the coefficients, i.e. the log-odds ratios:
(cf <- coef(fit))
# Exponentiating the log-odds ratios yields the 'odds ratios' (or the 
# probability of behaviour j relative to the reference behavior):
exp(cf)   
# This implies, for example, that the pups are almost 12 times as likely
# to be resting (intercept 3) than nursing (reference behavior)

# In this simple case, 'cf' directly corresponds to the 'linear predictor' (Xb),
# and the probability of behavior j for j={1,2,3,4} can be computed as shown
# in slide 7 (Mod 6 Generalized Linear Models part 2):
(p <- exp(cf) / (1 + sum(exp(cf)))) 
# This yields the probabilities of the four 'other' behaviors and we can
# simply compute the probability of the reference behavior as:
1-sum(p)
# We'll append the probability of the reference behavior to our vector
p <- c(p, 1-sum(p))    
names(p) <- dimnames(y)[[2]]
round(p*100,2)   # Average proportion of animals (%) exhibiting a given behavior

# Instead of computing these 'by hand', we can - as usual - simply extract
# the predicted values for each observation using the 'predict' function: 

# For example, to compute predicted values on the logit scale:
predict(fit)        
# Note that predicted values are the same for each observation because 
# we fit a simple "intercept" only model, which estimates a mean response 
# across observations. Predicted values on the logit scale are not 
# particularly useful and difficult to interpret, hence it is always
# useful to obtain predicted values on the probability scale (= estimated
# proportions or the probability of an animal exhibiting a given behavior)
predict(fit, type="response")  
fitted(fit)  # same
p   # compare to proportions computed "by hand" above

# Of course, the average proportions are rather trivial and can simply be 
# computed from the data by summing the total number of animals in each 
# column (i.e. all observations of a given behavior) and dividing by the 
# total number of observations:
colSums(y) / sum(y)

# That is only true for the rather trivial 'intercept only' model

# Let's fit a more interesting model we are interested in that serves as our
# "global" or "full" model and allows for:
#    - a linear trend in the proportion of each behavior over time
#           (for example, we may expect proportion of pups nursing to
#            decrease as they get older)
#    - a difference in the proportion of each behavior between morning, 
#      noon, and afternoon
#    - a change in the proportion of each behavior following branding
fit2 <- vglm(y ~ day + time + treatment, multinomial, data=brand)
summary(fit2)   
# Note estimated coefficients for each behavior (relative to reference)
# The coefficients are somewhat difficult to interpret but I provide 
# a basic interpretation below. 

# # Intercepts: 
#     The intercepts in this model correspond to the estimated log-odds ratios 
#     of the first 4 behaviors relative to the reference behavior on day 0
#     for the morning observations (8:00, first level of factor 'time') and
#     "pre-branding" (first level of factor 'treatment'). 
#     For example, Intercept 3 is substantially larger than 0 (1.438), 
#     suggesting that the proportion of resting animals in the morning on 
#     day 0 was significantly larger than the proportion nursing (reference level)

# # 'day' effect:
#    Day is a continuous variable and the coefficient corresponds to the slope
#    of the estimated log-odds over time. All four coefficients are positive
#    with relatively large z-values, suggesting that the proportions of the 
#    first four behaviors RELATIVE TO nursing all increase over time - this makes
#    sense and is as expected because the proportion of pups nursing is likely
#    to decrease as they get older!

# # 'time' effect:
#    The 'time' coefficients estimate the difference in the log-odds between
#    the 3 times of day. The first coefficient (time.L) is the difference 
#    between 12:00 (noon) and 8:00 (morning). This coefficient is negative for
#    all 4 behaviors, suggesting that the probability of a pup to exhibit any of 
#    these behaviors (relative to nursing) is lower at noon than in the morning.
#    Again, this may simply reflect a higher probability of nursing at noon 
#    (the reference level).
#    In contrast, in the evening the pups appear to be more likely to engage
#    in one of the other behaviors (log-odds ratio > 0, which implies an odds
#    ratio > 1).

# # 'treatment' effect:
#    These coefficients correspond to the log-odds ratios of the 4 behaviors
#    (relative to nursing, the reference behavior) of pups after branding
#    compared to before branding, suggesting that the proportion of all four 
#    behaviors, relative to nursing, is lower after branding or that the  
#    proportion nursing increases!

# We can examine the predicted proportions and will visualize them below
# after checking some alternative model formulations
fitted(fit2)    # predicted proportions

#### As an aside, you can fit the same model using the 'multinom' function
# in the package 'nnet':
library(nnet)
fit3 <- multinom(y ~ day + time + treatment, data=brand)
summary(fit3)
# Note that 'multinom' uses the first column as the reference behavior
# hence the interpretation of the coefficient changes, but the predicted
# proportions are (roughly) the same as you can verify by looking at 
# the differences (which are generally VERY small:
(diff <- fitted(fit2) - fitted(fit3))
range(diff)


#### Model comparison / selection
# First, we will check for interactions (thus actually making our 'full' 
# model even "fuller". We can include all 2-way interactions in the model 
# by 'squaring' the model formula (which simply expands the formula to 
# include all two-way interactions):
fit2i <- vglm(y ~ (day + time + treatment)^2, multinomial, data=brand)

# A short aside on formulas:
# Make a formula:
fo <- formula(y ~ (day + time + treatment)^2)
terms(fo)  # Examine 'terms' in the formula to see how it expands
# Thus, the formula is the same as:
formula(y ~ day + time + treatment + day:time + day:treatment + time:treatment)

# Let's take a look at the output:
summary(fit2i)

# Now the interpretation is even more complex and the simplest approach
# to testing interactions for significance is by fitting models with and
# without one or more interactions and comparing the (nested) models 

# Note that standard errors of all interaction coefficients are large relative 
# to the magnitude of the coefficients, hence approximate confidence intervals
# for most coefficients include 0, suggesting that interactions are not 
# significant. 

# We will use two approaches to compare models more formally: 
#  1. a likelihood-ratio test
#  2. comparison of AIC values

## 1. Likelihood ratio test:
#  The test is based on the ratio LAMBDA between the maximized likelihood of a nested 
#  model to that of the "full" model. The likelihood of the full model is always
#  larger than that of a reduced model (which has fewer parameters), hence the 
#  smaller LAMBDA is (< 1), the more likely we are to reject the reduced model in favor
#  of the full (more complex) model. 

#            LAMBDA = L_nested / L_full

# For the comparison, we use the negative log-likelihood, specifically -2log(LAMBDA),
# instead of LAMBDA because it has a known distribution, namely a  Chi-square 
# distribution with degrees of freedom:     df_full - df_nested
# or the difference in degrees of freedom between the full and nested model

# The likelihood-ratio test statistic is computed as:
# LRT = -2*log(likelihood ratio) = -2*(log(L_nested/L_full))
#     = -2*(log(L_nested) - log(L_full))
#     =   2*log(L_full) - 2*log(L_nested)

# There was no 'anova' method implemented for 'VGAM' models until recently
# but we can perform the likelihood ratio test easily enough:

# We can use the 'logLik' function to extract the log-likelihoods from 
# the fitted models to compute the test statistic LRT, where fit2i is
# the 'full' model and fit2 is the 'reduced' (nested) model:
(LRT <- 2*logLik(fit2i) - 2*logLik(fit2))

# The difference in the number of parameters between the two models is
# not included in the 'logLik' results for this model (unlike for most 
# models), but is straightforward:
length(coef(fit2))
length(coef(fit2i))
# The interaction terms add 20 parameters total (5 pairwise interactions
# with 4 parameters each!), hence we compare the LRT statistic to a 
# Chi-square distribution with 20 degrees of freedom. 

#------------------------------------------------------------------------------
# Exercise: Using 'pchisq' compute the probability that a Chi-square 
# distribution with 20 df is larger than or equal to LRT = 17.74:
1-pchisq(LRT, 20)  
#60%  that the dist is larger than LRT provided, which means we opt for the simpler model and we do not reject the null. 
#------------------------------------------------------------------------------

# Note that as of the most recent version of 'vglm', we can simply compute 
# the same test result using:
anova(fit2, fit2i, type = 1)

## Compare models with and without interactions based on AIC:
AIC(fit2); AIC(fit2i)  
# Both the LRt and the AIC values suggest that the model without interactions 
# is the better approximating model in this case (substantially smaller AIC)
# (Note, however, that this does not rule out the possibility that one
# of the interaction terms may be significant on its own, which we could
# test through adding interactions one at a time)

# Finally, we fit several simplified models without one or more of the 
# covariates (time and/or day) and compare models:
fit2.a <- vglm(y ~ day + treatment, multinomial, data=brand)
fit2.b <- vglm(y ~ time + treatment, multinomial, data=brand)
fit2.c <- vglm(y ~ treatment, multinomial, data=brand)

AIC(fit2); AIC(fit2.a); AIC(fit2.b); AIC(fit2.c)  
# The AIC best model includes all three covariates but the difference in 
# AIC between the full model and the model with time of day effect is small, 
# hence the evidence for a 'day' effect (linear trend) is weak!

# We can examine the summary output for individual 
# coefficients that estimate treatment effects:
summary(fit2)  
# Negative coefficients for the treatment terms imply that, relative to nursing,
# each of the other 4 behaviors decreased after branding 
# (or, more accurately, the log-odds decreased after branding)

# To formally check for a significant overall treatment effect, we refit the 
# best model after dropping the treatment effect and compare it to the full 
# model using a likelihood-ratio test or the AIC:
fit1 <- vglm(y ~ day + time, multinomial, data=brand)
# Likelihood ratio test statistic (- 2 times difference in log-likelihood)
L <- - 2*(logLik(fit1) - logLik(fit2))

# P-value for likelihood ratio test (one-sided test):
# (= probability that test statistic is larger than the observed statistic
#  under the null hypothesis that the simpler model (fit1) is true)
pchisq(L, df=1, lower.tail=FALSE) #if this # is small than often the simpler model is true! 

# We can also compare AIC values:
AIC(fit1); AIC(fit2)

# The P-value is not significant and the AIC value of the model without
# treatment is considerably smaller, hence we cannot reject the null hypothesis
# and accept the simpler model (fit1) as the more appropriate model

# Conclusion: There was no significant effect of branding on the 
# frequency of behaviors for sea lion pups after accounting for 
# diel differences and trends over time


##### Goodness of fit / diagnostics:
# We will examine goodness of fit for the full model even though the
# treatment effect was not quite significant

# Fitted values (combine with observed proportions for plotting):
f <- cbind(obs, predicted = stack(as.data.frame(fitted(fit2)))[,1])

xyplot(predicted ~ day | treatment, groups=activity, data=f, cex=1.3,
       scale=list(x="free"), auto.key=T)
xyplot(predicted ~ day | treatment + time, groups=activity, data=f, 
      type="b", scale=list(x="free"))

# Look at behaviors separately:
# Resting, Observed proportions:
xyplot(proportion ~ day | treatment, groups=time, data=f,  cex=1.2,
       scale=list(x="free"), auto.key=T, subset = activity == "rest")
# Resting, Predicted proportions
xyplot(predicted ~ day | treatment, groups=time, data=f, cex=1.2, type="b", 
       scale=list(x="free"), auto.key=T, subset = activity == "rest")
# The model (with treatment effect) suggests a moderate increase in the
# number of pups resting after branding, but the difference is not 
# significant as shown above!

# Nursing, Observed proportions:
xyplot(proportion ~ day | treatment, groups=time, data=f,  cex=1.2,
       scale=list(x="free"), auto.key=T, subset = activity == "nurse")
# Nursing, Predicted proportions
xyplot(predicted ~ day | treatment, groups=time, data=f, cex=1.2, type="b", 
       scale=list(x="free"), auto.key=T, subset = activity == "nurse")
# Note that the frequency of nursing decreases strongly over time, but seems
# to increase immediately after branding (not significiant, see above)

# Other behaviors (fitted values only):
#  Probability of 'walking':
xyplot(predicted ~ day | treatment, groups=time, data=f, cex=1.2, type="b", 
       scale=list(x="free"), auto.key=T, subset = activity == "walk")
#  Probability of being 'alert':
xyplot(predicted ~ day | treatment, groups=time, data=f, cex=1.2, type="b", 
       scale=list(x="free"), auto.key=T, subset = activity == "alert")
#  Probability of 'grooming':
xyplot(predicted ~ day | treatment, groups=time, data=f, cex=1.2, type="b", 
       scale=list(x="free"), auto.key=T, subset = activity == "groom")
# Apparent decrease in number of animals walking and being alert after 
# branding and an apparent increase in grooming

# While the treatment effect was not significant, the estimated change
# in behavior after branding is consistent with a stressful event
# (increases in resting, nursing and grooming)

# Residual diagnostics (full model)
par(mfrow=c(2,4))
plotvglm(fit2)   # (Limited) diagnostics for vglm models, see help file

r <-  resid(fit2, type="pearson")  # Extract Pearson residuals
par(mfrow=c(2,2))

# Residuals against fitted values (same as from 'plotvglm'):
for(i in 1:4) {
  plot(fitted(fit2)[,i], r[,i])
  abline(h=0, lty=2)
  title(names(brand)[i+3])
}

# All residuals plotted over time
for(i in 1:4) {
  plot(brand$day, r[,i])
  abline(h=0, lty=2)
  title(names(brand)[i+3])
}
# Note that there are 3 residuals for each day (morning, noon, evening)
# No major concerns for 'walking' and 'alert', but some outliers for 
# resting behavior (days 1-4) that we may want to explore.

# Negative outliers imply that fewer animals were resting (relative to nursing)
# than predicted by the model (observed < predicted)
# Identify outliers (absolute residual larger than 2:
j <- r[,3] < -2
brand[j,]   # not clear if these have anything in common?

# Daily residuals for a given time of day:
for(i in 1:4) {
  j <- brand$time == "12:00"   # Change to 12:00 and "18:00" to see other residuals
  plot(brand$day[j], r[j,i], type="b")
  abline(h=0, lty=2)
  title(names(brand)[i+3])
}

# Some evidence of serial correlation and apparent increase 
# in noon-time alertness over time!? 

# This could be captured through an inteaction term between 'time' and 'day':
fit4 <- vglm(y ~ day + time + treatment + day:time, 
             family = multinomial, data = brand)
AIC(fit2); AIC(fit4)
# No strong evidence for interaction (which requires 8 additional paramters)

#### A brief GAM add-on:
# We can fit the same model using a more flexible modeling approach, namely
# a generalized additive model that can account for any non-linearities in the 
# trend over time. The full model is easily fit using the 'vgam' function
# simply by replacing the linear term for 'day' in the model formula with 
# a smooth function:

fit.gam <- vgam(y ~ s(day) + time + treatment, family = multinomial, data = brand)
summary(fit.gam)
AIC(fit2); AIC(fit.gam)
# Comparison of AIC values does not indicate a significantly better 
# fit of the GAM over the GLM
# A quick look at the fitted values to see the non-linearity:
f <- cbind(obs, predicted = stack(as.data.frame(fitted(fit.gam)))[,1])
xyplot(predicted ~ day | treatment + time, groups=activity, data=f,
         type="b", scale=list(x="free"), auto.key = F)
xyplot(predicted ~ day | time, groups=activity, data=f, type="b", 
       scale=list(x="free"), subset=treatment=="pre-branding")
# There is some indication that there was a drop in resting (green) and an 
# increase in some of the other activities (in particular walking, blue) on 
# days 3 and 4 (pre-branding), but overall little evidence for any
# non-linearities that we need to worry about!


#### Note on "Hauck-Donner effect"
# Multinomial models are best used with data where there are not too many
# zeros. If probabilities of group membership are "too close to zero or one",
# significance tests and model comparisons become problematic. It is not
# always clear what "too close to zero or one" means, but the Hauck-Donner 
# effect is one way to measure whether there are problematic coefficients
# In this case, that does not appear to be the case, so we are safe.
