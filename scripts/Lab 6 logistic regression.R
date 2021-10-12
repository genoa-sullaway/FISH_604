##############################################

# Lab 6:  Generalized Linear Modeling
# Franz Mueter
# Last modified: 9/24/2019
# 
#####   Logistic regression

#   Two examples: Effects of toxin on moths (in-class example)
#                 Maturity of smooth oreo

##############################################


library(ggplot2)

# 1. Logistic regression example from Venables & Ripley (2002)
# This is the example from Module 6 on the effects of a toxin
# on male and female moths:

# This should be default settings, but just in case:
# (these do not affect model fits but do affect the 
# interpretation of model coefficients)
options(contrasts = c("contr.treatment", "contr.poly"))


# Read in data (based on example in V&R:

# Dose of toxin on log-scale:
ldose <- rep(0:5, 2)  
# Doses are quantified on log-scale (base 2), hence actual doses are:
# 2^0, 2^1, 2^2, 2^3, 2^4, 2^5

# Response variable (number of dead moths out of a batch of 20 moths)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))

# Combine in dataframe
moths <- data.frame(ldose, sex, numdead)

# Create a response 'matrix' for analysis (Number of 'successes' and 'failures')
DA <- cbind(numdead, numalive = 20 - numdead)
DA  # Examine structure of response

# Visualize data (proportions dying as a function of dose):
ggplot(moths, aes(ldose, numdead/20, color=sex)) +
  theme_bw() + labs(x="log(Dose)", y = "Proportion") +
  ggtitle("Proportion of moths dying") +
  geom_point(size=2) 


# Fit 'full' model (with interaction allowing for separate slopes):
moths.lr <- glm(DA ~ sex*ldose, family = binomial)

summary(moths.lr, cor = F)         
# z-test suggests that interaction is not significant
anova(moths.lr, test="Chisq")      
# Sequential chi-square test confirms this result


# More meaningful coefficients can be obtained by shifting the data such that the 
# intercept corresonds to a pre-specified dose. For example, we can compare the 
# intercept between males and females at dose 8 (= 3 on log-scale with base 2):
moths.lr2 <- glm(DA ~ sex*I(ldose-3), family = binomial)
summary(moths.lr2)   
# The coefficient for 'sexM' (i.e. difference between males and females) now 
# corresponds to the difference in the response between males and females 
# at a dose 2^3 = 8 of the toxin. 
log2(8)
#----------------------------------------------------------------------------------
### Exercise 1: 
# (a) Write an expression to compute the probabilities that a male moth and a female 
# moth will die from exposure to the toxin at a dose of 8 (=3 on the log2-scale)
# and at a dose of 16 (=4 on the log2-scale) (do NOT use 'predict')

# You will need the coefficients from the above model:
(cf <- coef(moths.lr2))   # Vector of length 4, saved a vector 'cf' for convenience
# and the equations on slide 28 of Module 6 part 1 ('Computing predicted values')
# Note that the independent variable in the model is 'ldose - 3' rather than 'ldose'

#----  insert your solution here -----

#for moth.lr2 you can look at the intercept, so take  because they just added the offset of 3 to the model, because that changes the intercept 
#if we did this for the original model, would use log2(8) because there is no offset
cf

# (b) Is there a signficiant difference in the probability of dying from a dose of 8 
# (ldose=3) between males and females? Back up your conclusion with an appropriate
# test, test statistic and p-value (all can easily be obtained from the output)

#----  insert your answer here -----
#are the males and female intercepts at different locations? 
anova(moths.lr, test = "Chisq")
#for this output we want to look at sex in the output, so yes there is a significant different. 
#second line is asking if the slope is significant kinda meaningless, is it different or not than 0. 
#third line with sex:ldose (the interaction) indicates that there is no significant difference in slope among sexes
#----------------------------------------------------------------------------------

# As an alternative to the t-test for testing whether the interaction coefficient
# is significantly different from zero, we can test whether the full model (with 
# interaction) provides a better fits than a reduced model without an interaction. 

# The model without the interaction term (equal slope by sex) is as follows:
moths.lr3 <- glm(DA ~ sex + ldose, family = binomial)
summary(moths.lr3)


#----------------------------------------------------------------------------------
### Exercise 2:
# Conduct a likelihood ratio test to compare the model with an interaction (Moths.lr) 
# to the model without the interaction term (moth.lr3).

# You will need to compute the likelihood ratio (LR) from the output. You can do so
# by extracting the likelihoods from each model using the 'logLik' function, 
# for example:
lr<-logLik(moths.lr) #this is log liklihood
lr3<-logLik(moths.lr3) 

# 
#neg log likelihood 

# this is a vector of length 1 that you can use in computations, although it prints 
# with some additional information, notably the number of parameters! 

# As an alternative, you can extract the deviance, which is related to the 
# log-likelihood (proportional to -2*logL):
deviance(moths.lr) 
deviance(moths.lr3) 
# While the likelihood here appears to drop some constants, hence the deviance 
# is NOT equal to -2*logL, you can compute the test statistic for the likelihood
# ratio test using either the difference in deviances or the likelihood ratio.
# To do so, use the equations on slide 15 of Module 6, part 1 and test whether 
# the magnitude of LR is larger than expected by chance. 

# Remember that, under the null hypothesis, the simpler model is the "true" model.
# The LR has a chi-square distibution with p-q degrees of freedom, where p and q are
# the number of parameters for the full and reduced model, respectively.

# Therefore, you can compute the probability of getting a LR statistic as large as
# or larger than the observed value using the cumulative density function for the
# chi-square distribution:
pchisq(q=LR, df=p-q)  
# where LR is the likelihood ratio and you will have to figure out p and q

# What can you conclude from the likelihood ratio test?
 
#-------- Insert code and answers here:
#ratiotest:
r.test = 2*(lr[1]-lr3[1]) #these are log liklihood values.  

p = 4
q = 3

#  1-chisq.stat 
    # means the probability larger than or equal to test statistic 
    # testing whether the interaction term is significant
    # our p-value is large () there is no indication that 
# reject Null hypothesis, the reduced model is better

# If null is true, then there is a 18% chance that you get a test statistic of that magnitude, not unusual, don't reject null.
# Null = simple model is best model 

1-  pchisq(q=r.test, df=p-q)  
  pchisq(q=r.test, df=p-q)  
#--------------------------------------------------------------------------------


# We can generally compare nested models in R using the 'anova' function, which 
# works with a wide variety of models and for most types of models has a choice
# of tests. Appropriate choices in this case are "Chisq" ("Analysis of Deviance")
# or "LRT" for "likelihood ratio test" (which are synonymous):
anova(moths.lr, moths.lr3, test = "Chisq")  
anova(moths.lr, moths.lr3, test = "LRT")  
# Note that an F-test is NOT appropriate in this case and results in a warning:
anova(moths.lr, moths.lr3, test = "F")  

# So far, we have only examined models with and without interaction, implying 
# separate slopes or a common slope for males and females, respectively. The 
# model with a common slope allows for differences in the intercept between 
# males and females, but we may also want to consider another obvious model 
# that fits a single logistic function to the data, regardless of sex:
moths.lr4 <- glm(DA ~ ldose, family = binomial)
summary(moths.lr4)
anova(moths.lr4, moths.lr3, test="LRT")
# A likelihood ratio test rejects the null hypothesis that the simpler model
# (without sex-specific intercept) is adequate, hence we may conclude that 
# model moths.lr3 is our overall best model.

# When comparing more than 2 or 3 models, I will typically use other selection
# criteria such as the Akaike Information Criterion (AIC) instead of pairwise 
# tests, which become difficult to interpret with more than a few models:
AIC(moths.lr, moths.lr3, moths.lr4)
# The AIC is computed from the likelihood and weighs the increase in likelihood
# against the additional parameters in a theoretically sound way based on 
# information theory. Smaller AIC values suggest a model that has more support
# in the data, although if differences are small (say less than 2), we may wish
# to choose the more parsimoneous model that has fewer paramters. 

# In this case, the AIC also argues for model moths.lr3, but the fit of the
# full model is very similar.

### We can easily extract and plot the predicted responses on the logit scale:
# I plotted the actual doses on the x-axis, but on a log-transformed scale
plot(2^ldose, log((numdead/20) / (1-numdead/20)), log="x", type="n",
     ylab = "log-odds", xlab = "Dose", cex.axis=1.2, cex.lab=1.5) 
text(2^ldose, log((numdead/20) / (1-numdead/20)), labels = as.character(sex))

# To add the fitted lines we compute predicted values over a range of x values
# Set up data frame for predicting male responses (logit scale) at pre-specified doses:
ld <- 0:5
Pred <- data.frame(ldose = ld, sex = factor(rep("M", length(ld)), levels = levels(sex)))
Pred
lines(2^ld, predict(moths.lr, Pred, type="link"), col = 3)

# Similarly, we set up data frame for predicting female responses (logit scale) 
# at pre-specified doses:
Pred$sex = factor(rep("F", length(ld)), levels = levels(sex))
lines(2^ld, predict(moths.lr, Pred,  type = "link"), lty = 2, col = 2)

title("Predicted log-odds ratios for full model with interactions")

### Finally, we plot predicted probabilities on the 'response scale' as 
# predicted probabilities of dying:
plot(c(1,32), c(0,1), type = "n", xlab = "dose",
     ylab = "prob", log = "x")
text(2^ldose, numdead/20, labels = as.character(sex))

ld <- seq(0, 5, 0.1) # x-values at which to compute predicted values for plotting
# Set up data frame for predicting male responses at pre-specified doses:
Pred <- data.frame(ldose = ld, sex = factor(rep("M", length(ld)), levels = levels(sex)))
Pred
lines(2^ld, predict(moths.lr, Pred, type = "response"), col = 3)

# Set up data frame for predicting female responses at pre-specified doses:
Pred$sex = factor(rep("F", length(ld)), levels = levels(sex))
Pred
lines(2^ld, predict(moths.lr, Pred, type = "response"), col = 2)

# Compare to fitted values from model without interaction (equal slopes on logit-scale):
lines(2^ld, predict(moths.lr3, Pred, type = "response"), col = 2, lwd=2, lty=2)
Pred$sex = factor(rep("M", length(ld)), levels = levels(sex))
lines(2^ld, predict(moths.lr3, Pred, type = "response"), col = 3, lwd=2, lty=2)
# The difference is fairly small and as suggested by the likelihood ratio test,
# the difference could easily arise by chance (p = 0.184), hence we would
# accept the simpler model (moths.lr3) as adequate 


### Evaluate goodness of fit:
par(mfrow=c(2,2))
plot(moths.lr3, which = 1:4)
# Some evidence for pattern in residuals that is not captured by 'dose'
# hence we may be missing an important covariate or the logistic model 
# is not quite appropriate for the data


###############################################################################

### 2. Maturity of smooth oreo (data courtesy of Andre Punt)
library(here)
#First, import data into R and explore:
oreo <- read.csv("data/oreo.csv")
oreo
names(oreo)
dim(oreo)
summary(oreo)

# Graphical data exploration:
plot(as.factor(Mature) ~ Age, data=oreo)     # Note default "stacked" plot when plotting a categorical variable  

# There is an obvious increase in the proportion of mature fish (light grey) with age
table(oreo$Mature, oreo$Sex)
# Overall (across ages), few of the females are mature (13 of 89)

plot(table(oreo$Mature,oreo$Sex), col=2:3)      # Comparing # of immature M & F and # of mature M & F oreos
plot(table(oreo$Sex,oreo$Mature), col=2:3)      # Comparing # of immature and mature females (ditto for males)
plot(Age~as.factor(Sex), data=oreo)           # Age distribution by sex
oreo$Mature
as.numeric(oreo$Mature)
oreo$Mature <- as.factor(oreo$Mature)
oreo$p <- as.numeric(oreo$Mature)-1    # Convert Yes/No variable to 0s and 1s for modeling
oreo

plot(oreo$Age, oreo$p)           # Note that many circles consist of overlapping data points!
sunflowerplot(oreo$Age, oreo$p)  # to visualize overlapping data points
                                 # the numnber of observations is indicated by number of "leaves"
plot(jitter(oreo$Age), oreo$p)   # Another option is to slightly "jitter" x-locations
plot(jitter(oreo$Age), jitter(oreo$p)) # or we can jitter both x and y, which shows all points but
                                       # can be misleading!

# Visualize maturity at age by sex using ggplot:
p1 <- ggplot(oreo, aes(Age, p, color=Sex)) +
  geom_jitter(width=0.5, height=0.05) +
  facet_wrap(~Sex, ncol=1) 
p1

# Fit the full model (with interaction, i.e. allowing for different slope and intercept by sex)
oreo.full <- glm(p ~ Age*Sex, family=binomial(link=logit),data=oreo)
oreo.full
summary(oreo.full)
oreo.red <- update(oreo.full, ~.-Age:Sex)		# Fit reduced model
summary(oreo.red)

# Compare the two models using an Analysis of Deviance:
anova(oreo.full, oreo.red, test="Chi")
AIC(oreo.full, oreo.red)

#----------------------------------------------------------------------------------
### Exercise 3: Fit a model that only includes age (oreo.age) and compare the 
# three models (oreo.full, oreo.red, oreo.age). Based on Analyses of Deviance 
# and/or AIC values, what is the "best" or preferred model for describing 
# variability in maturity with age? 

# ----------- Solution: -----------

oreo.age <- glm(p ~ Age, family=binomial(link=logit),data=oreo)
summary(oreo.age)

AIC(oreo.age,oreo.red,oreo.full)

#Reduced model (with age + sex) is the best model, but it is quite similiar to the full model. 
#We will choose reduced even though it is close because it is the more parsimonious model. 
#----------------------------------------------------------------------------------

# Create a grid of points for computing and plotting predicted values:
fit <- expand.grid(Age=5:53, Sex=c("F","M"))
fit

# Predicted values on logit scale:
fit$logit <- predict(oreo.full, newdata=fit)
fit
ggplot(fit, aes(Age,logit)) +
  geom_line() +
  facet_wrap(~Sex)

# Predicted values on (back-)transformed probability scale:
fit$p <- exp(fit$logit) / (1+exp(fit$logit))
fit

# Alternative for getting predicted probabilities:
predict(oreo.full, fit, type="response")
# You should make sure they are identical to the ones computed above!

# Plot fitted maturity schedule (probability of being mature) by sex
# in a single plot:
ggplot(fit, aes(Age,p,color=Sex)) +
  geom_line(size=2)
# Note that the lines are shifted relative to each other but increase 
# at the same rate (same slope on scale of linear predictor, or logit-scale)

# Finally, we examine the data and fitted model by sex in separate panels:
p1 + geom_line(aes(Age, p), data=fit)

## ggplot can also do the fitting for us and add confidence intervals
# Separate panels:
p1 + geom_smooth(method="glm", formula =  y ~ x, method.args = list(family = "binomial"))
# Single panel:
ggplot(oreo, aes(Age, p, color=Sex)) +
  geom_jitter(width=0.5, height=0.05) +
  geom_smooth(method="glm", formula =  y ~ x, method.args = list(family = "binomial"))

# Or, using 'visreg':
library(visreg)
visreg(oreo.full, "Age", by="Sex", scale="response", overlay=T, gg=T)

## Diagnostic plot
# These have to be taken with a grain of salt for binomial models:
par(mfrow=c(2,2))
plot(oreo.full, which=1:4)



# Data points with a large influence (as identified from diagnsotic plots:
oreo[68,]  # 26-year old male that is NOT mature
oreo[91,]  # 19-year old female that IS mature

#-----------------------------------------------------------------------------
### Exercise 4: 
# The influential male observation may affect the slope of the line for males 
# Re-fit the models without this influential point and re-examine the model 
# fits. Does removing the outlier change the overall conclusion
# about differences between males and females?

# ------------- Solution: -------------

oreo_smaller <- oreo[-68,]

oreo.age.small <- glm(p ~ Age, family=binomial(link=logit),data=oreo_smaller)
oreo.red.small <- glm(p ~ Age+Sex, family=binomial(link=logit),data=oreo_smaller)
oreo.full.small <- glm(p ~ Age*Sex, family=binomial(link=logit),data=oreo_smaller)

AIC(oreo.age.small,oreo.red.small,oreo.full.small)

summary(oreo.red)
summary(oreo.red.small)

# no removing it does not change overall conclusions about importance of sex
ggplot(oreo_smaller, aes(Age, p, color=Sex)) +
  geom_jitter(width=0.5, height=0.05) +
  geom_smooth(method="glm", formula =  y ~ x, method.args = list(family = "binomial"))

#-----------------------------------------------------------------------------




