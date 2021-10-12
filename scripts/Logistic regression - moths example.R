
##############################################

# Generalized Linear Modeling
# Franz Mueter
# Last modified: 9/26/2021
# 
#####   Logistic regression

#   Examples: Effects of toxin on moths (in-class example)
#   The first part is identical to the first ~50 lines
#   in Lab 6, where you will explore this logistic 
#   regression example in more detail

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
p1 <- ggplot(moths, aes(ldose, numdead/20, color=sex)) +
  theme_bw() + labs(x="log(Dose)", y = "Proportion") +
  ggtitle("Proportion of moths dying") +
  geom_point(size=2) 
p1

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


### Visualize model fit 

# Add predicted values on response scale (back-transformed to probabilities)
# to the original data frame:
X <- data.frame(ldose = rep(seq(0,5,length=51),2), sex = factor(rep(c("M","F"), each = 51)))
X$pred = predict(moths.lr, newdata=X, type="response")
X  # Examine 'X'

# Plot data and fitted model
p1 + geom_line(data=X, aes(x=ldose, y=pred, color=sex))



###############################################################################
# For illustration, we fit the model with "true" binary data, rather than a 
# two column matrix of the proportions dead/alive

# Expand data set to have one row per moth:
.rows <- rep(c(1:12), moths$numdead)
dead <- moths[.rows, c("ldose", "sex")]
.rows <- rep(c(1:12), 20-moths$numdead)
alive <-  moths[.rows, c("ldose", "sex")]
dead$dead <- 1
alive$dead <- 0

# Data frame with one row per moth:
(moths2 <- rbind(dead, alive))
nrow(moths2)  # Should be 12 * 20

# Fit logistic regression model using 0/1 response:
moths2.lr <- glm(dead ~ sex*ldose, data=moths2, family = binomial)
summary(moths2.lr)
# Compare coefficients (should be identical):
coef(moths2.lr)
coef(moths.lr)

# Visualize model fit with 95% confidence bands using visreg
library(visreg)
visreg(moths2.lr, xvar="ldose", by="sex",gg=T, scale="response", overlay=T, alpha=.05)

# The results are identical to:
visreg(moths.lr, xvar="ldose", by="sex",gg=T, scale="response", overlay=T, alpha=.05)


