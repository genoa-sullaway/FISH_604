######################################

# Module 3: Exploratory Data Analysis
#           part 2: Transformations
# Franz Mueter
# Last modified: 10 September 2019

######################################

# Example and exercise

# Import a dataset that shows plant species richness 
# and area for 17 Galapagos Islands
# (Data from Table 8.2 in Gotelli and Ellison 2004)
# We will use the data to explore the species area relationship
#   and the use of data transformations with a real dataset:

richness <- read.csv("Species richness.csv", row.names=1)
richness  # Inspect the data

# Modeling the species-area relationship:
# The species (S) - area (A) relationship is often modeled as:
#                  S = c * A^z
# where c and z are parameters

## Exercise 1:
# Show (analytically AND by plotting the raw and transformed data)
# how log transformations can deal with both linearizing this 
# relationship and dealing with heteroscedasticity

# ---------------------------------------------------------------------

# Solutions:

# Analytical: take logrithm of both sides to transform to linearity!

# log(S) = log(c * A^z)
#        = log(c) + log(A^z)
#        = log(c) + z*log(A)
#        = c' + z*log(A)    # which is a simple linear regression model
                            # with intercept c'=log(c) and slope z 

# Graphical: plot S against A, then plot S against log(A),
#            finally, plot log(S) against log(A)
#            Discuss patterns
# First, save number of species and area as S and A for brevity:
S <- richness$Num_Species
A <- richness$Area

plot(S ~ A)

# log-transform areas for more even spread along x-axis:
plot(S ~ log(A))  
# Note strong heteroscedasticity (variance increases with A)

# We can add a smoother to examine the relationship
scatter.smooth(log(A), S) 
# This suggests that the increase in S with area is not linear

# log-transform both the dependent and independent variable
plot(log(S) ~ log(A))  
# The relationship looks reasonably linear and there is little
# evidence of heteroscedasticity

# --------------------------------------------------------------------

# Theory suggests that we should log-transform S (see above),
# but we may want to check if there is statistical support for a
# particular transformation. The function 'boxcox' is based on 
# theory developed in a paper by Box and Cox (1964), two British
# statisticians who collaborated on this one paper because they 
# thought it would be great to have both of their names on a paper!

# The 'BoxCox' approach allows us to determine the transformation
# that results in residuals that are close to normal. The function
# uses a variety of different power transformations (defined by
# the value of 'lambda') on a variable y:
#          y.tf = (y^lambda - 1) / lambda 
# to find the lambda that makes the residuals from a regression
# involving y.tf as response variable to be as close to normal
# as possible.

# This transformation always maps a value of y=1 to the 
# transformed value y.tf = 0 and has the feature that 
#    y.tf = log(y)   as lambda approaches 0

# Using this transformation, we can assess a range of lambda
# values and compute the 'normal likelihood' based on some linear
# model. The normal likelihood is simply the log of the residual
# sum of squares from a model:
#     - n/2 * log(sum(r^2))
# where r is a vector of residuals r = y - y^ (y^ are the fitted values)

# Finding the lambda that maximizes the likelihood is therefore 
# the same as finding the lambda that minimizes the standard 
# deviation of the residuals (which is this case are standardized
# in a way such that they can be compared across models with
# different lambda values)

# The boxcox function plots the likelihood over a range of 
# lambdas and shows a 95% confidence interval for lambda
# based on maximum likelihood theory (which we'll cover later).

# The boxcox transformation is "conditional" on a (linear) model fitted to
# the data. The exact model doesn't matter much but you want to select a 
# model that fits the data reasonably well. For example in our case, the 
# plot suggest that a quadratic function may be adequate: 

scatter.smooth(S ~ log(A)) # Smooth fit for comparison
fit <- lm(S ~ poly(log(A),2))
# Add fitted line:
j <- order(A)  # to sort x- and y-values in same order for plotting
lines(log(A)[j], fitted(fit)[j], col=2)


library(MASS)
boxcox(fit)
# The result suggests a transformation with lambda close to 0.25 (fourth-root):
#    y.tf = y^0.25 = sqrt(sqrt(x))
abline(v=0.25, col=2)
# Note that while the transformation used is (y^lambda - 1) / lambda
# For convenience, we simply use the power transformation:
#        y.tf = y^lambda
# after we determine the appropriate lambda (subtracting 1 and
# dividing by lambda doesn't change the 'shape' but is needed for
# standardization in the algorithm)

# Although the result here suggests a fourth-root transformation,
# the confidence interval included zero, hence we may prefer the
# log-transformation (lambda=0) that matches the generally
# accepted form of the species-area relationship and has some 
# support in ecological theory!

# Hence lets log-transform both A and S:
plot(log(S) ~ log(A))   # Approximately linear and equal variance at each level of A


# Exercise 2: 
# Finally, use 'lm' to fit an appropriate linear model (a simple linear regression)
# to model log(species richness) as a function of log(island area):
# Use 'lm' and save resulting model fit as 'SA.fit'

SA.fit <- lm(log(S)~log(A))

summary(SA.fit)
par(mfrow=c(2,2))
plot(SA.fit)  # No pattern in residuals, no obvious violation of homoscedasticity or normality
par(mfrow=c(1,1))


# ---------------------------------------------------------------------------
# Plot fitted model:
plot(log(S) ~ log(A))
abline(SA.fit)
# ---------------------------------------------------------------------------

# Or, visualizing the increase in log(S) by area:
library(visreg)
visreg(SA.fit)

# And, finally, back-transforming the y-axis 
# to actual species numbers (richness):
visreg(SA.fit, trans=exp, add=T)
# Note that this does NOT apply a bias correction
# Note the asymmetrical confidence interval
# Note that residuals are not included because
#    of the transformation of the response
#    (residuals are only defined on log-transformed scale)

# However, we can manually add the fit to the data.
# As we (i.e. you) will be doing this sort of thing repeatedly,
# let's go through it step by step:

# 1. plot data on untransformed scale:
plot(A, S, pch=16, cex=1.4, ylim=c(0,600))

# 2. Set up a data frame that contains a range of values
#    for island area A. These will be the values of A
#    for which we compute predicted species richness
new <- data.frame(A=seq(0,6000,50))
#    (Note that this needs to be a data frame with one 
#    column for each of the explanatory variables because
#    that is what the 'predicxt' function requires)

# 3. Compute predicted values and standard errors for
#    the range of island areas A in our new data frame:
#    The 'predict' function computes these based on the
#    fitted model coefficients and the variance-covariance
#    matrix in 'SA.fit'. In this simple case we could also
#    easily compute these 'by hand'
p <- predict(SA.fit, newdata=new, se=T)
# The object 'p' is a list that contains the fitted values:
p$fit
# and the standard errors (among other things):
p$se.fit

# 4. Compute lower and upper values of the confidence
#    interval for each predicted value:
(lwr <- p$fit - 1.96*p$se.fit)
(upr <- p$fit + 1.96*p$se.fit)

# 5. Finally, add the predicted  values and lower and
#    upper confidence bands to the plot. Note that I
#    am simply exponentiating the predicted values
#    without bias correction, hence the fitted line
#    shows the predicted median species richness for 
#    a given island area:
lines(new$A, exp(p$fit), col=4, lwd=2)
lines(new$A, exp(lwr), col=4, lty=2)
lines(new$A, exp(upr), col=4, lty=2)

# Exercise 3: 
# Add an 80% confidence interval instead of the 95% 
# confidence interval that we added above
# Hint: Replace the multipliers +/- 1.96 with the appropriate
#   quantiles of a standard normal distribution, which you 
#   can very easily compute using 'qnorm':
qnorm(c(0.1, 0.9))
# These are the quantiles such that 80% of the 
# probability density of a standard normal distribution 
# falls between these two values

# We can then use these values as multipliers for 
# computing 'lwr' and 'upr':
lwr <- p$fit + qnorm(0.1)*p$se.fit
upr <- p$fit + qnorm(0.9)*p$se.fit
lines(new$A, exp(lwr), col=2, lty=2)
lines(new$A, exp(upr), col=2, lty=2)


