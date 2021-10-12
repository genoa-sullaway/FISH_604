###########################################################

# Linear Models in R

# Franz Mueter

# R code only, extracted from 'Module 5 Linear Models lab.Rmd' 
# using 'purl()'

###########################################################

library(MASS)			# attach MASS library for data and functions
mammals           # Look at the mammal data
plot(brain ~ body, data=mammals)
plot(log(brain) ~ log(body), data=mammals)
mammals.lm <- lm(log(brain)~log(body), data=mammals)
summary(mammals.lm)
names(mammals.lm)  	
mammals.lm$call		# the code that we can type to create this model
mammals.lm$coef		# parameter estimates (regression coefficients)
mammals.lm$fit		# fitted values
mammals.lm$res		# residuals
mammals.lm$rank		# number of parameters estimated (p=model degrees of freedom)
mammals.lm$df.res 	# residual degrees of freedom (=n-p)
coef(mammals.lm) 		# same as:  mammals.lm$coef
fitted(mammals.lm) 	# same as: mammals.lm$fitted
residuals(mammals.lm) 	# same as: mammals.lm$resid
# selects one of five different kinds of residuals
resid(mammals.lm, type="pearson")		# short alias for 'residuals'
names(summary(mammals.lm))       # Components of the summary object
summary(mammals.lm)$coef         # Table of 'coefficients' with std errors, etc
s <- summary(mammals.lm)
s$coe	  	# parameter estimates with std. errors and t test (Hyp: coef=0)
s$r.sq    # R-square
s$fstat	  # F-statistic ('value') with numerator and denominator d.f.
x <- log(mammals$body)
y <- log(mammals$brain)
fit0 <- lm(y ~ 1)
summary(fit0)
mean(y); sd(y)	
fit.slope <- lm(y ~ x - 1)
plot(x, y, xlim=c(0,max(x)), ylim=c(0,max(y))); abline(fit.slope)
fit1 <- lm(y ~ x)
summary(fit1)		
anova(fit0, fit1) 		# Compare two model fits using F-test
fit2 <- update(fit1, ~ . + I(x^2))	 
summary(fit2)
anova(fit1, fit2)
fit <- lm(y ~ poly(x, 2))
new <- poly(x, 2)
head(new)              # New, rotated variables
cor(poly(x,2))
coef(fit)
coef(fit2)
plot(fitted(fit), fitted(fit2))
abline(0,1) 	   # Add line with intercept 0, slope 1 (1:1 correspondence line)
termplot(fit, se=T, partial.resid=T)
plot(log(mammals$body), log(mammals$brain))
abline(mammals.lm)		# Add regression line

# Extract values for humans from data frame and log-transform, then highlight humans in plot: 
human <- log(mammals["Human", ])
human
points(human[1], human[2], pch=16, cex=2, col = 2)
text(human[1], human[2]+1, "me", cex=2)
pred <- predict(mammals.lm, newdata=data.frame(body=c(60, 500, 2000)))
plot(log(mammals$body), log(mammals$brain))
abline(mammals.lm)		# Add regression line
points(log(c(60, 500,2000)), pred, pch=16, cex=2, col=3)
plot(mammals$body, mammals$brain)
points(c(60,500,2000), exp(x), pch=16, cex=3, col=3)
par(mfrow=c(2,2))
plot(mammals.lm)		# Do you find any unusual patterns? Outliers?
par(mfrow=c(1,1))
plot(mammals.lm$fit, mammals.lm$res)
abline(h=0)
identify(mammals.lm$fit, mammals.lm$res, row.names(mammals))
(j <- row.names(mammals) != "Human")
mammals.lm2 <- update(mammals.lm, subset = j)
plot(log(mammals$body), log(mammals$brain))
abline(mammals.lm)			      # Original fit	
abline(mammals.lm2, col=2)		# Fit after removing humans
summary(mammals.lm)$coef
summary(mammals.lm2)$coef
coef(mammals.lm) - coef(mammals.lm2)
mammals.infl <- influence(mammals.lm)

# Change in coefficients after dropping humans (compare to results above!):
mammals.infl$coeff["Human", ]

# Change in coefficients for dropping any given observation (only first 5 are shown):
head(mammals.infl$coeff)
mammals.loess <- loess(log(brain)~log(body), data=mammals)
mammals.loess
summary(mammals.loess)  # Output is quite different from a linear model!
names(mammals.loess)    # Components of the model fit
fitted(mammals.loess)	# Fitted values
resid(mammals.loess)
# Plot the fitted values on log-transformed scale:
par(mar=c(5,4,1,1))
plot(log(mammals$body), fitted(mammals.loess), col=6)
# Or we can try to connect the fitted values with a line:
par(mar=c(5,4,1,1))
plot(log(mammals$body), fitted(mammals.loess), col=6)
lines(log(mammals$body), fitted(mammals.loess), col=4)
par(mar=c(5,4,1,1))
x <- log(mammals$body)
y <- mammals.loess$fit
plot(log(mammals$body), log(mammals$brain))
lines(x[order(x)], y[order(x)]) 	 # In order of increasing x!

?cabbages
names(cabbages)
cabbages
library(lattice)
dotplot(Date ~ VitC, data=cabbages) # Distribution of Vitamin C by date
boxplot(cabbages$VitC) # Range and distribution of all Vitamin C values

# Separate Vitamin C content by date, then plot their distribution by date:
VitC.by.date <- split(cabbages$VitC, cabbages$Date)
boxplot(VitC.by.date)

# which is actually the same as:
plot(cabbages$Date, cabbages$VitC)

# Perhaps Vitamin C content various with weight ("head weight"):
plot(VitC~HeadWt, data=cabbages)
scatter.smooth(cabbages$HeadWt, cabbages$VitC)

# The "scatter smooth" fits a LOWESS regression using loess and adds the result
# to a scatter plot, hence is the same as:
plot(VitC~HeadWt, data=cabbages)
(lo <- loess.smooth(cabbages$HeadWt, cabbages$VitC)) # Save output (predicted values)
lines(lo) # Add line to plot based on predicted values over range of x

# Separate Vitamin C contents by cultivar, then plot:
boxplot(split(cabbages$VitC, cabbages$Cult))

# Plot Vitamin C as a function of weight for different levels of Cult
coplot(VitC~HeadWt|Cult, data=cabbages)

# Add scatterplot smoother (LOWESS smoother)
coplot(VitC~HeadWt|Cult, data=cabbages, panel=panel.smooth)

# Plot Vitamin C as a function of weight by date, with LOWESS smoother
coplot(VitC~HeadWt|Date, data=cabbages, panel=panel.smooth, rows=1)

# Plot by date and Cult:
coplot(VitC~HeadWt|Date*Cult, data=cabbages) # OR, equivalently:
coplot(VitC~HeadWt|Date+Cult, data=cabbages)

# Exlore potential interactions between the two factors:
# (Check help file: '?interaction.plot'
with(cabbages, interaction.plot(Cult, Date, VitC))
with(cabbages, interaction.plot(Date, Cult, VitC))

cabbages.aov <- aov(VitC~Cult+Date, data=cabbages)
summary(cabbages.aov)
par(mfrow=c(2,2))
plot(cabbages.aov) # Diagnostic plots
par(mfrow=c(1,1))
plot(cabbages$HeadWt, resid(cabbages.aov), col=as.numeric(cabbages$Date), pch=16)
abline(h=0)
cabbages.lm <- lm(VitC ~ HeadWt, data=cabbages)
summary(cabbages.lm)
par(mfrow=c(2,2))
plot(cabbages.lm) # Residual diagnostics
par(mfrow=c(1,1))
plot(cabbages$Cult, resid(cabbages.lm)) # Residuals by cultivar
plot(cabbages$Date, resid(cabbages.lm)) # Residuals by planting date
cabbages.ancova <- lm(VitC~HeadWt+Cult*Date, data=cabbages)
summary(cabbages.ancova)  # Basic linear model summary with model coefficients
# There is no evidence that interaction term is signficant, hence we drop it:
cabbages.ancova <- update(cabbages.ancova, . ~ . - Cult:Date)
summary(cabbages.ancova)  # Basic linear model summary with model coefficients
anova(cabbages.ancova)    # Anova Table of results for main factors
par(mfrow=c(2,2))
plot(cabbages.lm) # Residual diagnostics
par(mfrow=c(1,1))
plot(cabbages$Cult, resid(cabbages.ancova)); abline(h=0)  # Residuals by cultivar
plot(cabbages$Date, resid(cabbages.ancova)); abline(h=0) # Residuals by planting date
plot(cabbages$HeadWt, resid(cabbages.ancova)); abline(h=0) # Residuals by weight
summary(cabbages.lm)$r.sq
summary(cabbages.ancova)$r.sq
summary(cabbages.aov)$r.sq
var(fitted(cabbages.aov)) / var(cabbages$VitC)

# Or, alternatively, we can fit the same ANOVA using lm() and extract the R^2:
cabbages.aov <- lm(VitC~Cult+Date, data=cabbages)
summary(cabbages.aov)$r.sq	# Note different output format compared to aov!
cabbages.0 <- lm(VitC ~ 1, data=cabbages) # Intercept only!
summary(cabbages.0) 		

# Full model including interactions between all explanatory variables:
cabbages.full <- update(cabbages.0, .~.+HeadWt*Cult*Date)
anova(cabbages.full)
add1(cabbages.0, scope = cabbages.ancova)
deviance(cabbages.0)
deviance(update(cabbages.0, .~.+HeadWt))
deviance(update(cabbages.0, .~.+Cult))
deviance(update(cabbages.0, .~.+Date))
AIC(cabbages.0)
AIC(update(cabbages.0, .~.+HeadWt))
AIC(update(cabbages.0, .~.+Cult))
drop1(cabbages.full, scope = cabbages.full)
cabbages.step <- step(cabbages.full)
summary(cabbages.step) # "Best" model based on step-wise regression
# ANOVA table for final model returned by 'step'
anova(cabbages.step)
# Compare null model, "best" stepwise model, and full model in terms of AIC:
AIC(cabbages.0, cabbages.step, cabbages.full)
library(visreg)
par(mfrow=c(2,2), mar=c(4,4,1,1))
visreg(cabbages.step)
New <- data.frame(Cult="c39", Date="d16", HeadWt=4.0)
predict(cabbages.step, New)

