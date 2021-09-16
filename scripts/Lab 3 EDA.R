#####################################################################

# FISH 604: Modern Applied Statistics for Fisheries
# Lab 3:  Data exploration, multicollinearity, hypothesis testing

# Franz Mueter
# Last modified: August 19, 2021

#####################################################################


######################## 1. Exploratory analysis

Beaufort <- read.csv("Beaufort.csv") 

### Examine the time series using basic scatterplots and trendlines:

par(ask=T)   # interactive mode (you'll have to advance plots by clicking)
par(mar=c(2,2,3,1))
for(i in names(Beaufort)[3:8]) {	# cycle through names in data frame
  # (environmental data only - vars 3-8)
  x <- Beaufort$Year		# extract years for plotting
  y <- Beaufort[,i]			# extract ith variable
  # Plot time series for variable i using both points and lines. 
  # Suppress default axis labels and add variable name as title:
  plot(x, y, type = "b", xlab="", ylab="")
  # Fit linear time trend and add line to plot:
  abline(fit.lm <- lm(y ~ x), col=4)
  slope <- signif(summary(fit.lm)$coef["x", "Estimate"],2)
  p <- round(summary(fit.lm)$coef["x", "Pr(>|t|)"],3)
  # Fit a LOWESS (= LOESS) smooth trend and add to plot:
  fit <- loess(y ~ x)	
  lines(x, fitted(fit), col=2)	# Add fitted line
  title(paste(i, "\n", "Linear trend: slope =", slope, "; p = ", p, sep=""))
}

par(ask=F)  # Return to non-interactive mode!

# Or change the layout to plot all graphs on single page 
# (re-run for loop after changing layout)
par(mfrow=c(3,2))  

# For comparison, plot all environmental time series in a single 
# panel using the 'matplot' function (parralel coordinates plot):
par(mfrow=c(1,1))  # change layout back to single figure
matplot(Beaufort[,-(1:2)], type="l")

# Well, that is not very useful! Let's scale (normalize) all variables 
# first so we can compare them on the same scale:
# (The default behavior for the 'scale' function is to subtracts means from 
# each column and then divide each column by its standard deviation)
matplot(scale(Beaufort[,-1]), type="l")
legend("topright",names(Beaufort)[-1], lty=1:7, col=1:7)

# These types of plots is where 'ggplot' really shines, but you have to reshape 
# the dataset to "long format" first. 
(Beaufort2 <- cbind(Beaufort[,1], stack(as.data.frame(scale(Beaufort[,-(1:2)])))))
names(Beaufort2) <- c("Year", "Index", "Variable")
library(ggplot2)
# Similar to matplot, legend added automatically:
ggplot(Beaufort2, aes(Year, Index, color=Variable)) + 
  geom_line() +
  geom_point()

# With 'ggplot' you can add layers one at a time:
# First save a "base plot" with multiple panels by variable
# We add to tbe base plot later:
p1 <- ggplot(Beaufort2, aes(Year, Index)) +
  geom_point() +
  facet_wrap(~Variable)
p1

# Add a simple linear regression using the function 'lm':
p1 + geom_smooth(method ="lm", se=F)

# You can also easily add a confidence band at the desired level:
p1 + geom_smooth(method ="lm", level=0.95)

# Or, we can add a loess smoother instead of a linear fit:
p1 + geom_smooth(method ="loess", level=0.95)

# Finally, add both loess smoother (with CI) and linear fit 
# (no CI, to avoid cluttering):
p1 + geom_smooth(method="loess", se=T, color="red", fill="orange") +
  geom_smooth(method ="lm", se=F, size=1.2)
# Note that in all cases (with the possible exception of 
# discharge in the Sag River) the linear fit falls well within
# the 95% confidence bands of the loess fits. This can be taken
# as lack of evidence for any non-linear time trends (except
# perhaps for discharge)

########################################################################
### 2. Explore pairwise correlations and address multicollinearity:

# You can generate a basic pairwise plot by simply using the full
# matrix of numerical variables with the 'plot' function:
plot(Beaufort[,-(1:2)])

# Some other useful ways of looking at pairwise relationships!
# You will need to install the 'corrgram' package first.
library(corrgram)
corrgram(Beaufort[,3:8])
corrgram(Beaufort[,3:8], lower.panel=panel.pts, upper.panel=panel.conf,
         diag.panel=panel.density)
corrgram(Beaufort[,3:8], lower.panel=panel.shade, upper.panel=panel.ellipse,
         diag.panel=panel.density)

### Quantifying multicollinearity

cor(Beaufort[,-c(1:2)])		# Pairwise correlations
cor.test(Beaufort$summer.airT, Beaufort$summer.SST)

library(Hmisc)
rcorr(Beaufort)   # You will get error because rcorr does not accept data frame

Env <- as.matrix(Beaufort[,3:8])  # Save environmental variables as a matrix 
(COR <- rcorr(Env))
COR$r   # Matrix of correlations
COR$P   # Matrix of p-value


#-------------------------------------------------------------------------------
# Exercise 1: Compute and test significance of Spearman rank correlations and 
#             compare results to the Pearson product-moment correlations. Are 
#             the same pairwise correlations significant at the 95% level?

# ... Your code and answers here

#-------------------------------------------------------------------------------

# Identify Principal Components:
pca <- princomp(x = Env, cor=T)  # Force use of correlation matrix
summary(pca)
plot(pca)   # Scree plot (variance explained by each PCplot(pca)
print(loadings(pca), cutoff=0)    # show all loadings
pca$scores    # Show all principal component scores!

# The Principal Components Analysis identifies a linear combinations of the
# underlying variables that captures most of the variability in the 'data cloud'
# (major axis or major mode of variability)
# Here, the first PC accounts for about 56% of the variability in the data
# and the loadings suggest that PC1 is positively correlated with temperature
# and negatively correlated with ice cover in the spring and fall. The discharge 
# and Wind variables have less influence on PC1. Thus PC1 primarily reflects 
# temperature conditions (positive values are associated with warm temperatures
# and low ice, negative values with cold temperatures and extensive ice)

# We can use the temperature and ice variables to construct a single variable
# (PC 1) that captures overall temperature conditions in the Beaufort Sea:
pca2 <- princomp(x = Env[,1:4], cor=T) # Include temperature & ice variables only
summary(pca2)
plot(pca2)   # Scree plot (variance explained by each PC
print(loadings(pca2), cutoff=0)    # show all loadings

# Add first PC to data frame:
Beaufort$Temp.PC1 <- pca2$scores[,1]  

# Model Cisco abundances as a function of PC1 (temperature), winds, and discharge, 
# using a multiple linear regression:
fit1 <- lm(Cisco ~ Temp.PC1 + Wind + Sag.discharge, data = Beaufort)
summary(fit1)
par(mfrow=c(2,2))
plot(fit1)
shapiro.test(resid(fit1))
# Note the relatively large outlier in row 8 (highest observed recruitment):
Beaufort

#######################################################################
### 3. Box-Cox transformations

# We check residual diagnostics to assess if a transformation of the response
# variable may be warranted.

# Distribution of residuals:
par(mfrow=c(1,1))
hist(r <- resid(fit1))   # Evidence of skewness?

# Test for skewness:
# install.packages("moments")   # uncomment to install if you don't have it yet!
library(moments)
skewness(r) # Normally distributed data has skewness = 0
# Test for significance:
agostino.test(r)  # marginally significant


# Use Box-Cox approach to select transformation, if any:

library(MASS)
# Fit Box-Cox model after dropping the non-significant variable. 
# Note addition of '1' to the response because of zeros in the response,
# which must be positive in the Box-Cox approach:
boxcox(lm(Cisco+1 ~ Wind + Sag.discharge, data = Beaufort))
# We can 'zoom in' on the range of interest:
boxcox(lm(Cisco+1 ~ Wind + Sag.discharge, data = Beaufort), lambda=seq(-0.6,0.4,by=0.05))

# The best transformation is close to lambda = 0 (corresponding to log-transformation)

# Fit model using log-transformation:
fit2 <- lm(log(Cisco+1) ~ Temp.PC1 + Wind + Sag.discharge, data = Beaufort)
summary(fit2)
par(mfrow=c(2,2)); plot(fit2)

# The diagnostics suggest no major violations of the assumptions 
# (constant variance, normality)

#-----------------------------------------------------------------------------
# Exercise 2: 
# Re-fit model after eliminating non-significant variables:

fit3 <-      # ... Your code and answers here

#-----------------------------------------------------------------------------

# A quick way to visualize the model fit:
termplot(fit3, se=T, partial.resid=T)

# A more flexible function to achieve this is provided in the 'visreg' package
library(visreg)
visreg(fit3)
# or the ggplot version:
visreg(fit3, gg=T)

# We can also add a linear regression with a 95% confidence interval
# directly in ggplot:
ggplot(Beaufort, aes(Wind, log(Cisco+1))) + 
  geom_point() +
  geom_smooth(method="lm")

# Note that this is a poor model because it predicts negative recruitment
# at wind speeds less than 0 (westerly winds)

###################################################################################
### 4. Outliers and robust regression: 

#----------------------------------------------------------------------------

# Exercise 3: Plot SST vs. spring ice conditions
#  Fit linear model (save as 'fit.lm') 
#  and add a fitted line to the plot (abline)

fit.lm <-       # ... Your code here


#----------------------------------------------------------------------------

win.graph(8,4)  		# Open a second graphics window
# quartz(8,4)       # for Mac users
par(mfrow=c(1,2))			# change layout (2 side-by-side plots)
plot(fit.lm, which=c(4,5))		# Show influence plots (see help for ‘plot.lm’)

# First plot shows Cook’s distances, which can be computed by
cooks.distance(fit.lm)

dev.off()				# Close graphics device and switch back to previous device

#----------------------------------------------------------------------------
# Exercise 4: Fit robust regression
library(MASS)

# ... Your code here

#----------------------------------------------------------------------------



#----------------------------------------------------------------------------
# Exercise 5: Plot fitted line after excluding outlier for comparison:

# ... Your code here

#----------------------------------------------------------------------------

