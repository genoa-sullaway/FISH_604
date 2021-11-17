############################################################### 

# Lab 12 - Dependent data
#          Spatial autocorrelation

# Last modified: November 11, 2021
# Last run with R v. 4.0.3

###############################################################


### A 'simple' spatial example:

# This portion of the lab works through the basics of spatial modeling using
# a simple example to estimate the mean density (catch-per-unit-effort) of 
# walleye pollock on the eastern Bering Sea shelf. Be warned that this is 
# (more or less) a semester-long course on spatial statistics condensed 
# into a single example but it should give you a flavor for how to deal 
# with spatial correlations.

# We have a random sample of pollock abundances from the eastern Bering Sea
# for two recent years (2009, 2010) and wish to estimate mean densities by 
# year and compare densities between years while accounting for spatial
# correlation in the samples

library(maps)
library(mapdata)
library(mapproj)
library(lattice)
require(ggplot2)
require(nlme)

# Import data
EBS <- read.csv("data/EBS shelf.csv")
str(EBS)
EBS$Year <- factor(EBS$Year)

ggplot(EBS, aes(Longitude, Latitude)) + geom_point() + facet_wrap(~Year)

# We first convert latitude & longitude to x/y coordinates that reflect 
# actual distances between each pair of stations. There are other projections
# we could use but the equal-distance Albers projection provides an easy choice:
# (I chose the parameters for the projection based on a 'rule of thumb' in 
# the help file for 'mapproject')
x <- mapproject(EBS$Longitude, EBS$Latitude, "albers", param=c(55.9, 60.8))
# This results in a list with (among other things) components 'x' and 'y'
# that are the coordinates for each data point in the Albers projection
# We'll add them to our data frame for later use (for estimating spatial 
# autocorrelation as a function of actual distance)
EBS$x <- x$x
EBS$y <- x$y

# Here's the spatial distribution of the sampling stations in the new 
# projection (x,y coordinates), which ensure that distances among points
# are proportional to the actual distances:
ggplot(EBS, aes(x, y)) + geom_point() + facet_wrap(~Year)

# We assume a log-normal distribution and estimate mean CPUE by year on 
# the log-transformed scale. For simplicity, we create a new variable:
EBS$logP <- log(EBS$pollock+1)
# If we assume the CPUE values in each year come from a simple random
# sample of independent observations, we could compare means between 
# years using lm with a factor variable for year:
f.lm <- lm(logP ~ Year, data=EBS)
summary(f.lm)
# Note that the intercept corresponds to the estimated mean log-density 
# in 2009 and the 'Year2010' parameter corresponds to the difference in 
# mean log-density between 2009 and 2010

# The question of interest is if there is a significant difference in 
# pollock density between 2009 and 2010! What can you conclude (if our
# independence assumption is correct)?

### Residual diagnostics. 
plot(f.lm, which=1:4)
# The residuals don't look great but at least there are no particularly 
# influential observations (Cook's distances < 0.02)

# Given the spatial nature of the data, a key diagnostic is a spatial 
# plot of the residuals to examine them for possible autocorrelation. 
# Let's do a quick graphical exploration using the 'bubble' function 
# in the spatial package 'sp':
library(sp)

# We examine autocorrelation for one year at a time because the spatial 
# autocorrelation, if present, should be evident within years only 
# (because the fish may have a quite different spatial distributions 
# between years and there is no reason to think that the CPUE residuals
# at two stations that are close in space but were sampled in different 
# years are correlated with each other.

EBS$r <- resid(f.lm)   # Extract residuals
# Create a spatial data object by specifying which variables correspond 
# to the spatial coordinates:
coordinates(EBS) <- ~ x+y
# Note that this creates a more complex object that includes the data 
# frame as a 'slot'. However, you can still use it like a data frame 
# for most purposes.
str(EBS)

# To examine residuals by year, we create a logical vector to use for 
# subsetting the data frame
(sub2010 <- EBS$Year == "2010")
bubble(EBS[sub2010,], zcol="r", key.entries = seq(-6,6,by=2))
# Note the strong spatial patterns of all positive residuals in some
# regions (middle & outer shelf, where pollock occur in high abundances) 
# and all negative residuals in shallow and deeper regions!

# A formal test for spatial autocorrelation is a little less straightforward
# than for a time series because we now operate in two dimensions. The 
# most common tool for examining spatial autocorrelations in the so-called
# "semi-variogram". The idea is to plot the squared difference in log(CPUE)
# between two samples against the distance separating the samples!

# This is called a 'Variogram' (for historical and mathematical reasons, it
# is really the Semi-Variogram', or half the Variogram, that is used). 

# For each pair of observations x1,x2, the corresponding semi-variogram 
# plots (x1-x2)^2/2 against the distance between observations x1 and x2.

# The variogram typically increases with distance to a plateau because 
# near-by observations tend to be similar, resulting in smaller squared
# differences on average, whereas far away observations tend to be very
# different (independent), resulting in larger squared differences.

# Let's compute the variogram for the 2010 data 'by hand' for illustration:
r <- resid(f.lm)[sub2010]   # Extract residuals for 2010

# Compute pairwise distances among locations based on the coordinates 
# 'x' and 'y':
d <- dist(coordinates(EBS)[sub2010,])
# d is a large (88x88) matrix of pairwise distances (for efficiency, only
# the values below the diagonal are returned!) but we can look at the 
# entire matrix if we convert it to a matrix:
View(as.matrix(round(d,3)))  # Look at entire matrix

# For the 'Variogram' function that we use below, we need to convert
# d to a vector (instead of a 'distance' object)
d <- as.vector(d)   
# The vector d contains only one set of pairwise distances (values 
# below the diagonal of the distance matrix). Hence it has length 
# n * (n-1) / 2, where n is the number of observations

# The 'Variogram' function from the 'nlme' package computes the 
# semi-variogram, given a variable (residuals in this case) and the 
# corresponding pairwise distances: 
SemiVar <- Variogram(r, d)
head(SemiVar, 10)  
# Each row of the result corresponds to a pair of observations where 
# the column 'variog' contains the semi-variogram values and the 
# column 'dist' the corresponding distances between samples. 

# The first row is the distance between observations 2 and 3:
SemiVar[1,]

# We can check this by directly computing the semi-variogram value as 
# (x1-x2)^2/2, where the x's are the residuals of pollock CPUE 
# (= differences between the 2010 mean log(CPUE) and log(cpue) 
#    at a given station):
r[1:2]  # Residuals at station 2 & 3
(r[1]-r[2])^2/2  # Compare to first row of 'SemiVar':

# 'SemiVar' is an object of class "Variogram", for which there is a default 
# plotting option that simply plots the values against distances with a 
# scatterplot smoother, which indicates an increase in the semivariogram 
# values with distance, as would be expected:
plot(SemiVar, xlim=c(0,0.1))
# The further two points are apart the less similar to each other they tend
# to be, although there is huge variability!
# (The inverse of this is that the closer two observations are together
# the more similar - i.e. positively correlated - they are). At zero 
# distance they may be perfectly correlated with each other (if every
# time you sample in the same location you get the same value) or they
# may be correlated with a correlation that is less than 1 (the 'nugget'
# effect) if the data are noisy at small distances (e.g. due to measurement
# error). The 'nugget' effect refers to that variability at small scales.

# Because of the large number of pairwise distances, it is typically more
# useful to plot binned values of the semi-variogram (which is the default
# for many spatial packages such as 'geoR')

# Let's create bins for the distances in increments of 0.005 (on the x/y 
# scale) and plot the distribution of the semi-variogram values in each
# as boxplots:
bins <- cut(SemiVar$dist, seq(0,0.1, by=0.005))  
plot(bins, SemiVar$variog)
# This shows more clearly that differences between observations taken within
# a small distance of each other tend to have low variance, i.e. are similar
# to each other!

# The modeling functions such as 'gls' and 'lme' fit a model to these 
# pairwise distances that quantifies how the semi-variogram increases 
# with distance. Some examples of these models are the exponential, 
# Gaussian, or spherical models.

# Let's look at these three models: 

# These are usually illustrated by plotting the implied correlation
# as a function of distance. Each model has a 'range' paramrter that
# describes the spatial extent of the correlation and a 'nugget' 
# parameter that reflects the correlation between observations at 
# zero distance. If r is the 'range' parameter and n is the 'nugget' 
# effect, we can compute assumed correlations by distance d as follows: 

# 1. Exponential correlation structure:
# The correlation between two observations that are a distance d apart 
# is exp(-d/r) when no nugget effect is present and is (1-n)*exp(-r/d) 
# when a nugget effect is assumed:
Exp <- function(d, r, n) {
  (1-n) * exp(-d/r)
}
curve(Exp(x, r=2, n=0), 0, 6, main="Sperical model", 
      ylim=c(0,1)); abline(h=0, lty=2)
curve(Exp(x, r=2, n=0.4), 0, 6, ylim=c(0,1), 
      add=T, col=4); abline(h=0, lty=2)
legend("topright", c("w/out nugget", "w/ nugget"), col=c(1,4), lty=1)

# 2. Gaussian correlation structure 
# The correlation between two observations that are a distance d apart 
# is exp(-(d/r)^2) when no nugget effect is present and is 
# (1-n)*exp(-(r/d)^2)  when a nugget effect is assumed:
gaus <- function(d, r, n) {
  (1-n) * exp(-(d/r)^2)
}
curve(gaus(x, r=2, n=0), 0, 6, main="Gaussian model", 
      ylim=c(0,1)); abline(h=0, lty=2)
curve(gaus(x, r=2, n=0.4), 0, 6, ylim=c(0,1),
     add=T, col=4); abline(h=0, lty=2)
legend("topright", c("w/out nugget", "w/ nugget"), col=c(1,4), lty=1)

# 3. Spherical correlation structure
# For distances d < r, the correlation is:  
#       1-1.5(d/r)+0.5(d/r)^3 
# when no nugget effect is present and:
#     (1-n)*(1-1.5(d/r)+0.5(d/r)^3) 
# when a nugget effect is assumed. If d >= r the correlation is zero. 
spher <- function(d, r, n) {
  ifelse(d<r, (1-n) * (1 - 1.5*(d/r) + 0.5*(d/r)^3), 0)
}
curve(spher(x, r=2, n=0), 0, 3, main="Sperical model", 
      ylim=c(0,1)); abline(h=0, lty=2)
curve(spher(x, r=2, n=0.4), 0, 3, ylim=c(0,1), 
      add=T, col=4); abline(h=0, lty=2)
legend("topright", c("w/out nugget", "w/ nugget"), col=c(1,4), lty=1)

# We can estimate these (and other) variogram models directly as part of 
# model fitting using the 'gls' function, which finds the maximum 
# likelihood estimates of the parameters of a model relating the 
# response to the explanatory variable (in our simple case a model 
# with two separate means for 2009 and 2010) AND the parameters 
# of the variogram model at the same time.

# We include a "nugget" effect in the model in case there is measurement 
# error at very short distances (that is, repeated measurements taken at 
# the same exact point do not have identical values but have 
# some non-zero variance due to small-scale variability).

# To include residual correlation structures in a model, we can use one of 
# the spatial variogram functions, such as 'corExp' for an exponential model,
# 'corGaus' for the Gaussian model and 'corSpher' for the spherical model

# These functions have a 'value' argument to specify starting values for the 
# range and or nugget, and a 'form' argument that specifies the coordinates 
# (here x and y, or we could use lat/long or another sensible coordinate 
# system) that are used to compute pairwise distances among observations.
# If the data are grouped, such as by year in our case, we need to specify 
# the grouping structure so distances are not computed across groups but
# only within a given year (~ x+y | Year). The default for the 'nugget'
# argument is FALSE (no nugget effect), hence we need to set it to TRUE
# to estimate a nugget effect:
f.exp1 <-gls(logP ~ Year, data=EBS, correl=corExp(form= ~ x+y | Year, nugget=T))
summary(f.exp1)
# The output includes the year effect, as well as estimates of the 'range' 
# and nugget parameters. The former corresponds to the distance over which 
# autocorrelation extends (i.e. roughly where the exponential function starts
# to flatten out) and the latter to the variance at zero distance.
# Note that the estimated nugget effect is very small and we may not need it.
# The range parameter is on the scale of the coordinate system, but we can
# 'translate' it into (statutory) miles based on the latitudinal extent of
# the study region 

#-------------------------------------------------------------------------
### Exercise 1: 
# Compute the 'range' (distance in miles) over which walleye pollock CPUE 
# tends to be correlated by: (1) determining the latitudinal range of the 2010 
# data, keeping in mind that one degree latitude corresponds to 60 nautical
# miles, and (2) determining the y-range in the projected (x,y) coordinate 
# system, and (3) calculating the distance in (statutory) miles corresponding
# to the range parameter.
#1nm = 1.15078 stat mile 

ebs_df <- data.frame(EBS) %>% filter(Year == "2010")
# diff in latitudes, convert to nm by * 60
latdiff<-(62.00832-55.00455) * 60
 #420.2262 nm 

# diff in y coord
ydiff<-0.6746428-0.5504710
#  -0.1241718

# proportion to get equivalent in y coord.
x = (ydiff* 60) / latdiff
  # 0.01752335/latdiff*ydiff
# 0.01752335 range param
# 0.01773883 = 60nm? 

x = (0.01752335 *60)/0.01773883

1.15078 *x #change to stat miles

#68.20806 miles?
#yay got it right! 
#-----------------------------------------------------------------------

# To test if the nugget effect improves the model fit, we can simply fit 
# a model without the nugget parameter (n=0) and compare the two models:
f.exp2 <-gls(logP ~ Year, data=EBS, correl=corExp(form= ~ x+y | Year, nugget=F))
summary(f.exp2)
AIC(f.exp1, f.exp2)  
# The result suggests that the model without a nugget effect is adequate!

# We can easily visualize the estimated variogram model as follows:
plot(Variogram(f.exp2), ylim=c(0,1.1))
# Note that the y-axis ranges from 0 to 1 with 0 being no difference in 
# observed values and 1 representing the maximum difference 
# (corresponding to the "sill" parameter in some variogram models) 
# The variogram is essentially 1 minus the result from our 'Exp'
# function for the exponential correlation with disctance.

# Other spatial packages such as 'geoR' are much more flexible and allow 
# estimation of the range, nugget, and "sill" (i.e. the maximum value of 
# the semi-variogram)

#-------------------------------------------------------------------------
## Exercise 2
# Based on these models, which account for spatial autocorrelation in CPUE,
# what can you conclude about differences in pollock CPUE between 2009
# and 2010? Compare to the result from the linear model! (f.lm)

#---------------------------------------------------------------------

# Although the exponential model is probably adequate, we can fit other 
# correlation structures and compare, for example, a Gaussian and a 
# spherical model as follows:
f.Gaus <-gls(logP ~ Year, data=EBS, correl=corGaus(form= ~ x+y | Year, nugget=F))
summary(f.Gaus)
plot(Variogram(f.Gaus), ylim=c(0,1.1))
# This fit does not look quite a good but these are often hard to judge

f.Spher <-gls(logP ~ Year, data=EBS, correl=corSpher(form= ~ x+y | Year, nugget=F))
summary(f.Spher)
plot(Variogram(f.Spher), ylim=c(0,1.1))
# The spherical model seems clearly inadequate in this case as the 
# resulting semi-variogram values should be roughly centered on 1

# This tends to happen quite frequently and it is important to check 
# your fits!! It is not uncommon for these models to fail to converge. 
# In that case, don't give up! You can provide starting values for the 
# parameters of the autocorrelation function.

# For example, if the spherical model fails to converge, we could use the 
# estimate of the range parameter from the exponential model above as a 
# starting value. This can be specified by the 'value' argument of the 
# corSpher() function, which takes a single value for the range parameter 
# if nugget = F, or a vector of two values if nugget = T:
# (see help file for 'corSpher')
f.Spher <-gls(logP ~ Year, data=EBS, correl=corSpher(value = 0.0175, 
                                         form= ~ x+y | Year, nugget=F))
summary(f.Spher)
plot(Variogram(f.Spher), ylim=c(0,1.1))

# This spherical model does not fit particularly well based on a viusal 
# assessment alone, hence we may want to include a nugget for this model 
# (with starting value for range from previous fit):
f.Spher2 <-gls(logP ~ Year, data=EBS, 
               correl=corSpher(value=0.033, form= ~ x+y | Year, nugget=T))
summary(f.Spher2)
plot(Variogram(f.Spher2), ylim=c(0,1.1))   # Much better!

# We can compare these (non-nested) models using AIC, which clearly favors
# either the exponential model without a nugget or the spherical model with
# a nugget effect:
AIC(f.exp2, f.Gaus, f.Spher2,f.Gaus.nug)

# For parsimony we may want to chose the exponential model, but the choice
# has little effect on parameter estimates and their uncertainty. As long  
# as you include a "reasonable" spatial correlation structure, you
# should be OK to use the model for drawing inferences.
f.Gaus.nug <-gls(logP ~ Year, data=EBS, correl=corGaus(form= ~ x+y | Year, nugget=T))
summary(f.Gaus)
plot(Variogram(f.Gaus), ylim=c(0,1.1))

#------------------------------------------------------------------------
### Exercise 3:
# Compare the estimated differences in mean CPUE between 2009 and 2010 
# from these models and the associated standard errors with those from 
# the simple linear model!
summary(f.lm)
summary( f.Gaus )
#------------------------------------------------------------------------


# We can easily include additional covariates in the model that may help 
# explain the observed spatial patterns in density, such as a dome-shaped 
# temperature effect, reflecting an optimum temperature for walleye pollock:

# For model comparisons, you need to fit all models using 'ML':
f.exp3 <-gls(logP ~ Year + poly(Temperature,2), data=EBS, 
             correl=corExp(form= ~ x+y | Year, nugget=F), method="ML")
summary(f.exp3)
# Check the variogram:
plot(Variogram(f.exp3), ylim=c(0,1.1))  # looks reasonable

# re-fit 'Year' only model using ML for model comparison:
(f.exp4 <- update(f.exp2, method="ML"))
summary(f.exp4)

# Compare models via AIC:
AIC(f.exp3, f.exp4)
# Or use a likelihood ratio test:
anova(f.exp3, f.exp4)
# This suggests that the model with a temperature effect fits 
# substantially better than one without! 

# Hence, we can conclude that the local CPUE of pollock appears to be 
# substantially affected by temperature (assuming the same relationship 
# between temperature and CPUE in both 2009 and 2010), but accounting 
# for differences in temperature does not affect our conclusion that there 
# is no significant difference in CPUE between 2009 and 2010 
# (see 't-table' from summmary output). Note also that the estimate of the 
# range, hence the extent of spatial autocorrelation, is smaller after
# including the temperature effect (compare estimates of 'range' 
# from the summary outputs). This is to be expected because including 
# temperature (which is itself spatially autocorrelated) in the model  
# accounts for some of the autocorrelation in the residuals from the
# simpler model and therefore removes it from the new residuals.

# Let's examine the temperature effect a bit more by setting up a data 
# frame that includes a range of temperature values for 2010 (we only 
# need one year because the temperature effect is assumed to be the 
# same for both years):
new.df <- data.frame(Year = factor("2010"), Temperature = seq(-1.7, 6.3, by=0.1))
# Note that you need to make the 'year' variable  a factor, to match the 
# original data frame (EBS), otherwise 'predict' will not work. 

# Compute and save predicted values and plot result:
new.df$pred <- predict(f.exp3, newdata=new.df)
plot(pred ~ Temperature, data=new.df, type="l", ylab="log(pollock CPUE)", col=2)
# The temperature effect is quite substantial (note range of y-values) and we
# can compute a 'pseudo' R2 value (analogues to other linear models) as 1 minus
# the portion of the variability in log(CPUE) that remains unexplained. The 
# latter is estimated by dividing the sum of squared residuals from the model 
# by the sum of squared residuals from the 'null' model. The null model is a 
# model that includes an intercept (overall mean) only, so we can simply remove 
# the mean and use the sum of the squared residuals from the mean for the 
# denominator:
1 - sum(resid(f.exp3)^2) / sum((EBS$logP - mean(EBS$logP))^2)
# Thus, temperature accounts for nearly a quarter of the variability in CPUE,
# which is quite a bit for highly variable CPUE data!

# Note that this does not account for patterns in the data due to spatial 
# autocorelation, which "explains" additional variability in CPUE. Unfortunately,
# the actual spatial pattern that is 'induced' by autocorrelation is not easy to 
# extract from the 'gls' output ('geoR' has functionality to do so) because the
# focus of 'gls' is on the main effects in the model (i.e. Year and Temperature)
# rather than the spatial pattern due to autocorrelation.

# We could compare the magnitude of the temperature effect to the estimated 
# differences in mean CPUE between 2009 and 2010 by superimposing the 
# estimates for 2009 on those for 2010:
new.df <- data.frame(Year = factor("2009"), Temperature = seq(-1.7, 6.3, by=0.1))
new.df$pred <- predict(f.exp3, newdata=new.df)
lines(pred ~ Temperature, data=new.df, type="l", ylab="log(pollock CPUE)", col=4)
# Note that the estimated CPUE at a given temperature is virtually identical 
# between the two years hence the two lines are almost indistinguishable. 
# Of course the Year effect was not significant in the model, so this is to 
# be expected. Note that we also assumed that the shape of the temperature 
# effect is identical between the two years (no interaction)

##################### THE END
