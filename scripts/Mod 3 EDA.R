######################################

# Module 3: Exploratory Data Analysis
# Franz Mueter
# Last updated: August 6, 2021

######################################
#From Kevin S: 
#Here is an alternative qq function I used in another class "qqPlot".    >library(car); qqPlot(x, distribution="gamma", shape = coef(fit)["shape"], rate = coef(fit)["rate"])

# We will use a few data objects that I saved previously as 
# an R object. Unless you opened R by clicking on the file
# 'Mod 3 data.RData', attach it now to get access to the 
# data object that we will use today:
attach("Mod 3 data.RData")
# I will use this approach occasionally to share data or a number 
# of (any kind of) R objects for class or lab. 
search()  # Note the attached '*.RData' workspace in your 'search list'
objects(2)  # Check objects in the workspace
# (Should have 3 objects: GAK1, GAK1.SST, and pollock.cpue)

# Whenever you encounter a new data set, you should quickly explore its structure.
# Usually I take a look at the data using 'str' and, if it is a data frame,
# I may examine the first few rows and perhaps some other relevant aspects.
# You should get in the habit of doings this using whatever tools are appropriate
# for a given object, as it can avoid a lot of mistakes!
str(GAK1)
head(GAK1, 20) 
# Note that the sampling day is given as decimal 'year' and that temperature,
# salinity and density / density differences are shown for discrete depths
# from 0 to 250 m. Since this is a data frame, we can get a useful summary
# of the variables in the data frame using the 'summary' function:
summary(GAK1)
# This lists some useful statistics and would also list any missing values

# The second data frame is a reduced vesion of the first one:
str(GAK1.SST)  # Note that depth column is not included because this only 
               # includes surface values (depth == 0)
head(GAK1.SST)  # Note also the added 'dec.year' column (for plotting)

# Finally, the third data set, which is simply a vector with 1264 values
# corresponding to the catch-per-unit effort (kg/ha) of walleye pollock in 
# Bering Sea bottom trawl samples collected over three years
str(pollock.cpue)  # a vector with 1264 numerical values
hist(pollock.cpue) # Take a quick look at the range and distribution
                   # of values. We'll get back to these data later


### Illustrate 'co-plots':
# Plot GAK1 data for 4 different depth ranges and year ranges:
coplot(salinity ~ temperature | depth * year, data = GAK1, number = c(4,4))

# Tidyverse version 
library(tidyverse)
GAK1 %>%  group_by(DEPTH=cut(depth, 4, labels=c("surface", "mid1", "mid2", "deep")), 
                   YEAR=cut(year,3, labels=c("early", "mid", "late")))  %>% 
ggplot() + 
  geom_point(aes(temperature, salinity)) +
  facet_grid(YEAR ~ DEPTH)


### Illustrate model fitting and visualizing the resulting fitted values

library(lattice)
# Plot seasonal temperature patterns by year
xyplot(temperature ~ dec.year | Year, data = GAK1.SST)

library(nlme)
# Fit a sine wave to data for each year separately with year-specific
# means, amplitudes, periods, and offsets:
START <- list(fixed = c(7.5, -5, 1, 0.8))    # Starting values

# This is a complex non-linear, mixed-effects model that may take a little while to run:
fit <- nlme(temperature ~ mu + a * sin(period * 2 * pi * dec.year + k), data = GAK1.SST, 
	fixed = mu + a + period + k ~ 1, random = mu + a + k + period ~ 1 | Year, start = START)
# Do not worry about any of these details. This is simply intended for illustration

# We can use 'nlme' functionality to easily make a data frame that contains both the data
# and fitted values in a format suitable for plotting. 
fits <- augPred(fit, primary = ~ dec.year)
fits 		# data frame of observed and predicted values
plot(fits)   # Plot of observed and predicted values
# In the 'tidyverse', the 'fortify' function serves a similar purpose but is 
# much more versatile. Throughout class, we will compute fitted values from 
# a variety of models and overlay the fitted values on a plot of the observed 
# values, either 'by hand' or using built-in functions where available.
# That is a very important but sometimes neglected part of any data analysis


### Some graphical approaches to exploring data

## 1. To illustrate a histogram, density estimates, and a normal distribution
# fit to the data, we will use a randomly generated data set (using 'rnorm')
data <- rnorm(200)   # Simulate a sample
hist(data, prob=T, col=5, ylim=c(0,0.6))
# We use the density function to approximate the distribution of the data. The 
# function uses a kernel density estimator. See the help file:
?density
# There are a number of different options, but we'll stick with the default 
# (Gaussian) kernel. The output includes two components, x and y, that correspond
# to predicted densities (y) over equally spaced values of the data (x). Therefore,
# we can directly overlay the density on an existing plot using 'lines' (which, by 
# default, uses components 'x' and 'y' in the input object, if available):
lines(density(data), lwd=3, col=2)

# Next, we generate a sequence of x-values spanning the range of our data, compute 
# the corresponding normal densities for a normal distribution with a mean and 
# standard deviation estimated from the data, and then overlay the normal curve:
x <- seq(-3, 3, by=0.1)
mu <- mean(data); sd <- sd(data)
lines(x, dnorm(x, mu, sd), lwd=3)

## 2. Illustrate a Q-Q plot
x <- rnorm(20)
qqnorm(x, pch=16)
qqline(x, col=2)     # goes through first and third quartiles

## Exercise: Construct a q-q plot from 'scratch' 
# You can do this in R or Excel.
# If you work in R, upload the R code you used to generate the plot
# (preferably as simple text). If you do it in Excel, upload your sheet.

# To get started, pick some data (15-20 values) or simulate random data:
x <- rnorm(15, 10, 4)

#  Compute and plot the sample quantiles (of x) against theoretical quantiles
#  from a normal distribution. Steps:
#    1) compute f_i (see definition, slide 20 in 'Mod3 EDA (1).pdf'), i.e. sort
#       the 'x' vector and compute f_i
#    2) compute the probability density under a normal distribution (q_i) 
#       corresponding to the quantiles 'f_i' using 'dnorm' (see help file)
#    3) plot f_i against q_i



#################################################################################
# Let's look at some real data: Examine a typical distribution of pollock 
# catch-per-unit-effort:
x <- sort(pollock.cpue)

# Check for normality:
qqnorm(x)     # Does not look normal at all (not on straight line)
qqline(x, col=2, lwd=2)

# Check whether data are approx log-normal: 
qqnorm(log(x))		# log-normal fits better, but not too good!
qqline(log(x), col=2, lwd=2)

# Evaluate distribution using histogram:
hist(x)       # Obvious right-skew
hist(log(x))  # Better, but obvious left-skew

# A good distribution to consider for right-skewed data is the gamma distribution
# Let's fit a gamma distribution to the pollock data and assess the fit:
library(MASS)
fit <- fitdistr(pollock.cpue, "gamma")
fit  # Note the two estimated parameters of the gamma distribution

# Construct q-q plot:
#  Sample quantiles:
f.i <- ((1:length(x))-0.5)/length(x)
#  gamma quantiles corresponding to f_i
q.i <- qgamma(f.i, shape = coef(fit)["shape"], rate = coef(fit)["rate"])
plot(q.i, x, pch=16, xlab = "Theoretical quantiles", ylab = "Sample quantiles")
# Line through 10th and 90th quantiles (the default for q-q plots in R):
x1 <- qgamma(c(0.1, 0.9), shape = coef(fit)["shape"], rate = coef(fit)["rate"])
x2 <- quantile(x,c(0.1,0.9))


abline(lm(x2~x1))  # quick way to add straight line that goes through these points

# Gamma approximation perhaps "reasonable", but a bit 'long-tailed' 
# (sample values in the upper tail of the distribution are larger than 
# they "should be" if they were from a gamma distribution). 




########################################################################
# Illustrate Leverage and Influence (see 'Mod3 EDA (2).pdf'):

# Read in some simulated data with two outliers for illustration
x <- scan(nlines=6) 
-1.671768246 -2.034975515 -1.767349106 -1.442374181 -1.098283085
-0.792424332 -0.658611128 -0.792424332 -0.543914095  0.048687238
0.564823883  1.386819281 -0.142474482 -0.008661278 -0.199822998
-0.123358310  0.086919582  0.392778335  1.080960529  0.889798808
0.583940055  0.736869432  1.042728185  1.310354593  1.405935453
1.539748658  9.931748186  9.855283498

y <- scan(nlines=5)
-2.26377604 -1.92783486 -1.50027335 -1.53081346 -1.13379206 -0.64515034
-1.04217174 -1.43919314 -1.10325196 -1.31703271 -0.70623056  0.08781223
1.67589782  1.21779621  0.66807428  0.08781223  0.24051277  0.17943256
1.43157696  1.49265718  1.49265718  1.09563578  1.12617589  0.97347535
0.75969460  0.94293524 10.04388727 -0.18704873


# Plot data and highlight outliers
plot(y~x, type="n", cex.lab=1.5, cex.axis=1.3)
points(y~x, subset=-c(27,28), cex=1.5)
points(y~x, subset=27, pch=17, cex=1.8)  
points(y~x, subset=28, pch=16, cex=1.8)

# High leverage point (triangle) and influential point (solid circle)
# are highlighted!

# Show linear fit with high leverage point included, but without
# influential point (solid line), without either outlier (dashed line),
# and with influential point excluded, but without high leverage point
# (dash-dotted black line)

with.lev.pt <- lm(y~x, subset=-28)
abline(with.lev.pt, lty=1, col=4, lwd=2)

without.lev.pt <- lm(y~x, subset=-c(27,28))
abline(without.lev.pt, lty=2, col=4, lwd=2)

with.infl.pt <- lm(y~x, subset = -27)
abline(with.infl.pt,lty=4, lwd=2)


