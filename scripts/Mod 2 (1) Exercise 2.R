############################################################################################

# FISH 604: Modern Appied Statistics for Fisheries
# In class exercise: 2a. Estimate variance using delta method
#                    2b. Simple bootstrap

# (Exercises with solutions)

# Franz Mueter

############################################################################################

#### 2a. Compute the variance of the inverse of a variable

# Let X be the number of salmon returning to a particular stream to spawn
# Assume we know that, on average, approximately X_bar = 500 salmon return 
# to spawn, with a standard deviation of +/- 100 fish.

# We are interested in some "per spawner" ratio
# and wish to know the mean and variance of 1/X

# Use the Delta method to estimate the mean and variance of 1/X:

# First, we can check if linearity is a reasonable approximation near the mean:
plot(200:800, 1/(200:800), type="l", xlab = "x", ylab="y")

mu.x <- 500     # Assumed known
s2.x <- 100^2   # Assumed known
abline(v=mu.x, col=2)  # Assumed mean

# Assumed 95% confidence interval IF x were normally distributed:
abline(v=c(mu.x + 1.96*sqrt(s2.x),mu.x - 1.96*sqrt(s2.x)) , col=2, lty=2)  

# Over the range of likely x-values, the transformed variable (1/x)
# shows some curvature, but not too extreme, and it curves in a 
# smooth manner that a Taylor expansion can reasonably approximate

#----------------------------------------------------------------------
#### Exercise 2a.1:

# (a) Using the formula from Mod 2 Basics(1) for the Delta method, 
# compute the expected mean ("better approximation") and variance 
# of 1/X. We will need the first and second derivatives of 1/X

# How and why does the expected mean from the delta method differ 
# from the value of 1/mu.x?

# Expected mean of Y = 1/X
# (a)
(mu.y <- 1/mu.x + s2.x/mu.x^3)
# Replot everything:
plot(200:800, 1/(200:800), type="l", xlab = "x", ylab="y")
abline(v=mu.x, col=2)  # Assumed mean
abline(v=c(mu.x + 1.96*sqrt(s2.x),mu.x - 1.96*sqrt(s2.x)) , col=2, lty=2)  
abline(h=mu.y, col=4)
# The estimated mean (blue line) is larger than 1/mu.x (intersection of red and black lines) 
# because of the concave upwards curvature in 1/X, which implies that the second derivative
# (and hence the second term above) is positive. This is as it should be!

#----------------------------------------------------------------------


# (b) Estimate the standard deviation of 1/X, given E(X)=mu.x and Var(X)=s2.x

# Variance of Y = 1/X
s2.x / mu.x^4

# Standard deviation of Y = 1/X
(sd.y <- sqrt(s2.x / mu.x^4))
# Replot everything:
plot(200:800, 1/(200:800), type="l", xlab = "x", ylab="y")
abline(v=mu.x, col=2)  # Assumed mean
abline(v=c(mu.x + 1.96*sqrt(s2.x),mu.x - 1.96*sqrt(s2.x)) , col=2, lty=2)  
abline(h=mu.y, col=4)
abline(h=c(mu.y + 2*sd.y, mu.y - 2*sd.y), col=4, lty=2)

# --------------------------------------------------------------------


### Exercise 2a.2:
# Assuming X is normally distributed with mean 500 and SD = 100  
# examine how well the approximation actually works as follows:

# 1. Simulate 10,000 random draws from a Normal(500, 100) and save as X (use 'rnorm()'!)
# 2. Compute Y = 1/X and compute the mean and standard deviation of Y
# 3. Evaluate how well the delta method approximates the "true" mean and variance
#    (assuming that the mean and and SD estimated in (2) approximate the true values)
# 4. Plot a histogram or density estimate of the full distribution of 1/X (simulated values)
# 5. Add horizontal lines at the mean and at +/- 2 standard deviations
#    For a normally distributed random variable, the interval x_bar +/- 2 SD contains roughly 
#    the central 95% of observations. 
#    Does the interval for 1/X provide an appropriate 95% "confidence interval"

X <- rnorm(10000, 500, 100)  # 10,000 simulated draws from a normal random variable
mean(X); sd(X)  # Verify that X has correct mean & SD
Y <- 1/X
mean(Y)  # Estimated mean of Y based on 10,000 simulated values
mu.y     # Expected mean based on delta method (pretty good!!)
sd(Y)    # Estimated standard deviation of Y based on 10,000 simulated values
sd.y     # The estimate of standard deviation from the delta method

# In the simulation, we generate 10,000 random normal values (X) and the inverse 
# 1/X for each simulated value of X. This gives us an approximation of the "true" 
# mean [mean(Y)] and the standard deviation sd(Y) of 1/X. The simulated mean is
# very close to the mean estimated using the delta method (mu.y). The simulated 
# standard deviation is typically somewhat larger than the standard deviation 
# estimated by the delta method: 0.0004. Repeated simulations (or increasing the 
# number of simulated values) suggests that the 'true' standard deviation is 
# approximately 0.0005 or about 25% larger than that estimated by the delta method.
# We could improve the estimate by using more terms of the Taylor expansion that
# include higher-order derivatives

# Examine distribution of 1/X and evaluate confidence interval
hist(Y, freq=F, ylim=c(0,1100))
lines(density(Y), col=4, lwd=2)
abline(v=mean(Y), col=2, lwd=2)
abline(v = mean(Y) + c(2,-2)*sd(Y), col=2, lty=2)

# The interval appears to provide a reaonable 95% confidence interval altough it
# may be preferable in this case to constuct an asymmertical confidence interval
# because the distribution of 1/X is strongly right-skewed and far from normal
# This could be accomplished using a bootstrap, as shown below
# --------------------------------------------------------------------




#############################################################################################
# 2b. A simple bootstrap example to estimate the variance of the median:

# Returning to the spawner example, assume we have repeated aerial estimates of 
# spawner abundances and we wish to use the data to estimate a more robust measure
# of 1/X and it's variability, for example the median and it's variance

# Observations X:
X <- c(455, 583, 643, 596, 424, 545, 498, 412, 461, 513, 472, 616, 645, 655, 593,
       491, 632, 467, 525, 517, 491, 589, 546, 615, 546, 702, 540, 534, 573, 564)

Y <- 1/X
hist(Y)

median(Y)  # Estimate of median from data

# Bootstrap estimates of median and its variance
out <- numeric(999)  # Set up a vector for the output (see for loop below)
for(i in 1:999) {
  Y.boot <- sample(Y, replace=T)
  out[i] <- median(Y.boot) 
}
hist(out)
mean(out)  # Bootstrap estimate of median
sd(out)    # Bootstrap estimate of standard error of median

# IF we assume the median follows a normal distribution, an approximate 
# 95% confidence interval is given by +/- 2 sd:
CI1 <- mean(out) + c(-2,2) * sd(out)
CI1
abline(v=CI1, col=2, lwd=2)

# A (possibly better) confidence interval that avoids the normal assumption 
# could be based on the lower 2.5th and upper 2.5th percentiles of the 
# bootstrap distribution of the median (i.e. that part of the distribution 
# that contains the central 95% of the bootstrapped values):
CI2 <- quantile(out, c(0.025, 0.975))
CI2
abline(v=CI2, col=4, lwd=2)
# In this case, the two (equally valid) confidence intervals agree quite well!

