###################################################################################

# FISH 604 Modern Applied Statistics
# Module 2: Illustrate expected value and variance of a product and sum
# Franz Mueter
# Last modified: July 31, 2021

###################################################################################

# Example: Estimate total number of crab caught in
#          personal use Dungeness crab pot fishery


# Number of pot lifts in three bays in SE Alaska:
N = c(A=86, B=123, c=58)

# Probability of catching one or more crab in a given pot lift by bay
# (estimated previously):
p = c(A=0.3, B=0.5, C=0.8)

# Average number of crab and standard deviation per trap, given crab are caught 
# (estimated from previous sampling):
mu = c(A=6, B=8, c=8)
sd = c(A=1.8, B=2.2, c=3)

# Let Xi = Ni * fi be the number of successful pot lifts in Bay i
# Let Yi be the number of crab caught in a succeessful pot lift in Bay i

# Expected number of of crab caught in a given bay:
#     E(i) = E(Xi * Yi)
#          = E(Xi) * E(Yi)      
#          = E(Ni * fi) * mui
#          = Ni * E(fi) * mui
#          = Ni *  pi   * mui

(C = N * p * mu)

# Total catch:
sum(C)

# Variance of total catch:
# var(Ci) = var(Xi * Yi)
#          = E(Yi)^2 * var(Xi) + E(Xi)^2 * var(Yi) + var(Xi) * var(Yi)
# with: E(Xi) = Ni * pi; E(Yi) = mui; var(Xi) = Ni * pi * (1-pi); and var(Yi) = sdi^2
# var(Ci) = mui^2 * (Ni * pi * (1-pi)) + (Ni*pi)^2 * sdi^2 + Ni * pi * (1-pi) * sdi^2

# Variance by bay:
(VAR = mu^2 * (N * p * (1-p)) + (N*p)^2 * sd^2 + N * p * (1-p) * sd^2)
# Standard errors by bay:
sqrt(VAR)

# Standard errors of the total number of crab
sd.C <- sqrt(sum(VAR))

# 95% Confidence Interval for total catch assuming normal distribution:
cat(round(sum(C) - 1.96*sd.C), "to",  round(sum(C) + 1.96*sd.C), "Dungeness crab")


###########################################################################
# The next section illustrates some basic simulations, which are sometimes
# useful simply to check your calculations and to get a better feel for the
# variability in derived / estimated quantities

# Simulations are of course also critical to many analyses such as
# randomization tests and ecological modeling

## A basic simulation of the average number of crab that might be caught in  
## Bay A, and its variance, under the assumptions outlined above:

# First, simulate the number of pots (out of N[1] = 86) that successfully catch crab 
# given a probability of p[1] = 0.3 that a pot catches any crab. 
# Here we generate 1000 random draws, assuming a 30% (p[1]) chance of catching anything

# We use the binomial distribution: Each out of the 1000 draws simulates N[1] = 86
# random flips of a coin (0/1) and tallies the number of heads or "successes" (1)
(posC.sim <- rbinom(1000, N[1], p[1]))

# Examine the distribution of the number of successful pot lifts, which gives you a 
# visual summary of the expected variability of successful pot lifts:
hist(posC.sim)

# Look at the average number of pots with crab from 1000 random "realizations"
# This should be fairly close to p[1] * N[1] = 0.3 * 86 = 25.8
mean(posC.sim)   

# Now we simulate the number of crab per pot, given that a pot catches crab
# and assuming that the number of crab caught is approximately normally 
# distributed with the mean and standard deviation previously estimated
(N.sim <- round(rnorm(1000, mu[1], sd[1])))  # Simulated number of crab where present
hist(N.sim)

# Total number of crab in Bay A from 1000 simulations:
C.sim <- posC.sim * N.sim  # Simulated total number of crab -- this works because vectors are the samelength, so the vectors are just multiplied element by element
hist(C.sim)
mean(C.sim)     # Average of simulations (should be close to p*n*mu)
var(C.sim)      # Simulated variance

##### Example: Variance of a sum of two correlated random variables, x and y:

# In the previous example, we assumed that the catches in each bay are independent
# of the catches in the other bays (probably a reasonable assumption), which is a 
# necesary assumption for adding the variances:
#         var(C) = var(C_A) + var(C_B) + var(C_C)

# The next section shows you how to simulate two correlated random variables,
# which easily extends to more than 2 variables

# We then compare the variance of the two simulated variables to the
# theoretical variance based on the formula you saw in class [Mod 2 Basics(1)]:
#     var(X + Y) = var(X) + var(Y) + 2*cov(X,Y)

# You will need the following package:
library(mvtnorm)

# generate two correlated random variables
sig2.1 <- 4  # Assumed variance of variable 1
sig2.2 <- 3 # Assumed variance of variable 2
sig.12 <- 2  # Assumed covariance of variables 1 and 2

# Set up the variance-covariance matrix
(sig <- matrix(c(sig2.1,sig.12,sig.12,sig2.2),2,2))

# Exercise 1a: Compute the corresponding correlation matrix
#  (i.e. standardize variances and covariances appropriately)
# based on definition of correlation in Mod 2 Basics(1)

# Insert answer here and submit completed script file
# OR, preferably, simply paste code as text into Canvas 

# Simulate data and estimate variance, sample size 1000:
X <- rmvnorm(1000, sigma = sig)
dimnames(X)[[2]] <- c("x","y")
X

#FORMULA: var(X + Y) = var(X) + var(Y) + 2*cov(X,Y)
# Variance of the sum of the two simulated variables (x+y):
var(X[,"x"] + X[,"y"])  

# Compare to "true" variance:
sig2.1 + sig2.2 + 2*sig.12

# Exercise 1b: Repeat with a larger simulated sample size.
#              With increasing sample size, does the simulated
#              variance converge on the "true" variance

# Insert answer here and submit completed script file
# OR, preferably, simply paste code as text into Canvas 

X <- rmvnorm(1000000, sigma = sig)
dimnames(X)[[2]] <- c("x","y")
X

#FORMULA: var(X + Y) = var(X) + var(Y) + 2*cov(X,Y)
# Variance of the sum of the two simulated variables (x+y):
var(X[,"x"] + X[,"y"])  

# Compare to "true" variance:
sig2.1 + sig2.2 + 2*sig.12
