############################
#
# FISH 604
# Module 2 (part 2): Discrete distributions
# Franz Mueter
# 
############################

# Probability distributions - Discrete distributions

# Data (Average Dec-Mar sea-surface temperatures in the eastern Bering Sea, 1900-2008:
SST <- scan(nlines=8)
0.740888872 0.576666654 0.863333314 0.220888884 0.446666657 0.172888885 0.317555548 -0.06911111 0.16933333 0.713333317 0.319555548 -0.127555553 -0.087333331 -0.189333329
0.264444439 0.772666649 0.230666662 0.295333327 0.486888878 0.272444438 -0.362888881 -0.663555541 0.991777756 0.559111099 0.399333324 0.386666658 0.4971111 0.65044443 1.08644442
2.015111066 0.914888868 0.863555536 0.279555549 1.220666639 0.830666648 0.716222206 1.009999977 1.232444417 -0.006444444 -0.002444444 -0.085555554 1.563777743 2.384444391 
0.439111101 1.163777752 1.704666629 1.445555523 1.496888855 1.163333307 0.65844443 0.316222215 0.808888871 0.712222206 0.229111106 -0.35977777 0.13733333 -0.226444439 0.645333319
0.423555546 0.633999986 0.463111101 0.375111103 0.588666654 1.03533331 0.644222208 0.539777766 1.143777752 0.907111091 0.506222211 0.973111089 1.085111087 -0.047555554 -0.096888887 
0.732444428 0.429333324 0.359555548 -0.135555553 1.179555529 1.333777748 1.798222182 0.919555535 0.817777759 1.525111077 0.857777759 1.221111084 2.118444397 1.456888856 1.796666627 
1.800888849 1.130888864 1.243333306 1.074222198 1.05133331 0.94711109 1.899555513 0.713777762 1.653111074 1.399555524 1.515333299 0.986888867 0.102888887 2.572444387 1.050222199 
2.271999949 1.514666633 1.84511107 0.816666648 0.883777758 0.075777776

######### Empirical probability density functions
hist(SST, freq=F)            # Histogram to approximate empirical pdf:
d <- density(SST, n = 100)   # Non-parametric density estimate to approximate empirical pdf 
lines(d,col=2)

# Use density estimate to draw empirical pdf
plot(d, type="n", xlab="Winter SST, Bering Sea", main="Empirical Probability Density") 
polygon(d, col="lightgreen")

# Highlight area corresponding to: Pr(x < 2)
polygon(c(d$x[d$x<2],2), c(d$y[d$x<2],0), col="grey")
# Highlight area corresponding to: Pr(x is between 1 and 2)
polygon(c(1, d$x[d$x<2 & d$x>1],2), c(0,d$y[d$x<2 & d$x>1],0), col="black")


###########################################################################
# Distributions in R:
# Most distributions have four functions associated with them
# that serve the following basic purposes:
# The function start with 'r', 'd', 'p', and 'q', followed by the (abbreviated) 
# name of the distribution, for example for the normal distribution:

# rnorm()  # Generate random numbers from normal distribution
# dnorm()  # Compute density (pdf) of normal distribution
# pnorm()  # Compute cumulative pdf of   "         "
# qnorm()  # Compute quantiles of normal distribution 

# We will explore these in more detail in the next lab

################## Bernoulli 
# Can be computed using binomial distribution with size = 1
# the 'rbinom' function is a random number generator for this distribution
rbinom(10, 1, prob=0.3)  # 10 random draws from a Bernoulli(0.3)

# Probability density for Bernoulli(0.3)
barplot(dbinom(c(0,1), 1, 0.3), names = c(0,1))  

# Estimated mean and variance from 1000 Bernoulli trials:
(X <- rbinom(1000, size=1, prob=0.3))
mean(X) # 'True' mean (used as basis for simulation) = p = 0.3
var(X)  # True variance = p *(1-p) = 0.21

################## Binomial
# Compute probability of getting k successes in n trials
n <- 10  # Number of trials
k <- 4   # Number of successes
p <- 0.4  # Probability of success

# Pr(k 'successes'), computed from definition of binomial pdf (Slide 10):
choose(n,k) * p^k * (1-p)^(n-k)  # binomial probability distribution
# same as: 
dbinom(k, size=n, prob=0.4)   # Probability of 4 successes
# Density functions computes probability density for any possible value
# of the random variable (k in this case)

# Probability of getting 0, 1, 2, ..., or 10 successes in 10 trials: 
dbinom(0:10, 10, prob=0.4)
# Plot probability density:
barplot(dbinom(0:10, 10, prob=0.4), nam=0:10, ylab = "Probability")

x <- rbinom(100, 10, prob=0.4)  # 100 random draws from 'binom(10,0.4)' distribution
barplot(table(x)) # Frequency distribution of random draws
# Re-run these two lines repeatedly to observe typical variability in distribution

## Normal approximation to binomial (for moderate n and p):
n <- 40   # Sample size ("trials")
p <- 0.1  # Probability of success
# First, compute and plot binomial probabilities:
y1 <- dbinom(0:30, n, prob=p)   # Probability of 0,1,2,3,...,30 successes
plot(0:30, y1, type="h", lwd=4, xlim = c(-3,25))  # Binomial probabilities
# Next, we compute and overlay the normal density (pdf) for a normal 
# distribution with the same mean (n*p) and standard deviation {sqrt(n*p*(1-p))}
# as the binomial random variable
x <- seq(-5,15,length=100) # x-values (quantiles) at which to compute density
y2 <- dnorm(x, mean = n*p, sd = sqrt(n*p*(1-p)))  # normal densities
lines(x, y2, col=2, lwd=2)  # Normal probabilities 

## Normal approximation to binomial (for "large" n):
n <- 500   # Sample size ("trials")
p <- 0.1  # Probability of success
y1 <- dbinom(0:100, n, prob=p)   # Probability of 0,1,2,3,...,30 successes
plot(0:100, y1, type="h", lwd=4, col="lightblue")  # Binomial probabilities
x <- 0:100 # x-values (quantiles) at which to compute density
y2 <- dnorm(x, mean = n*p, sd = sqrt(n*p*(1-p)))  # normal densities
lines(x, y2, col=2, lwd=2)  # Normal probabilities 
    # X axis here is the number of successes. 


# Poisson distribution (Example):

# Number of Sebastes species in 49 bottom trawl hauls 
# conducted around Kodiak Island
sp <- scan(nlines=2)
0 1 0 1 1 1 2 1 1 1 1 0 3 2 3 1 0 1 3 2 2 0 3 0 0 2 5 1 1 0 3 
1 0 1 2 1 0 1 1 0 1 1 2 2 0 1 3 0 4

# Number of hauls with k=0,1,2,.. species: 
(nr.sp <- table(sp))

# Convert to proportions (empirical probabilities):
(y <- nr.sp/sum(nr.sp))
y[c("6","7")] <- 0   # add zeros for plotting
y

# Plot probabilities
x <- barplot(as.vector(y), ylab = "Probability", cex.lab=1.6, col="slateblue")

# Distance between center of bars (for plotting):
d <- diff(x)[1]

# Fit Poisson distribution to data (maximum likelihood estimate):
library(MASS)
fit <- fitdistr(sp, "Poisson")
fit$est   # Estimated lambda!

# Predicted probabilities based on Poisson with rate lambda:
pred <- dpois(0:7, fit$est)
lines(x-d/2,pred,type="s",lwd=4)


