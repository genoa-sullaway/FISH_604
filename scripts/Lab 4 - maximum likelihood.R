######################################################################

# Lab 4 (part 2): Maximum likelihood approach to regression: 
#                 Poisson regression example

# Franz Mueter   Last modified: 9/7/2021

######################################################################


# In this lab exercise we will work through both a "brute-force" approach to
# maximum likelihood estimation (MLE) and using the function 'glm' to estimate 
# parameters via MLE

# The 'brute-force' approach conducts a grid search over a range of reasonable
# parameter values to find the combination of regression parameters that maximize
# the likelihood of the observations. 

# We will use the Poisson regression example that I illustrated in class
# (Mod 4 Estimation II - Maximum Likelihood.*, slides 14-21)

# The example looks at the rate at which sharks are caught as bycatch in 
# a tuna longline fishery. We have observations from 100 longline sets
# that counted the number of sharks per (standardized) set. Each set also 
# recorded the  depth at which the gear was set and we are interested in 
# exploring if and how the bycatch rate differs with depth (so we may be 
# able to minimize bycatch by fishing at depths with lower bycatch rates) 

# Based on an examination of the data we model the number of sharks per 
# (standardized) longline set as a function of depth with a quadratic term
# to allow for a dome-shape relationship (as evident in the data)

# Because the response ('sharks') consists of relatively low counts (including
# some zeros), which cannot be negative, a regression model with a normal
# error distribution is NOT appropriate here (it could result, for example, in
# negative estimates of bycatch). Moreover, if we look at the data (see below)
# we note that there is a higher variability in the counts when the counts are
# large. This suggests that we can use a Poisson regression. Recall that a Poisson
# variable can take positive integer values (including zero) and that the 
# variance of a Poisson increases in proportion to the mean count (higher 
# counts have a larger variance).

# Therefore we assume a Poisson distribution for the shark data and model 
# the rate parameter 'lambda' (corresponding to the average bycatch rate) 
# as a quadratic function of depth (This model may be described as a 
# Quadratic Poisson regression)

# First, read in and plot the data:
depth <- c(412,429,387,271,100,50,405,456,450,349,202,265,325,381,485,
	390,313,107,121,470,257,111,451,114,136,261,49,31,182,339,252,10,
	339,217,299,257,405,277,429,261,369,336,382,226,479,114,261,435,
	164,294,377,275,379,56,215,63,367,131,412,28,384,106,484,66,14,342,
	159,321,378,31,32,331,40,390,301,282,292,159,445,284,125,252,359,
	221,91,120,364,51,270,69,180,463,372,241,201,407,172,198,328,230)
sharks <- c(1,3,4,6,2,2,3,2,2,3,3,4,5,2,3,4,6,1,3,1,7,6,1,5,4,7,4,4,9,4,
	3,2,1,9,2,6,6,9,1,5,2,7,4,7,1,7,4,3,6,4,0,3,3,2,3,4,4,2,3,2,8,10,
	0,3,2,1,11,3,3,6,4,7,3,6,6,5,5,6,4,6,5,5,1,11,10,3,5,5,10,3,5,1,2,
	2,4,6,4,3,8,3)

plot(depth, sharks)
# Note the apparent peak in bycatch rate at intermediate depth (hence the
# quadratic form), the relatively low counts, and the higher variance when
# counts are higher (ranging from ~0-2 sharks for deep sets where bycatch
# rates are low to ~2-11 sharks at intermediate depths where rates are high)


### Model:
# We fit a Poisson regression of the average bycatch rate 'lambda' as a 
# quadratic function of depth (depth and square of depth as predictors, to 
# allow for a dome-shaped relationship). 

#         lambda = a + b1*depth + b2*depth^2

# This is equivalent to a quadratic regression in the general linear model case
# (where we would assume that the observed counts follow a normal distribution
#  with mean lambda and some variance)

# However, here we do NOT assume normal residuals, but we assume that
# the observed counts follow a Poisson distribution with mean lambda


### Data
# The response variable consists of integer counts (number of sharks per 
# longline set) and our only explanatory variable in this simple example is 
# depth (with linear and quadratic term)


### Fitting the model
# We can fit the model using a 'generalized linear model' (GLM) approach, which 
# we will cover in more detail later. A Poisson regression is one example of
# a GLM and is appropriate here because rare events such as catching one
# or more (or no) sharks in a longline set are often Poisson distributed.

# To fit the model we use the 'glm' function. We will cover this function 
# and details about this model in the module on "Generalized Linear Models".
# Because we fit a quadratic model, let's first create a new variable that 
# represents squared depth:
depth2 <- depth^2

# To fit the model using glm we use the exact same kinds of formulas that we use
# when fitting linear models via 'lm' (after all, GLMs are linear models too).
# The difference is that we now assume a Poisson distribution instead of a normal
# distribution for the data. This is specified by the 'family' argument! You can 
# ignore the 'link' for now - we will cover what that means later!
fit.pois <- glm(sharks ~ depth + depth2, family=poisson(link="identity"))
summary(fit.pois)

# The output looks very similar to the output from 'lm' and the components 
# should mostly be familiar: We have three parameter estimates corresponding 
# to 'a' (Intercept) 'b1' (depth) and 'b2' (depth2) in our model (see above)

# There are some differences that we will discuss later, most notably that we 
# now use a 'z' test (instead of a t-test) for testing the null hypotheses 
# (H0:a=0, H0:b1=0, or H0:b2=0). This is because we assume that the variance 
# for the Poisson is known (= lambda), rather than estimated, as in a 
# regression with normal errors, where we estimate the variance from the 
# residuals, which introduces additional uncertainty

# The z-tests suggest that both the linear and quadratic terms are highly 
# significant (that is, they are significantly different from zero. Of primary 
# interest is the quadratic term, which suggests that there is indeed a dome-
# shaped  relationship and that bycatch rates of sharks are significantly 
# higher at intermediate depths! (The fact that 'b1' is also significantly 
# different from zero cannot be separated from the quadratic term in this 
# case - if it were not significant we would still include it IF the 
# quadratic term is significant).

# We can extract the estimated coefficients (a, b1, and b2, respectively):
coef(fit.pois)


#############################################################################
### Fitting the same model using a grid search (the 'brute force' approach)

# To illustrate the maximum likelihood approach, we will find the parameter
# values that maximize the likelihood (or, rather, minimize the negative log-
# likelihood) by doing a grid search, that is we compute and plot the negative 
# log-likelihood over a grid of possible values that the regression parameters 
# b1 and b2 may reasonable take (for illustration, we assume for the moment 
# that the intercept 'a' is known)

# First, we write a function to compute the negative log-likelihood (logL) for 
# a Poisson variable (See 'Mod 4 Estimation II - Maximum Likelihood.*' slides
# 7 & 17 for the appropriate likelihood formulas)

# The function computes logL as a function of the slope parameters b1 and b2 
# (here included as a single vector of length 2 that we call 'b'), the 
# intercept ('a'), the observed counts (a vector of counts 'k' of length n) 
# and the observed depth variables (vectors x1 and x2 for depth and depth^2). 
# Our function first computes the expected count 'lambda' (as a quadratic
# function of depth), then computes the corresponding negative log-likelihood 
# as a function of lambda and k. Both 'lambda' and 'logL' will be vectors of 
# length 'n' with one value for each observation. The last line of the function
# sums up the individual log-likelihoods for each observation to get the total 
# log-likelihood. The sum on the last line is what is returned when you use 
# the 'lik.pois' function with data as you'll see below 
lik.pois <- function(b, a, k, x1, x2) {
	lambda <- a + b[1]*x1 + b[2]*x2
	logL <- lambda - k*log(lambda) + log(factorial(k))
	sum(logL)
}

# For example, we can compute the total log-likelihood for a=2.5, b1=0.035, 
# and b2=-0.00005 as follows:
lik.pois(b=c(0.035,-0.00005), a=2.5, k = sharks, x1=depth, x2=depth2)

# Similarly, we can calculate the likelihood for any desired combination 
# of the parameters a, b1, and b2

# To compute the likelihood over a range of trial values, we set up a grid of 
# reasonable values for b1 and b2, ranging from 0.02 to 0.04 for b1 and 
# from -0.00008 to -0.00004 for b2 (I restricted the parameter space over
# which to search based on prior knowledge of the likely range of values, 
# which you may have to determine through trial and error)
(b1 <- seq(0.02, 0.04, length=50))
(b2 <- seq(-0.00008, -0.00004, length=50))

# b1 and b2 are vectors. We can use 'expand.grid' to create a 2-column 
# matrix that includes each combination of b1 and b2:
b <- expand.grid(b1,b2)
b # Examine the output to see how 'expand.grid' works!
names(b) <- c("b1", "b2")

# We can plot all of the parameter combinations (the 'grid') at which we 
# will evaluate the log-likelihood:
plot(b, pch=16, cex=0.4)  
# Each point in the plot corresponds to one unique combination of b1 and b2
# (or one row of the matrix 'b')


# Now we are ready to use our function (lik.pois) to compute the negative 
# log-likelihood for each combination of b1 and b2 (i.e. for each row 
# in the matrix b)

# For simplicity (to restrict the problem to maximizing the likelihood over 
# two parameters rather than three), we fix 'a' at the ML estimate that we 
# obtained above (there is really no reason not to estimate 'a' - I'm only 
# fixing it for illustratio so that I can plot the likelihood surface 
# in 2-D space for b1 and b2)

a.est <- coef(fit.pois)[1]  
a.est  # intercept 'a' as estimated by glm above

# We can compute the negative log-likelihood for each parameter combination 
# (row) in b using the 'apply' function (a very useful function that you 
# should familiarize yourself with!): 'apply' uses our function 'lik.pois' 
# with the trial values for b1 and b2 in a given row and the values for the 
# other quantities as specified (k, x1, x2, a) to compute the negative 
# log-likelihood for each row:
lik <- apply(b, 1, lik.pois, k=sharks, x1=depth, x2=depth2, a = a.est)
lik # negative log-likelihood for each row in b!

# Note that there are parameter combinations that result in "NaN" - these are
# non-permissible combinations that result in a negative rate parameter (lambda)
# and therefore in non-sensible results.

# Remember that our goal is to find the parameter combination (b1,b2) that 
# minimizes the negative likelihood. You could scroll through the results to 
# find it (good luck) or we can plot the likelihood surface and find the 
# minimum graphically or numerically

# To do so, we first convert the vector of likelihoods to a matrix, then use
# the 'image' and 'contour' functions to draw a color contour plot:
z <- matrix(lik,50)   
image(b1,b2,z, col=heat.colors(100))
contour(b1,b2,z,add=T)
# The resulting surface (colors) shows the values of the negative log-likelihood 
# at each combination of b1 and b2 over the selected range of b1 and b2
# (excluding the impossible parameter combinations)

# If you wish, you can add a legend using the 'image.plot' function in the 
# 'fields' library:
library(fields)
image.plot(b1,b2,z, col=heat.colors(100))
contour(b1,b2,z,add=T)

# The contours show parameter combinations that have the same negative 
# log-likelihood, here shown at intervals of 20 likelihood units

# We can superimpose our maximum likelihood estimate that we obtained by glm:
points(coef(fit.pois)[2],coef(fit.pois)[3], pch=16, cex=2)


# We covered a large range of likelihood values, hence it is hard to clearly  
# see the low point in the likelihood surface that minimizes the negative 
# likelihood (large 'reddish area'). For a clearer display, let's limit our 
# view to only the likelihood values near the minimum:  
z[z > 230] <- NA  # Arbitrarily set all likelihood values over 230 to missing
image(b1,b2,z, col=heat.colors(100))
contour(b1,b2,z,add=T)
# This makes it a litte easier to visualize the "valley" in the likelihood 
# surface. This "valley" runs from the upper left to the lower right. There 
# are a number of # combinations of b1 and b2, that all result in a negative 
# log-likelihood less than 212.

# Note also that the parameter estimates for b1 and b2 are negatively 
# correlated. That means we will get a very similar fit (similar likelihood)
# if we increase b1 by a small amount while simultaneously decreasing b2 
# (i.e. we stay on the "valley floor")

# The "steepness" of the valley relates to how certain we are about the 
# parameter values. The steeper the likelihood surface around the minimum, 
# the smaller the confidence intervals for the parameters will be.

#------------------------------------------------------------------------------------------
# Exercise 1:
# Compute the magnitude of the correlation between parameters 'b1' and 'b2'
# based on our original Poisson model fit (fit.pois)

# To do so, you need the variance / covariance matrix, which is estimated by 
# the numerical algorithm based on the steepness of the gradients in the 
# vicinity of the parameter values that maximize the likelihood. 

# As an aside, the gradients in the likelihood surface are captured by the 
# "Hessian matrix", which is the second derivative of the log-likelihood 
# function with respect to the parameters. The covariance matrix that 
# quantifies parameter uncertainty is the inverse of the negative of the 
# Hessian matrix [some algorithms will only return the Hessian matrix!]

# Recall from part 1 that you can extract the covariance (=var-cov) matrix 
# from a model object using the 'vcov()' function

# Find the covariance between the parameters b1 and b2 (for depth and depth^2,
# respectively) using the definition of the correlation between two variables
# x and y, which is given by:
#     cov(x,y) / sqrt(var(x)* var(y))
#  (see Module 2(part 1) on Expectations and Variances)

# Does the estimated correlation agree with the patterns in the likelihood 
# surface for the two parameters that we examined above?

#------------------------------------------------------------------------------------------

# Let's find the combination of b1 and b2 among all of the "trial" values 
# that minimizes the negative log-likelihood (We can use a logical vector 
# to identify the row in 'b' where 'lik' is equal to the minimum):
(b.hat <- b[which.min(lik),])
# This is the row (that is, the combination of b1 and b2) that had the 
# smallest value for the negative log-likelihood and hence the values in
# this row are the Maximum Likelihood Estimates for b1 and b2:
b.hat   

# We can compare the results from our grid search to the estimates for b1  
# and b2 that we obtained using the 'glm' function:
coef(fit.pois)

# Note that the coefficients differ slightly because we searched only over
# a fixed grid of discrete values

# Finally, we plot the fitted model(s):
plot(depth, sharks, xlab="",ylab="", cex.axis=1.5)

x <- 1:500  # Set up x variable to span range of depths
b1 <- as.numeric(b.hat[1])   # Need to convert to a scalar to avoid matrix multiplication below!
b2 <- as.numeric(b.hat[2])
# Compute estimated mean catches (catch rate lambda as a function of depth):
y.hat <- a.est + b1*x + b2*x^2
lines(x, y.hat, lwd=3, col=4)

# Alternatively, you can simply use the curve function:
plot(depth, sharks, xlab="",ylab="", cex.axis=1.5)
curve(a.est + b1*x + b2*x^2, add=T, col=4, lwd=3)


###############################################################################
# Use numerical algorithms to find Maximum Likelihood Estimates
# (Gradient search)

# In some cases, the grid search is very useful, in particular to explore 
# complicated likelihood surfaces. More often, if we can't use 'glm' (or some 
# other model function) to estimate parameters (as is the case with any 
# "non-standard" models), we can use a number of different functions to find 
# the parameter combination that minimizes the negative log-likelihood (or,
# equivalently, maximizes the likelihood)

# One such function is 'optim()'. By default, it finds the minimum of a 
# specified function (also called the 'objective function')

# To work with 'optim' (for estimating all three parameters simultaneously) 
# we need to rewrite the likelihood function slightly. Note that 'b' is now 
# a vector of length 3 consisting of values for 'a' (b[1]), 'b1'(b[2]), and 
# 'b2'(b[3]):
lik.pois2 <- function(b, k, x1, x2) {
	lambda <- b[1] + b[2]*x1 + b[3]*x2
	logL <- lambda - k*log(lambda) + log(factorial(k))
	sum(logL)
}

# As for many numerical methods, we need to provide starting values for the 
# three parameters. The algorithm first computes the negative log-likelihood 
# for these starting values and then proceeds to change the values slightly 
# using some scheme until the 'objective function' (the negative log-likelihood 
# in this case) cannot be reduced any further (i.e. reaches its minimum):
b.hat <- optim(c(2, 0.1, 0.1), lik.pois2, k=sharks, x1=depth, x2=depth2)
b.hat
# The result is a list that includes the parameter estimates ('par') and the
# value of the objective function (logL = 211.4).

# You can compare the estimates from 'optim' and 'glm'
b.hat$par
coef(fit.pois)

# Note that the estimated coefficients are not identical although both 
# 'optim' and 'glm' maximize the likelihood over all three parameters. 
# This is due to differences in the estimation procedure used. The 
# negative log-likelihood values are very similar:
b.hat$value
-logLik(fit.pois) # Extract (negative) log-likelihood from 'glm' fit!
# Note that 'glm' achieved a slightly lower negative log-likelihood

# We can compare the fitted lines for the two estimates:
bb <-b.hat$par # extract parameter values
# Add line for 'optim' solution:
curve(bb[1]+ bb[2]*x + bb[3]*x^2, add=T, col=2, lty=2, lwd=2)


# The difference is due to the fact that 'optim()' found a local minimum 
# or it stopped too early (stopping rules can be modified using different 
# control settings, see help file for 'optim')

# To guard against local minima, we can use a number of different 
# starting values. Let's try one:
optim(c(2.5, 0.05, 0.001), lik.pois2, k=sharks, x1=depth, x2=depth2)
# This is much closer to the 'glm' solution, but still has a slightly 
# larger value for logL

# Note that if we use the 'glm' results as starting values, 'optim()' 
# immediately converges to the same solution obtained by glm:
optim(coef(fit.pois), lik.pois2, k=sharks, x1=depth, x2=depth2)

# Hence be careful when using 'optim' or similar algorithms!

# The modeling function such as 'glm' have a variety of ways to ensure
# that they converge to a true minimum, but there are no guarantees!!

#------------------------------------------------------------------------------
# You may want to try some very different starting values for 'optim()' to see 
# if you can get the algorithm to converge! What is the log-Likelihood of the 
# 'worst' solution you can find that still seemingly converges? 
# What are the corresponding parameter estimates?
#------------------------------------------------------------------------------
