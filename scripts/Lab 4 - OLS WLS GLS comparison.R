#########################################################

# Lab 4 (part 1): Least-squares regression: OLS WLS GLS comparison

# Franz Mueter   Last modified: 9/03/2021

# Comparing ordinary, weighted, and generalized least squares

# requires data in 'Lab 4 OLS-WLS-GLS example.csv'

# See also spreadsheet 'OLS WLS GLS comparison.xlsx'
#Group: Queen salmon
#########################################################

# This script compares three approaches for fitting a model 
# via least squares (LS) regression

# 1. OLS: Ordinary Least Squares is used when all regression assumptions are
# met, i.e. observations are "independent and identically distributed"
# As the name (LS) implies, a model is fit to data by finding the parameter 
# values that minimize the residual sum of squares
#    (See: Mod 4 Estimation I.pdf,  slides 9-11)

# 2. WLS: Weighted Least Squares is used when the variances are NOT equal.
# The approach works by putting more weight on observations that have a smaller
# variance and less weight on observations that have a larger variance (i.e. we 
# are less certain about their values). If variances are known, optimum weights
# to use in the regression are the inverse of the variance of each observation.
# Often we don't know the variances but can estimate them as a function of the 
# mean response, e.g. the standard deviation or the variance may increase linearly
# with the mean. We can then use optimum weights based on this relationship.
# Models are fit to data by finding the parameter values that minimize the 
# weighted residual sum of squares
#   (See: Mod 4 Estimation I.pdf,  slides 13-21)

# 3. GLS: Generalized Least Squares is used when both the equal variance assumption
# and the independence assumption may be violated (or when only the independence 
# assumption is violated). 

# The approach works by setting up a variance-covariance matrix that accounts for 
# the dependence structure in the residuals. There are many possible forms for this 
# dependence. For example, time series data may have a dependence of successive data 
# points on each other (first-order autocorrelation) with no correlation among 
# observations more than one time step apart. Other examples include spatial 
# autocorrelation that decreases as a function of the geographical distance 
# between observations.

# If an appropriate dependence structure can be identified and accounted for, models
# are fit to the data by finding parameter values that minimize the "generalized residual
# sum of squares", which are squared residuals modified by their correlation structure.
# The correlation structure itself is parameterized (e.g. as a function of distance) and
# these parameters are estimated along with the usual regression parameters.
#   (See: Mod 4 Estimation I.pdf,  slides 22-29)

# We will use the 'nlme' package for generalized least squares regression 
#  (functions gls, gnls). We will use the same package later for fitting mixed
#   effects models. 'nlme' stands for 'non-linear mixed-effect' but the package
#   can fit a variety of regression models)
library(nlme)  

########## Import data to work with:
dat <- read.csv("Lab 4 OLS-WLS-GLS example.csv")
dat
# The data comprise a short time series linking a response (y) that was measured at 
# times t=1 through t=25 to a single predictor variable (x). The data file already
# includes appropriate weights that reflect our confidence in the individual measurements
# (smaller variance = more confidence = larger weights)

# We already encountered the data earlier (script 'Mod 4 non-linear OLS example.R').
# Let's plot the response against the predictor variable:
plot(y~x, data=dat, pch=16)

# Looking at the data again shows a curvilinear relationship between 
# y and x that may be fit by a power model: y = a*x^b, where the power 
# coefficient b is less than 1 (based on experience):

# An example of such a model is shown here (with arbitrary parameters):
curve(5*x^0.5)

# We will fit models using the three different approaches to these data, 
# starting with Ordinary Least Squares regression

###################################################################################
########## OLS: 
# Here we use the non-linear least-squares function 'nls' to fit a non-linear model
# of the form y = a*x^b to the data by minimizing the residual sum of squares.
# Note that we need to specify starting values for parameters 'a' and 'b'
# (in the form of a list):

# (See also 'Mod 4 non-linear OLS example.R')

fit.OLS <- nls(y ~ a*x^b, data=dat, start = list(a=2.5, b=0.1))
summary(fit.OLS)
# We can extract the coefficients:
(cf <- coef(fit.OLS))

# The regression output should look familiar and includes the estimates for the
# parameters along with their ("asymptotic") standard deviations and a t-test
# that simply tests the null hypotheses H0: a=0 and H0: b=0. Not surprisingly, 
# both 'a' and 'b' are significantly different from zero (which is not very 
# informative as the null hypotheses a=0 and b=0 are not useful!

# Let's check the fit by visualizing the fitted curve:
plot(y~x, data=dat, pch=16)
curve(cf[1]*x^cf[2], add=T, lwd=2, col=2)
# Here cf[1] and cf[2] correspond to a and b in our model


###################################################################################
########## WLS 
# Next we fit the same model using a weighted least squares regression. Note that
# the residuals from the previous fit do not show a pattern that would suggest 
# unequal variances, as can be seen in a standard diagnostic plot of residuals
# against fitted values:
plot(fitted(fit.OLS), resid(fit.OLS)) 

# However, we have prior knowledge suggesting differences in our confidence in each
# observation and we assumes that the Weights in the data frame are inverse 
# proportional to the variance of each observation. For example, they could represent 
# sample sizes if each observation in y is the mean of a number of measurements  

# We can use the same 'nls' function to fit the model using WLS. The function 'nls',
# like most modeling functions, has a 'weights' argument. Beware that the 'weights' 
# are defined differently in, for example, lm and nls compared to gls and gnls
# The functions lm and nls use weights that are INVERSE PROPORTIONAL to the known or
# assumed variance: w = 1/s^2), whereas weights in 'gls' and 'gnls' are PROPORTIONAL
# to the assumed variances

# Let's re-fit the model but this time giving different weights to each observation
# (based on the weights ('Weights') included in the data frame, which - again - are
# inverse proportional to the variance of each observation)
fit.WLS <- nls(y ~ a*x^b, data=dat, start = list(a=2.5, b=0.1), weights = Weights) 
summary(fit.WLS)
# The results are very similar to those from the OLS regression above but note that
# the standard errors of the regression coefficients 'a' and 'b' are slightly smaller.
# That is not uncommon if we chose our weights wisely to down-weight observations
# that we have less confident in (implying that these observations are "noisier"
# with a higher variance)!

# We can also fit the model using the function 'gnls', which we will use below to fit
# the model via Generalized Least Squares (GLS, see below). The results using the default
# settings for'gnls' (assuming independent residuals) are the same as in 'fit.WLS' above,
# but note that the weights are specified differently:
fit.WLS2 <- gnls(y ~ a*x^b, data=dat, start = list(a=2.5, b=0.1), weights = ~ 1/Weights)	
summary(fit.WLS2)
# Note that this is identical to using the 'varFixed' structure - see '?varClasses':
gnls(y ~ a*x^b, data=dat, start = list(a=2.5, b=0.1), weights = varFixed(~1/Weights))

# Be careful with weights in any regression package to make sure you know how weights
# need to be specified. As noted earlier, the functions 'gls' and 'gnls' use the assumed 
# or known variances (= 1/Weights) instead of the corresponding weights. This is because
# gnls (& gls) can use a variety of modeled variance that may be a function of the
# predicted values or of other covariates. If given as a formula, it assumes 
# that the variance of the observations is proportional to the fixed vales on the 
# right hand side of the equation (see help files for 'gls' and 'varClasses')

# For example, if we have reason to believe the variance depends on the 
# x-values (or some other covariate), we can specify a power function for 
# the variance structure such that the weights are a function of x:
fit.WLS3 <- gnls(y ~ a*x^b, data=dat, start = list(a=5, b=0.1), 
                 weights = varPower(form = ~x))	
summary(fit.WLS3)

fit.WLS3.test <- gnls(y ~ a*x^b, data=dat, start = list(a=4, b=0.5), 
                 weights = varPower(form = ~x))	
summary(fit.WLS3.test)

# The parameter estimate for the power parameter is -0.314, but its 
# confidence interval includes 0 (which would imply equal weights):
intervals(fit.WLS3)
# hence there is little evidence that the variance changes with x

# You can examine the "implied" variance weights as follows:
vf <- varPower(value = -0.314, form = ~x)
vf <- Initialize(vf, dat)
plot(dat$x, varWeights(vf))  
# where the weights range 'moderately' from about 1.8 to 3.4

# These weights are set up to be inverse proportional to the estimated
# standard deviations, which are therefore:
1/varWeights(vf)

# In general, be careful with weights and only include them if you
# have a good reason and either use 'external' (known) weights or 
# estimate weights based on a known relationship that is well
# supported by the data

# See also the script 'Mod 4 WLS & variance structures.R'

#########################################################################################
########## GLS 
# Finally, we re-fit the model using both weights and a first-order autoregressive 
# correlation structure. 

# Our observations represent a time series and the residuals from the above model
# fits show some evidence of autocorrelation over time:
plot(dat$time, resid(fit.WLS), pch=16, type="b"); abline(h=0)

# We can formally test for autocorrelation using the Durbin-Watson test (covered later),
# which is implemented in the package 'lmtest':
library(lmtest)
dwtest(resid(fit.WLS)~1)  
# Note that 'dwtest()' takes a formula or linear model object as its first argument.
# Since we already fit a non-linear model to the data, we simply test the residuals
# from that model using a '1' on the right hand side (i.e. "intercept only" model, 
# although in this case the intercept is zero by definition). 
# Note that the autocorrelation is not significant (p  = 0.134)

# Nonetheless, "to be safe" we may want to fit a model that includes autocorrelated 
# residuals. This can be implemented in both 'gls' (for linear models) and 'gnls' 
# (for non-linear models) using the argument 'correlation'. This argument can take 
# a variety of functions that represent different "correlation classes". Here we 
# assume the residuals have a first-order autoregressive structure (i.e., the residual
# at time t depends on the residual at time t-1 in the way shown in class: 
# Mod 4 Estimation I.pdf, slide 25)
fit.GLS <- gnls(y ~ a*x^b, data=dat, start = list(a=2.5, b=0.1), 
    weights = ~ 1/Weights, correlation = corAR1(form = ~ time))	
summary(fit.GLS)
# The output shows the familiar coefficient estimates, their standard errors, and t-tests
# for testing the null hypotheses that a=0 and b=0, respectively (here the p-values are
# rounded to 0). In addition, we get an estimate of the autoregressive coefficient, which is
# rather small and is not significantly different from 0 as we can see if we look at the 95%
# confidence intervals for all parameter estimates (including residual std. error):
intervals(fit.GLS)
# Note that the confidence interval for the first-order autoregressive coefficient 'phi'
# includes zero, consistent with the Durbin-Watson test that there is no significant
# autocorrelation at the 95% confidence level (see above)

# Note that the results differ somewhat from the 'Solver' (Excel) solution in the 
# spreadsheet ('Lab 4 OLS WLS GLS.xlsx') because different algorithms are used! 
# In particular, the 'phi' parameter differs quite a bit. In general, the phi parameter 
# is difficult to estimate, has large uncertainty (large confidence intervals), and 
# can vary substantially depending on the algorithm and its assumptions.

# In this example, 'gnls' converges on the same solution from a variety of 
# starting values, including "bad" values! 
# (Try fitting the model above (line 166/67) using 'bad' starting values
# for the 'a' and 'b' parameters)

# Plot fit:
plot(y ~ x, dat)
curve(coef(fit.OLS)[1]*x^coef(fit.OLS)[2], col=2, add=T)
curve(coef(fit.WLS)[1]*x^coef(fit.WLS)[2], col=3, lty=2, lwd=2, add=T)
curve(coef(fit.GLS)[1]*x^coef(fit.GLS)[2], col=4, lty=3, lwd=2, add=T)
legend("bottomright", c("OLS", "WLS", "GLS"), col=2:4, lty=1:3)

# Note that the weighting makes a difference (compare green WLS fit to red OLS fit), 
# whereas including autocorrelation has little influence on the fit (compare blue GLS
# fit to green WLS fit). This is not surprising, because the level of autocorrelation 
# was small (and not significant at 95%) as we saw above!

# This concludes comparison of OLS, WLS, and GLS. For the mathematical details of how
# the algorithms work, please look at the spreadsheet 'OLS WLS GLS comparison.xlsx'


#######################################################################################
# Construct confidence intervals for the parameters of the model

# Assuming "asymptotic normality" (meaning that the distribution of parameters will 
# tend towards a normal distribution, given large enough sample sizes), we construct
# approximate 90% and 95% confidence intervals (CI) for the parameters for each model. 

# First, we compute CIs directly from the output of the model:
#    theta +/- z * SE(theta)
# Here, 'theta' is a parameter of the model ('a' or 'b'), z is a multiplier corresponding
# to the appropriate quantile of the normal distribution that will result in the desired
# confidence interval, and SE(theta) is the standard error of the parameter. 

# You may (should?) recall that the multiplier 'z' is 1.96 for a 95% CI, meaning that
# a normally distributed random variable with mean theta has a 95% probability of 
# falling within the interval: theta +/- 1.96 * SE(theta). 

# More precisely, in the 'frequentist' interpretation, we say that with 95% probability
# the true parameter value falls within the interval. 

# Recall that you can compute the multiplier for any normal confidence interval 
# using the 'qnorm()' function, which returns the quantile (q) of the standard 
# normal distribution corresponding to a given cumulative probability that random 
# variable x (with mean 0 and standard deviation 1) is less than or equal to q:
qnorm(0.025) # Pr(x <= -1.96) = 0.025
qnorm(0.975) # Pr(x <= 1.96) = 0.975

# In addition to the multiplier, we need the variance (or standard error) of the parameters
# to compute the confidence interval. We can extract the variance-covariance matrix 
# for the parameters of almost any model using the 'vcov()' function, for example:
vcov(fit.OLS)
# The diagonal of this matrix has the VAIRANCES for the two model parameters. 
# We can extract the diagonal using:
diag(vcov(fit.OLS))

# Note that the 'summary' output also includes the same standard errors of the parameters. 
# We can extract those standard errors from the output using the 'coefficients'
# component of the summary output, for example:
summary(fit.WLS)$coef
# This is a matrix from which you can simply extract values, for example by name:
summary(fit.WLS)$coef["b","Std. Error"]

# These are roughly the same as the square root of the variances:
sqrt(diag(vcov(fit.OLS)))
# The difference is due to 


#----------------------------------------------------------------------------------------
# Exercise 1:
# Given the above guidance on computing z and extracting parameter variances, write simple
# expressions to compute the 90% and 95% confidence intervals for the parameters of the 
# model for each of the three approaches! Write general code so you can simply replace
# 'fit.OLS' with 'fit.WLS' or 'fit.GLS' to compute the corresponding confidence intervals.

ci<-function(model, model_name,parameter, ci_level){ 
                    theta <- summary(model)$coef[parameter, "Estimate"]
                    se <- summary(model)$coef[parameter,"Std. Error"]
                    
                        upper <- theta + ci_level * se 
                        lower <- theta - ci_level * se
                        ci<- cbind(model_name,parameter,lower,upper)
                        return(ci) }
#90%
qnorm(0.05) # Pr(x <= -1.96) = 0.025
qnorm(0.95)

#95%
qnorm(0.025) # Pr(x <= -1.96) = 0.025
qnorm(0.975)

ci(model = fit.OLS,model_name="fit.OLS.model", parameter = "a", ci_level = 1.96)
ci(model = fit.OLS,model_name="fit.OLS.model", parameter = "b", ci_level = 1.96)


ci(model = fit.OLS,model_name="fit.OLS.model", parameter = "a", ci_level =  1.644854)
ci(model = fit.OLS,model_name="fit.OLS.model", parameter = "b", ci_level =  1.644854)


ci(model = fit.WLS,model_name="fit.WLS.model", parameter = "a", ci_level = 1.96)
ci(model = fit.WLS,model_name="fit.WLS.model", parameter = "b", ci_level = 1.96)


ci(model = fit.WLS,model_name="fit.WLS.model", parameter = "a", ci_level =  1.644854)
ci(model = fit.WLS,model_name="fit.WLS.model", parameter = "b", ci_level =  1.644854)

#do it differently for GLS because variance is in the var cov matrix:

  theta_a <- summary(model)$coef["a", "Estimate"]
  se<- diag(vcov(fit.GLS))
  ci_level = 1.96
  upper <- theta_a + ci_level * se[1] 
  lower <- theta_a - ci_level * se[1]
  ci_GLS_95_a<- cbind(lower,upper)

  
  theta_b <- summary(model)$coef["b", "Estimate"]
  se<- diag(vcov(fit.GLS))
  ci_level = 1.96
  upper <- theta_b + ci_level * se[2] 
  lower <- theta_b - ci_level * se[2]
  ci_GLS_95_b<- cbind(lower,upper)
  
#90%
  theta_a <- summary(model)$coef["a", "Estimate"]
  se<- diag(vcov(fit.GLS))
  ci_level = 1.644854
  upper <- theta_a + ci_level * se[1] 
  lower <- theta_a - ci_level * se[1]
  ci_GLS_90_a<- cbind(lower,upper)
  
  
  theta_b <- summary(model)$coef["b", "Estimate"]
  se<- diag(vcov(fit.GLS))
  ci_level = 1.644854
  upper <- theta_b + ci_level * se[2] 
  lower <- theta_b - ci_level * se[2]
  ci_GLS_90_b<- cbind(lower,upper)
  

  # Exercise 2:
# Compare the width of the confidence intervals between the WLS (weights only)
# and GLS fits (weights + autocorrelation parameter).
# How do they differ? How might you expect them to differ based on the 
# example in class? (Mod 4 Estimation I - OLS WLS GLS.df, slide 29)
# component of the summary output, for example:

ci(model = fit.WLS,model_name="fit.WLS.model", parameter = "a", ci_level = 1.96)
ci_GLS_95_a

ci(model = fit.WLS,model_name="fit.WLS.model", parameter = "b", ci_level = 1.96)
ci_GLS_95_b

# The width of the WLS is larger and has greater uncertainty than the GLS, these differ because the GLS accounts for non-independence in 
# the response variable and by accounting for this increases our certainty in the estimates. 


# Exercise 3:
# For several classes of models in R you can simply use the 'confint()' function 
# to compute confidence intervals for any desired parameter and at a specified 
# significance level. Check the help function and compute 
# CIs using 'confint'! 

# How do the 'confint' intervals compare to those you computed above?
ci(model = fit.OLS,model_name="fit.OLS.model", parameter = "a", ci_level = 1.96)
ci(model = fit.OLS,model_name="fit.OLS.model", parameter = "b", ci_level = 1.96)

?confint
confint(fit.OLS, level = 0.95)
confint(fit.WLS, level = 0.95)
confint(fit.GLS, level = 0.95)

#These barely differ, we think the difference is due to rounding in the ci_level (our function rounds it). 

# Exercise 4 (optional, if time allows):
# To quantify uncertainty in the fitted values, we would ideally construct 
# confidence bands around the fitted line based on the uncertainty in predicted 
# values. For most models, we can use the 'predict' function to compute 
# predicted values and their standard errors, which can then be used to 
# construct confidence bands (as you saw for a linear model in an earlier 
# example in 'Mod 3 Transformations.R'.
# For non-linear models, confidence intervals for parameters and confidence
# bands for predicted relationships can be misleading and incorrect as they
# are based on linear model theory. Therefore, computing standard errors for
# predicted values from 'nls' has not even been implemented in the correpsonding
# 'predict.nls' funtion and we would typically use different approaches for 
# such models, for example bootstrapping (as you will see later). An alternative
# is to use the delta method, which is implemented in the package 'propagate'
# for models fit via 'nls'
library(propagate) 

# The package has a function 'predictNLS' (which is documented on the author's blog:
# https://rmazing.wordpress.com/2013/08/14/predictnls-part-1-monte-carlo-simulation-confidence-intervals-for-nls-models/
# The help file for predictNLS has the basics but is a bit scrambled.

# Use 'predictNLS' with either fit.OLS or fit.WLS (or both) to construct
# simulated confidence bands and overlay the results on the scatterplot
# of y vs x

# Because this is not very well documented (unpublished), I typically prefer 
# either a Bayesian or a bootstrap approach for non-linear models and we will 
# do so later in the semester



