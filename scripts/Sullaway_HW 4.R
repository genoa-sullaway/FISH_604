##################################################################################

# SULLAWAY HWK 4

# Assignment: Based on the formulas given in the instructions, compute the model-
#           averaged estimate and the variance of the intercept of the Ricker model
#           (i.e. the Ricker productivity parameter for this salmon stock) 

###################################################################################
library(here)
library(tidyverse)
# Make sure to run through the lab 10 solutions first to get
# the object 'Branch.fits'. If you have the following file and 
# required data file (Branch sockeye.csv) in your current
# working directory, you can simply 'source' the solutions script: 
source("scripts/Lab 10 Model selection & averaging - Solutions.R")
# You may need to hit 'Enter' a couple times to advance the graphics
# Note that I 'turned off' the bootstrap as there is no need 
# to run that for the homework.


### Problem 1:

# Extract and save Ricker- a parameters:
seq <- c(1:7)
intercepts <- list()
for (i in 1:length(seq)) {
  df<- coef(Branch.fits[[i]])
  intercepts[[i]] <- df[1]
}
#change to dataframe
intercepts_df<-as.vector(purrr::map_df(intercepts, dplyr::bind_rows))

# Compute model-averaged parameter estimate:

#model averaged theta:: 
theta_hat<- intercepts_df %>%
        cbind(W.aicc) %>%  # ** if I want to do weights for aicc use column 6
        mutate(theta = `(Intercept)`* W.aicc) %>%
        summarise(sum = sum(theta))

#####ANSWER: Model averaged estimate for theta is 1.1. 
theta_hat
###############################################################################
### Problem 2:

# Extract and save standard errors of each coefficient a_i:
 seq <- c(1:7)
 var <- list()
 for (i in 1:length(seq)) {
   df<- vcov(Branch.fits[[i]])
   var[[i]] <- df[1,1]
 }
 
# Compute model-averaged variance / standard error (Equation 2)
 var_df<- unlist(var) %>%
   data.frame() %>%
   rename(variance = ".") %>%
   cbind(intercepts_df) %>%
   cbind(W.aicc) %>%
   mutate(theta_hat = theta_hat) %>%
   dplyr::mutate(var_theta =  W.aicc*(sqrt(variance + (`(Intercept)`-theta_hat)^2))) %>% 
   dplyr::summarise(sum = (sum(var_theta)^2))
 
#####ANSWER: Model averaged variance is 0.02
var_df

###############################################################################
### Problem 3:  Model averaging using 'modavg' or 'model.avg':
library(AICcmodavg)
?modavg
 Modnames <- names(Branch.fits)
 output<-modavg(parm = "(Intercept)", cand.set = Branch.fits, modnames = Modnames)
 
# I got 1.09 (rounds up to 1.1) when I did the parameter averaging by hand and I got 1.1 using modavg

#Do comparison it using model.avg
library(MuMIn)
model.avg(Branch.fits)

#####ANSWER: 
#I got the same answer using model.avg and modavg.

###############################################################################
### Problem 4:

# 95% confidence interval computed "by hand"
lower <- theta_hat - 1.96*sqrt(var_df) 
upper <- theta_hat + 1.96*sqrt(var_df) 

## ANSWER: 95% CI by hand: 0.8170338, 1.380605


# 95% confidence interval as computed by 'modavg'
output$Lower.CL 
output$Upper.CL

#ANSWER: 95% CI with  modavg is (0.8163925, 1.381246)

###############################################################################
### Problem 5:
# Model-averaged parameter estimates, standard errors, and confidence intervals:

output_lag2<-modavg(parm = "SST.t2", cand.set = Branch.fits[c(1,2,3,4)], modnames = Modnames)
output_lag3<-modavg(parm = "SST.t3", cand.set = Branch.fits[c(1,2,3,4)], modnames = Modnames)

#####ANSWER:
# Model averaging indicates that 95% of the time, the mean estimate for SST lag 2 will be between 0.11 and 0.65, modeling averaging returns and estimate of 0.38. The standard error is 0.14.  
# Model averaging indicates that 95% of the time, the mean estimate for SST lag 3 will be between -0.48 and 1.17, modeling averaging returns and estimate of 0.35. The standard error is 0.42. 
 
# What can you conclude about the effects of SST on the recruitment of 
# Branch River sockeye salmon?

#####ANSWER:
# We can conclude that SST time lagged by 2 years is important in predicting Branch river sockeye recruitment and we see a positive relationship between SST-lag 2 and recruitment.  
# However, SST lag 3 is not as an important of a predictor in recruitment because the confidence interval crosses 0. 

