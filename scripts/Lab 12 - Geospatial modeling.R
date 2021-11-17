##################################################################################

# FISH 604: Modern Applied Statistics
# Lab 12: A geo-spatial modeling example (kriging)

# Written by Franz Mueter, with functions provided by Margaret Short (UAF)
# Last modified: November 14, 2021

# Warning: This is a fairly advanced method and lab
#          Only apply this to your own data if you know what you are doing.  
#          This is likely to require a fair amount of additional reading 
#          on variograms and spatial data analysis

##################################################################################

# The example fits several spatial models to interpolate or smooth environmental
# data over an irregular spatial region based on a gridded sampling design

# The general steps in the analysis are as follows. These are suggestions for 
# how to approach an analysis of this sort:

# 1. Compute and plot an empirical variogram. The variogram is a description of
#    how neighboring observations covary. If you suspect a simple trend over 
#    the region, estimate variogram for trend = 'cte', '1st' and '2nd' 
#    (no trend, linear trend, quadratic trend, respectively). That means that
#    a trend of that order will be removed first before computing variogram

# 2. Select the trend (or not) that results in the most reasonable looking 
#    variogram. Don't pay too much attention to the observatons at the largest 
#    distances, which tend to be noisy as few data pairs (distances between 
#    data points) contribute to these observations and they typically have 
#    little influence on the results. What's important is the behavior near 
#    the origin (nugget, range, and how steeply it increases). 

# 3. Get reasonable starting values for parameters of the semi-variogram 
#    function (visually or trial and error): nugget, range, sill

# 4. Fit models to semi-variogram using either weighted least squares (variofit) 
#    or restricted maximum likelihood (likfit) - you can compare different 
#    autocorrelation functions if fit via REML using the AIC/BIC (because they 
#    all have the same structure). 

# 5. If you want you can look at simulation envelopes to get a feel for how 
#    uncertain the variograms really are

# 6. Compare models visually and select a "reasonable" model

# 7. Fit a geospatial model using the chosen autocorrelation function!
  
  
######### Example:

# Estimate and visualize patterns in Sea-Surface Temperature (SST) 
# in the northern Bering Sea, 2010

library(geoR)
library(deldir)
library(tcltk)
library(maps)
library(mapdata)
library(mgcv)
library(fields)
library(sp)
library(here)
source("scripts/Get borders function.R")

#####################################################################
### Fit geospatial data to environmental data (Surface temperature)
#   collected in the Northern Bering Sea in 2010:
env10 <- read.csv("data/NBS env data 2010.csv")
str(env10)

# Extract coordinates (long/lat) and surface temperatures:
dat10 <- as.matrix(env10[,c("longitude","latitude","Temp.surf")])

plot(dat10[,1], dat10[,2])   # Station locations
plot(dat10[,1], dat10[,3])   # Temperature by longitude
plot(dat10[,2], dat10[,3])   # Temperature by latitude

# Find duplicates in the location. If there are any, we need to 'jitter' 
# those locations slightly or use the mean of observations at the same 
# location because cannot have zero distanced between any two samples
any(duplicatedxy(dat10[,"longitude"], dat10[,"latitude"]))

# Determine borders of the region (convex hull around data points)
borders<- get.borders(dat10[,1:2], frac = 0.1)

plot(borders[c(1:nrow(borders), 1),],type="l")
points(dat10[,1], dat10[,2], col = "red") #check that points are within the borders

# Set up a 'geodata' object that defines the coordinates (first 2 columns) and
# the 'response' variable (column 3) and adds the convex hull as a component
# of the object (the object will be a list of class 'geodata'):
geodat <- as.geodata(dat10, coord.col = 1:2, data.col = 3)
geodat$borders <- borders

# A 'geodata' object has a plotting method ('plot.geodata', see help file)
# that generates some exploratory plots by mapping quartiles of the response 
# variable, plotting  the response against x and y coordinates, respectively,
# and plotting a histogram and density plot of the data values.

plot(geodat)

######################################################################
## Examine empirical and theoretical variograms:

# Empirical variogram. Circle size is proportional to the number of pairwise
# distances / semi-variogram values:
dat10.variog <- variog(geodat, estimator.type = "modulus", trend = "1st", 
                       uvec=seq(0, 15, length=20), bin.cloud=T)
plot(dat10.variog, bin.cloud=T)
plot(dat10.variog, main = "Robust est. of gamma(h)", pts.range = c(1,2), type='b' )

# Estimate variogram "by eye" if needed:
# Get good starting values if needed, or guess
# trial.fit<- eyefit(dat10.variog)  # Opens interactive sliders to pick parameters
# trial.fit 

# Estimate variogram using weighted least squares with selected starting values:
# The starting values ('init.cov.pars') are the sill and range parameters
# of the variogram, respectively
dat10.var.fit <- variofit(dat10.variog,
                        ini.cov.pars=c(12, 10),
                        cov.model= "exponential",
                        fix.nugget=FALSE, nugget=0, max.dist=15)

# This fits the exponential model to the binned data points using the number 
# of pairs in each bin as weights. 

# Output:
dat10.var.fit
# In the output, 'tausq' is the nugget effect ('nugget variance'), 'sigmasq'
# is the sill and 'phi' is the range paramter
lines(dat10.var.fit, col=4, lwd=2)

# Compute and plot simulation envelope
env <- variog.model.env(geodat, obj.variog = dat10.variog, 
                        model.pars = dat10.var.fit, nsim=199)  
lines(env$u, env$v.lower, col=4, lty=2, lwd=2)
lines(env$u, env$v.upper, col=4, lty=2, lwd=2)
# Note the extremely wide simulation envelope, which shows how poorly the
# empirical variogram is estimated, which is not uncommon.

# The simulated values are generated using a Gaussian random field model 
# at the data locations, given the parameters of the variogram model. The 
# empirical variogram is then computed for each simulation at the same 
# distances as the original variogram. The envelopes are computed by 
# taking, at each distance, the maximum and minimum values of the 
# variograms for the simulated data.

# Estimate variogram using restricted maximum likelihood
(ini <- dat10.var.fit$cov.pars)  # sill and range parameters
(nug <- dat10.var.fit$nugg)      # 

dat10.reml <- likfit(geodat, ini.cov.pars=ini, fix.nugget=F, trend = "1st", 
                     cov.model = "exponential", nugget=nug, lik.method="REML")
summary(dat10.reml)
plot(dat10.variog)
lines(dat10.reml, col=2, lty=2)
# This results in a very poor fit, which is not uncommon. The model estimates
# a very larger range and sill, which often indicates a spatial trend that
# is not accounted for, hence we may want to fit and remove a 2nd order 
# (quadratic) trend from the data using "trend = "2nd'": 
# This assumes that the mean for the Gaussian Random Field is a second order
# polynomial of the x and y coordinates:
# ??(x)= beta0 + beta1*x + beta2*y + beta3*x^2 + beta4*y^2 + beta5*x*y

dat10.reml2 <- likfit(geodat, ini.cov.pars=ini, fix.nugget=F, trend = "2nd", 
                     cov.model = "exponential", nugget=nug, lik.method="REML")
summary(dat10.reml2)
lines(dat10.reml2, col=4, lty=2)

# This looks much better and the AIC is much lower:
AIC(dat10.reml, dat10.reml2)



# Grid for predictions and plotting:
dxdy <- 0.25  # grid spacing in degrees
my.grid <- pred_grid(c(min(geodat$borders[,1]), max(geodat$borders[,1])),
                     c(min(geodat$borders[,2]), max(geodat$borders[,2])),
                     by=dxdy)
plot(my.grid)

# Kriging to interpolate data:
# Use parameters from REML fit)
kr.obj <- krige.control(type.krige="OK", trend.d = "2nd", trend.l = "2nd",
                           cov.model="exponential",
                           cov.pars=ini, nugget=nug)  
# Note that the 'trend' argument allows the inclusion of covariates to 
# model the mean trend as a function of other covarites (rather than lat/lon)

# Identify all locations inside the convex hull:
(gr.in <- locations.inside(my.grid, geodat$borders))
plot(gr.in)

# Make spatial predictioins over the grid:
kr.results <- krige.conv(geodat, locations=gr.in, krige=kr.obj)
# Extract standard errors and predicted values
kr.sd <- sqrt(kr.results$krige.var)
kr.pred <- as.matrix(kr.results$predict)


##### plot results:
# a. Extract longtiude and latitude values for plotting
lons <- seq(min(geodat$borders[,1]), max(geodat$borders[,1]), dxdy)
lats <- seq(min(geodat$borders[,2]), max(geodat$borders[,2]), dxdy)

# Identify which points in the overall grid are inside the 'borders':
in.border <- point.in.polygon(my.grid[,1], my.grid[,2],
                geodat$borders[,1], geodat$borders[,2])

# Set up matrices to hold output (predicted values and std. errors:
my.results <- my.sd <- matrix(NA, nrow=length(lons), ncol=length(lats))

# Replace the values that are inside the borders (1) with predicted values / sd:
my.results[in.border!=0] <- kr.pred
my.sd[in.border!=0] <- kr.sd

# Exclude points too far from data:
v1 <- rep(lons, length(lats))
v2 <- rep(lats, each=length(lons))
ex.tf <- exclude.too.far(v1, v2, geodat$coords[,1], 
             geodat$coords[,2], dist = 0.05)
my.results[ex.tf] <- my.sd[ex.tf] <- NA

# Plot 2010 predicted values:
par(mar=c(2,2,1,1))
image(lons, lats, my.results, main="", col=tim.colors(50), cex.main=1.5)
contour(lons, lats, my.results, add=T, labcex=1.2)
map('worldHires',fill=T,add=T, col="grey")
points(dat10[,1], dat10[,2], cex=0.5)

# Plot 2010 standard errors:
image(lons, lats, my.sd, main="SD", col=rev(heat.colors(50)), cex.main=1.5)
#contour(lons, lats, my.sd, add=T, labcex=1.2)
map('worldHires',fill=T,add=T, col="grey")
points(dat10[,1], dat10[,2], cex=0.5)




