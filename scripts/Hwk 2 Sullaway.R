##########################################################

# FISH 604
# Homework 2 script file
# Chatham sablefish growth
# Franz Mueter - last updated Sep 12, 2019

##########################################################

library(tidyverse)

############ Prob. 1: Data import 
sable <- read.csv("sablefish.csv", as.is=F)
head(sable)	      # Look at first few rows
hist(sable$Age)   # Age composition
table(sable$Age)  # Same in table format
summary(sable)    # Basic summary of each variable
is.factor(sable$Sex)


############ Prob. 2: Data exploration:

### 1. Plot Length-at-age, grouped by sex:
plot(Length ~ Age, data=sable, subset = Sex == "Female", col=2, 
     xlab = "Age", ylab = "Length (mm)")
points(Length ~ Age, data = sable, subset = Sex == "Male", col=4)

# ggplot version
p1 <- ggplot(data = sable, aes(x= Age, y= Length))
p1 +  geom_point(aes(color = Sex))

## Jittering coordinates to show points that overlap 

# base graphics:
plot(jitter(Length) ~ jitter(Age), data=sable, subset = Sex == "Female", 
     col=2, xlab = "Age", ylab = "Length (mm)")
points(jitter(Length) ~ jitter(Age+0.5), data = sable, 
       subset = Sex == "Male", col=4)

# ggplot version:
p1 + geom_jitter(aes(color = Sex))


### 2. Plot Length-at-age by year, grouped by sex:
sable$Year <- factor(sable$Year)

# From 'lattice' package:
xyplot(jitter(Length) ~ jitter(Age) | factor(Year), data=sable, groups = Sex)

# Females only, base graphics, with scatterplot smoother (spline):
sub <- sable$Sex == "Female"
scatter.smooth(jitter(sable$Age[sub]), jitter(sable$Length[sub]), col=2, 
               xlab = "Age", ylab = "Length (mm)")
# ggplot, females only
ggplot(data = subset(sable, Sex == "Female"), aes(x = Age, y = Length)) +
  geom_jitter() +
  geom_smooth(method = "loess")

# lattice graphics, scatterplots and scatterplot smoothers, by year
xyplot(jitter(Length) ~ jitter(Age) | factor(Year), data=sable, groups=Sex, auto.key=T)
xyplot(jitter(Length) ~ jitter(Age) | factor(Year), data=sable, groups=Sex, 
       xlim = c(0,50), ylim = c(450, 1000), auto.key=T)
xyplot(Length ~ Age | factor(Year), data=sable, groups=Sex, panel = "panel.superpose", 
       panel.groups = "panel.loess", auto.key=T)

# ggplot version, with loess smoother and 95% confidence bands:
p1 + geom_smooth(method = "loess", aes (color = Sex)) +
  facet_wrap(~Year)
# scatterplot + loess smoother:
p1 + geom_jitter(aes(color = Sex), pch=1) +
  geom_smooth(aes(color=Sex), method="loess", level=0.95) +
  facet_wrap(~Year)


########### Prob. 3: Fit Ludwig van Bertalanffy (LVB) 
###########          growth model to sablefish data 

# LvB growth model function to compute predicted values:
LvB <- function(a, k, L.inf, a0) {
  L.inf*(1-exp(-k*(a-a0)))
}

# Try out some starting values and superimpose on scatterplot:
ST <- c(k = 0.05, L.inf = 'insert a value' , a0 = 'insert a value') 
# Make sure to pick sensible values for L.inf and a0, then add line:
plot(jitter(Length) ~ jitter(Age), data=sable, col=2)
lines(1:80, LvB(1:80, k = ST[1], L.inf = ST[2], a0 = ST[3]), lwd=3)

# Fit the model using 'nls()' with the regression equation
# 'nls' finds parameters that minimize the sum of squared
# differences between the left-hand side (observed lengths)
# and the right-hand side (predicted lengths) of the formula
fit <- nls(Length ~ LvB(Age, k, L.inf, a0), data = sable, start = ST)
summary(fit)
coef(fit)

# Visualize the overall model fit
plot(jitter(Length) ~ jitter(Age), data=sable, col=2)
cf <- coef(fit)

# Add the fitted model to the scatterplot...


##### Prob. 4: Diagnostics
r <- resid(fit)  # Extract residuals (make sure 'fit' is the combined model fit)
cf <- coef(fit)
r2 <- sable$Length - LvB(sable$Age,k=cf[1],L.inf=cf[2], a0=cf[3])


# Example plot:
boxplot(r ~ Year, data = sable); abline(h=0, col=2) 


##### Prob. 5: Analyses by sex

