##########################################################

# FISH 604
# Homework 1 script file

# Expected values and variances: Sockeye salmon fecundity

# Genoa Sullaway
 
##########################################################
library(here)
library(tidyverse)
library(moments)
#### Problem 1:

### Import data:
f.dat <- read.csv("data/sockeye fecundity.csv", head=T)
str(f.dat)           # Examine data
table(f.dat$Lake)    # Number of observations by lake

### (a) Plot histograms and boxplots

#---- Enter your code below and briefly annotate as needed ----
#plot histrogram of fecundity
hist(f.dat$Fecundity) 
#Fecundity looks like it is normally distributed, with a slight right skew. 

#plot histogram, group fecundity by lake. Looks like a slight difference in fecundity distributions by lake, Lake 2 looks like it has higher fecundity. 
hist.lake<-ggplot(f.dat, aes(x=Fecundity, color=Lake)) +
  geom_histogram(fill="white") +
  theme_classic()
print(hist.lake)

#plot boxplot of data - again we see slight right skew with some outliers around 6000. 
ggplot(f.dat, aes(y=Fecundity)) +
  geom_boxplot(fill="white") +
  theme_classic()

#similarly it is a little more apparent that there is a higher median fecundity in Lake 2. 
ggplot(f.dat, aes(x=Lake, y=Fecundity, color=Lake)) +
  geom_boxplot(fill="white") +
  theme_classic() +
  theme(legend.position = "none") 
#--------------------------------------------------------------

# For further examining the data, it is useful to split the 
# fecundity vector into two components ('split' creates a list)
(f.dat2 <- split(f.dat$Fecundity, f.dat$Lake))


### (b) Compute mean fecundities (expected values) and their standard errors:

#---- Enter your code below and briefly annotate as needed ----

f.dat.sum<-f.dat %>%
  group_by(Lake) %>%
  summarise(mean_fecundity = mean(Fecundity), sd.Fecundity= sd(Fecundity), n=length(Fecundity), 
            se=sd.Fecundity/sqrt(n))

print(f.dat.sum)
#--------------------------------------------------------------


### (c) Compare means 
# Check out help files for 't.test' and 'oneway.test'
# ?t.test
# ?oneway.test

#---- Enter your code below and briefly annotate as needed ----
t.test(f.dat$Fecundity~f.dat$Lake, var.equal=FALSE)

#Results: Mean fecundity for sockeye salmon in lake 2 (Mean = 3958.024, SE= 46.1) is greater than fecundity for sockeye in 
#Lake 1 (Mean = 3389, SE= 44.2) (p<0.05). 

#I think separate estimates of fecundity and variance should be reported results should be reported because mean fecundity varies by lake additionally the samples are independent, pooling the results across lake would cover up inherent variability in the sockeye population.
 
#--------------------------------------------------------------

###############################################################
#### Problem 2:

### Input values needed for expected values and variance:

# Proportion of females by lake:
p1 <-  0.58
p2 <- 0.52

# Expected numbers of total spawners:
E.N1 <- 124360
E.N2 <- 73608

# Standard errors:
s.N1 <- 23245
s.N2 <- 10120

#convert standard errors to sd
#sd = sqrt(SE^2 * N)
sd.N1 <-sqrt(s.N1^2 * 265)
sd.N2 <- sqrt(s.N2^2 * 254)

# Expected mean fecundity
E.f1= 3389
E.f2= 3958
  
#Sfi, standard error of fecundity
s.f1 = 44.2
s.f2 = 46.1

#std_deviation of fecundity 
sd.f1=720
sd.f2=734

#Var = SE^2 * N
#need to convert standard error to std dev or variance 

### (a) Expected number of eggs (in millions):
# Note that pi is a constant, and Ni and fi are independent

#---- Enter your code below and briefly annotate as needed ----
#E(Eggs)=E(E1+E2)=E(p1 N1 f1 + p2 N2 f2)
E.Eggs1=p1*E.N1*E.f1 
E.Eggs2=p2*E.N2*E.f2
E.Eggs = E.Eggs1 + E.Eggs2

# We expect 395 million eggs across both lakes (395,941,544). 
#--------------------------------------------------------------



###  (b) Variance:

#---- Enter your code below and briefly annotate as needed ----
# I feel like there are two ways to do this... the way franz wrote it out where Var(X+Y) = Var(X) + Var(Y) or the formula
        # var(E1 + E2) = var(E1) + var(E2) =E(X^2) -E.fi = ((E.Eggs1^2) -E.fi1) + ((E.Eggs2^2) -E.fi2)

# where Var(X+Y) = Var(X)+Var(Y)+2Cov(X,Y) <--- *** this = 0 if variables are uncorrelated. So OG formula is the same! aha. 
# E.var.eggs=(E.Eggs1 -E.f1)^2 + (E.Eggs2 -E.f2)^2
# sqrt(E.var.eggs)

#ORRR THIS WAY: ---- check later to make sure everything is inputted correctly - very hard to keep track
# var(E1 + E2) 
# = var(E1) + var(E2)
# = var(E(p1*N1*f1) + (p2 N2*f2)) # p = constant, f and E = random variables
# = p1^2 * E(N1*f1) + p2^2 * E(N2*f2)

#breaking down the above formula ^^

varN1=sd.N1^2 #sd.N1^2 *E.N1^2 
varf1=sd.f1^2 #sd.f1^2*E.f1^2
E.N1f1= E.f1^2*varN1+ E.N1^2*varf1+ varN1*varf1

varN2= sd.N2^2 #s.N2^2 *E.N2^2
varf2= sd.f2^2 #s.f2^2*E.f2^2
E.N2f2 = E.f2*varN2+E.N2*varf2+varN2*varf2

E.var.eggs = p1^2 * E.N1f1 + p2^2 * E.N2f2

#I still get really big numbers!! ?!?!
 

#--------------------------------------------------------------



### (c) 95% confidence interval, assuming normal distribution:

#---- Enter your code below and briefly annotate as needed ----
#Remember that 95% of a normal distribution with mean μ and standard deviation σ is contained within the interval {μ - 1.96σ, μ + 1.96σ}.

# 95% Confidence Interval for total catch assuming normal distribution:

cat(round(E.Eggs - 1.96*sd.C), "to",  round(E.Eggs + 1.96*sd.C), "Sockeye Eggs")

#--------------------------------------------------------------


### Problem 3
#3a #(a) In one sentence, define skewness (3 pts)
#---- Enter your code below and briefly annotate as needed ----
 #Skewness is the degree to which the distribution of your data are uneven.  
#--------------------------------------------------------------

#3b Based on the definition in (a) or examining the formula for skewness, what is the expected skewness (average across many samples) of a sample from a normally distributed population? (Justify your answer). 
#You may also want to check your answer by repeatedly generating large random samples from a normal distribution (use 'rnorm()') and computing their skewness
#(use, for example, function 'skewness()' in package 'moments'). (3 pts).

#--------------------------------------------------------------
# The skewness of a normal distribution should be 0 because the cubed sum of individual observations from the mean will be zero or very close to 0.

#check 
N.sim <- rnorm(100000)
hist(N.sim)
moments::skewness(N.sim)

#Confirmed - the skewness of a simulated normal distribution is very close to zero. 
#--------------------------------------------------------------





#3(c) Compute the skewness for the two sets of fecundity measurements in lake 1 and lake 2 (separately) and compare qualitatively.
#Is there any evidence of skewness? (3 pts)
#--------------------------------------------------------------
#calculate skew for each lake
lake.skew <- f.dat %>%
  group_by(Lake) %>%
  summarise(moments::skewness(Fecundity))

print(lake.skew)
#Skew calculations for both lakes are larger than 0, indicating a right skew. 
#However, the right skew is larger for fecundity in Lake 1 (0.961) compared to Lake 2 (0.382). 
#--------------------------------------------------------------






#3(d) Follow the bootstrap example estimate the uncertainty in the skewness measure for each of the two lakes and construct 
#95% confidence intervals based on the 2.5th and 97.5th percentiles of the bootstrap distribution. 
#Does the 95% confidence interval include the expected value for a normal distribution? 
#What can you conclude about skewness of fecundity in the two lakes? (10 pts)
#--------------------------------------------------------------
lake.skew <- f.dat %>%
  group_by(Lake) %>%
  summarise(moments::skewness(Fecundity))

#there is a quicker way to do this with apply but I was getting weird answers 

#f<-sapply(f.dat2, function(x) sample(x,replace = T))

# for(i in 1:999) { #loop for bootstrapping
#   f<-sapply(f.dat2, function(x) sample(x,replace = T),simplify = "array")
#   
#   out[i]<-sapply(f, function(x) skewness(x),simplify = "array")
# }

lake1.out <- numeric(999)  #vector for the output
lake2.out <- numeric(999)
for(i in 1:999) { # loop for bootstrapping
  lake1 <- sample(f.dat2[[1]], replace=T)
  lake1.out[i] <- moments::skewness(lake1)

  lake2 <- sample(f.dat2[[2]], replace=T)
  lake2.out[i] <- moments::skewness(lake2 )
   }

mean(lake1.out) #take mean of bootstrapped skewness outputs
mean(lake2.out)

#calcualte confidence intervals for the bootstrapped outputs 
cat(mean(lake1.out) - 1.96*sd(lake1.out), "to",  mean(lake1.out) + 1.96*sd(lake1.out), "Lake 1 Skewness 95% CI")
cat(mean(lake2.out) - 1.96*sd(lake2.out), "to",  mean(lake2.out) + 1.96*sd(lake2.out), "Lake 2 Skewness 95% CI")

#The 95% CI does not include the expected value for a normal distribution, which would be zero (though the lower end of the confidence interval for lake 2 is close to 0).
#So we can conclude that the fecundity for both lakes has a right skew and that Lake 1 has a larger right skew than Lake 2.
#--------------------------------------------------------------






#3e
#(e) Finally, use the 'agostino.test()' (package: moments) to test for skewness. 
#Compare conclusions to those in (d) (5 pts)
#--------------------------------------------------------------
agostino.test(as.numeric(L[[1]]$Fecundity))
agostino.test(as.numeric(L[[2]]$Fecundity))
print(lake.skew)

#The calculations provide similiar outputs, thus using both the skewness calculation and the agostino test I can conclude that the data do have right skew. 
#--------------------------------------------------------------







