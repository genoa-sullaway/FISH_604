#########################################################

# Lab 1: Introduction to R / refresher course
# Franz Mueter

# Last modified: August 3, 2021

#########################################################
library(here)

x <- 1:10			# creates sequence 1, 2, 3, .., 10
x				# view result
y <- rnorm(10)			# creates 10 random normal numbers
y	
ls()				# to see objects currently in Data directory
args(ls)			# to see arguments that ls() can take
ls				# to see function definition

Mean(x)			# results in error message!

rm(x, y)			# removes x and y
# q()				# to exit! Alternative: Use File menu

# Getting help
help(rnorm)			# to get help file for function 'ls'
?rnorm  			# alternative
args(rnorm)			# see what arguments it takes

# Expressions and assignments
3 + 4					# Expression, output directed to screen

sqrt(3/4) * (5 - pi^2)		# Expression, note different operators
x <- rnorm(20)			# Assignment operator: output assigned to x 
x					          # view result
mean(x)				# Statistical expression 
m <- mean(x) 
v <- var(x)
m / sqrt(v)
sqrt(x)				# Note missing values in output (where x < 0)
rm(m, v)

sd				# simple function
glm				# very complex function
args(mean)			# find out what arguments a function takes

x[4] <- NA			# replace the 4th element of x by a missing value
x
mean(x)				# fails due to missing value
mean(x, na.rm =T)			# specify na.rm argument (na.rm=F by default)
mean(x, T)			# specify argument by position (without using name)
mean(x, na.rm=T, trim = 0.2)
 
x <- 1:5
sum(x);   prod(x)
cumsum(x) #cummulative sum
cumprod(x) #cummulative product 
cummax(x) 
cummin(x)
pmin(1:20, 20:1)
pmax(1:20, 20:1)
?pmin				# to learn more about pmin\pmax

round(2.3857)
round(2.3857, 2)
round(23857, -2)
floor(1.8)
ceiling(1.8)
ceiling(-1.8)
floor(-1.8)
trunc(-1.8);  trunc(1.8)	# Use semi-colon to separate expressions on one line

mean <- mean(x) 		# Can cause problems later
find("mean")		# To find each occurrence of object mean
search()			# show search list where objects are located, most contain
				# built-in (and essential) R functions
rm(mean)
?sink
sink("test.txt")			# direct output to file
"This is a test"
3 + 4					# Expression, output directed to file
mean(x)				# Another expression
sink()					# redirect output to screen
# Go to your working directory and check for file 'test.txt'

x <- rnorm(20)
mode(x)				# to get mode of x
x > 0
mode(x > 0)
5 * (x > 0)			# logical vector acts like numerical vector: 0,1
sum(x > 0, na.rm=T)		# useful for getting counts of elements meeting a
				# specified condition
y <- c("this", "is", "a", "character", "vector")		# concatenate elements
y
mode(y)
plot(0:6, 0:6, type="n")
text(1:5, 5:1, y)
c(1:5, rep(6, 5), seq(7, 20, by=2))	# Another example of concatenating

is.logical(is.character(y))		# returns TRUE

# Vectors
x <- 1:20
y <- rnorm(20)
plot(x, y)				# Plots y against x because x and y are
					# simple vectors. Uses plot.default method
x <- c(rep(1,10), rep(2,10))		# create vector with 2 values
plot(x, y)
x <- factor(x)				# convert to a factor (categorical variable)
plot(x,y)
x <- sample(1:20)			# 1:20 in random order
fit <- lm(y~x)				# linear model, regression of y on x
fit					# print out a short form of the model object

class(fit)				# show class
plot(fit)				# diagnostic plots: uses plot.lm
methods(plot)				# list all methods written for plot, produces 
# 					  different type plot for each class of objects
plot					# print function definition
plot.lm					# print plot.lm
args(plot.lm)				# show arguments that plot can take
plot(fit, 2)				#  plot only one of the graphs
summary(fit)				# uses summary.lm to extract info
plot(resid(fit))
resid(fit)
plot(coef(fit))
fitted(fit)

iris				# built-in data frame
attributes(iris)
names(iris)

# Data input \ output
GAK1.dat <- read.csv("data/GAK1.csv")
# Alternative if you have to browse for your file:
GAK1.dat <- read.csv(file.choose())

names(GAK1.dat)
GAK1.dat$depth
plot(GAK1.dat$depth, GAK1.dat$temperature)
GAK1.dat$depth <- factor(GAK1.dat$depth)
plot(GAK1.dat$depth, GAK1.dat$temperature)

# Vectors, subscripting, and missing values
x <- 1:20				# easy way to generate sequence of integers
x <- seq(0,10, length=20)
names(x)
names(x) <- paste("Element", 1:20, sep="")
x

y <- runif(20, 1, 3)			#this generates 20 values between 0 and 1
y[4] <- NA				# insert missing value at entry 4
y					# note missing value (NA)
sqrt(y)					# results in missing value in output, true for most 
# 		 			  functions that operate element-by-element
plot(y, type="l")		 	# ignores missing value, segment of line not drawn
plot(na.omit(y), type="l")
mean(y)				# results in NA
mean(y, na.rm=T)			# requires you to specifically remove missing 

cbind(x,y)				# combine two vectors of equal length into matrix
cbind(1:3, 1:20)			# cbind will happily combine vectors of equal length
					# (Shortest vector is simply repeated!)

var(cbind(x,y))			# again, by default fails because of missing value
var(y, use = "pair")		# another argument for how to handel missing values

names(x) <- letters[1:20]
x[c("a", "h", "m")]

# Matrices: 
x <- matrix(rnorm(20), nrow=5)
dim(x)
dimnames(x) <- list(NULL, letters[1:4])
x
x > 0					# matrix of logical values, all elements 0 (is a logical operator)
sum(x>0)				# total number 0, logicals are 0 (F) and 1(T) in computations
means <- apply(x, 2, mean)		# apply works on matrices (by row or column), in the second entry, a 1 indicates rows, and a 2 indicates columns
?apply
apply(x, 1, function(y) sum(y > 0))	# accepts any other function
sweep(x, 2, means)			# subtract means
?sweep

mean(x)					# scalar result
cor(x)					# matrix of pairwise corrrelations
var(x)					# var-cov matrix
row(x)					# useful in modeling, matrix computations
col(x)
lower.tri(x)
x[lower.tri(x)]				# extract elements of lower triangle
diag(x)
diag(5) #creates a 
diag(diag(x)) #creates a all zero matrix aside from the diagonal. 

X <- runif(30) #run if creates vector of 30 numbers between 0 and 1. 
dim(X) <- c(10,3)			# Alternative to matrix for creating matrices
y <- 1:10
t(X) %*% y #matrix mulitplication
crossprod(X, y)				# more transparent and efficient, same as %*% 
x <- 1:10
outer(x,y)				# outer product
x %o% y					# alternative to 'outer'

# compute arbitrary function over a grid of values. Example:  f(x,y) = cos(y) / (1+x2)
f <- function(x, y) cos(y) / (1+x^2)
a <- b <- seq(-4,4,length=50)
z <- outer(a, b, f)		
image(a, b, z)
contour(a,b,z, levels=seq(-1, 1, by=.25), 0, add=T)

# Arrays
x <- array(rnorm(20*4*3), c(20,4,3))		# create 3-dimensional array 
dimnames(x) <- list(1980:1999, c("Spring", "Summer", "Fall", "Winter"), 
	paste("Station", 1:3, sep="")) 
x
x[c("1980", "1985"), "Spring", 1:2]

apply(x, c(2,3), mean)			# long-term means by season and station
apply(x, 3, mean)				# long-term means by station

# Lists:
# (Data frames are actually stored as lists and act like lists)
x <- 1:20
y <- 10 + 3*x + rnorm(20, 0, 5)
plot(x,y)
fit <- lm(y~x)
abline(fit) #add a line to plot 
is.list(fit)
fit						# uses print.lm which does NOT print all contents of the list!
names(fit)					# prints names of all components of the list
fit[1]						# can be subscripted like a vector, extracts component of length 1
fit[[1]]					# proper way to extract components of a list!
fit$coef					# alternative for named lists
fit[["coef"]]				# alternative (useful when writing functions)
coef(fit)					# special function to extract a particular component


