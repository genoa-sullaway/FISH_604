library(tidyverse)
library(here)

#calculate mean for each year, and then get % increase
#calculate SD for start and ending time point: 
sd_1_sq = sd()
sd_2_sq = sd()
#difference of two SD's
SD = sqrt(sd_1_sq +sd_2_sq)
#SE= SD/sqrt(sample size)

POP <- read.csv("data/Pacific Ocean Perch.csv")

pop_lm<-lm(Biomass/1000~Year, data=POP)
summary(pop_lm)

new.dat <- data.frame(Year = 2013)
predict(pop_lm, newdata = new.dat, interval = 'confidence')

new.dat <- data.frame(Year = 1984)
predict(pop_lm, newdata = new.dat, interval = 'confidence')

#confint(pop_lm)

plot(pop_lm$residuals, col = "red")
#library(car)
# res <- leveneTest(Biomass/1000~Year, data=POP)
# res

preditct_pop<-predict(fit.lm, se=T)
sum(preditct_pop$se.fit)
#for every increase in year the biomass goes up by 26715 (the slope) of the lm, the survey period is 29 years,
  #so for the survey period that means that the biomass went up by 774.735.

#634  1501 

#Std error of year = 5.917, so the annual change in biomass can vary by 5.917

#(I will caveat this answer by saying I dont think a linear model is the best way to answer this question)


5.917*29=171

