rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(lattice)
library(ggpubr)
library(effects)
library(visreg)
library(DHARMa)
library(MASS)
library(lme4)
library(car)
library(png)
library(cowplot)
library(visreg)
library(pscl)
library(glmmTMB)


# Load and check data ----
data <- read.csv2(here("data", "richness in pairs.csv" ))

data <- data %>% 
  mutate(Region = as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(data)

model.rich <- glmmTMB(rich ~ Lupin + Region + (1|Site), family= "poisson", data=data)
summary(model.rich)
mod_dharma1 <- model.rich %>% simulateResiduals(n=1000)
plot(mod_dharma1)
testDispersion(mod_dharma1)
testDispersion(mod_dharma1, alternative = "less")
visreg(model.rich, "Lupin", scale="response", rug=FALSE, line = list(col="black"), xlab= "Lupin present", ylab="Species richness / quadrat") 

#test whether there are differences between regions to be sure. Final model is above
model.rich2 <- glmmTMB(rich ~ Lupin + Region + (1|Site), family= "poisson", data=data)
summary(model.rich2)
mod_dharma1 <- model.rich2 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
visreg(model.rich, "Lupin", scale="response", rug=FALSE, line = list(col="black"), xlab= "Lupin present", ylab="Species richness / quadrat") 

# test autocorrelation of model residuals (only if requested, but here is a start)

library(ncf)

# data are in data frame df, with X and Y columns containing longitude and latitude
# data were analysed with model m

# join coordinate data and model data

coords <- read.csv2(here("data", "sites with coordinates.csv" ))

coords <- coords %>% 
  mutate(N = as.numeric(N)) %>% 
  mutate(E = as.numeric(E))
str(coords)


ncf.spl.res <- spline.correlog(coords$E, coords$N, resid(model.rich, type  = "pearson"), resamp=1000, na.rm = T, latlon= F)
plot(ncf.spl.res) # seems like i need the coords in latitude longitude format or it wont work (these are in SWEREF 99TM)