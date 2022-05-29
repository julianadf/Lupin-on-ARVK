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

# with the random subset data:
# Load and check data ----
raw.data <- read.csv2(here("data", "Raw data.csv" ))

data.1 <- raw.data %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(Parnr)) %>% 
  mutate(Lupin.NoLupin = as.factor(Lupin.NoLupin)) %>% 
  group_by(County, Site, parnr, Lupin.NoLupin, Species) %>%
  summarise(Occ= sum(Occurrence))
data.1

# use the randomly pre-selected pairs to select data
# Load the selected pairs
selected.pairs <- read.csv2(here("data", "selected pairs per site.csv"))

selected.pairs <- selected.pairs %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr))
str(selected.pairs)

# Use the selected pairs to subset the entire dataset
final.selection <- semi_join(data.1, selected.pairs, by= c("Site", "parnr"))

data <- final.selection %>% 
  mutate(County = dplyr::recode(County, "Örebro"= "Värmland")) %>%
  rename(Lupin = Lupin.NoLupin) %>%
  mutate(Lupin = case_when(Lupin == "Lupin" ~ "Yes", TRUE ~ "No")) %>% 
  mutate(County = as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  group_by(County, Site, Lupin) %>% 
  summarise(rich = n())
str(data)

model.rich <- glmmTMB(rich ~ Lupin * County + (1|Site), family= "poisson", data=data)
summary(model.rich)
car::Anova(model.rich, type= "III")
mod_dharma1 <- model.rich %>% simulateResiduals(n=1000)
plot(mod_dharma1)
testDispersion(mod_dharma1)
testDispersion(mod_dharma1, alternative = "less")
visreg(model.rich, "Lupin", scale="response", by= "County", rug=FALSE, line = list(col="black"), xlab= "Lupin present", ylab="Species richness / plot") 

# with all data:
# Load and check data ----
data <- read.csv2(here("data", "richness in pairs.csv" ))

data <- data %>% 
  mutate(Region = as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(data)

model.rich <- glmmTMB(rich ~ Lupin * Region + (1|Site), family= "poisson", data=data)
summary(model.rich)
mod_dharma1 <- model.rich %>% simulateResiduals(n=1000)
plot(mod_dharma1)
testDispersion(mod_dharma1)
testDispersion(mod_dharma1, alternative = "less")
visreg(model.rich, "Lupin", scale="response", by= "Region", rug=FALSE, line = list(col="black"), xlab= "Lupin present", ylab="Species richness / quadrat") 

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