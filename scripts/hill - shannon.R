rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(permute)
library(geosphere)
library(iNEXT)
library(ggplot2)
library(png)
library(grid)
library(cowplot)
library(magick)
library(glmmTMB)
library(DHARMa)
library(visreg)

raw.data <- read.csv2(here("data", "Raw data.csv" ))

data.1 <- raw.data %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(Parnr)) %>% 
  mutate(Lupin.NoLupin = as.factor(Lupin.NoLupin)) %>% 
  group_by(County, Site, parnr, Lupin.NoLupin, Species) %>%
  summarise(Occ= sum(Occurrence))
data.1

# Turn into wide format
data<- data.1 %>% 
  mutate(County = dplyr::recode(County, "Örebro"= "Värmland")) %>% 
  group_by(County, Site, parnr, Lupin.NoLupin, Species) %>% 
  summarise(abund = sum(Occ)) %>% 
  rename(Lupin = Lupin.NoLupin) %>% 
  mutate(Lupin = case_when(Lupin == "Lupin" ~ "Yes", TRUE ~ "No"))
data

df <- data %>% 
  group_by(County, Site, parnr, Lupin) %>% 
  spread(., key="Species", value = "abund") %>% 
  replace(is.na(.), 0)
df 

# Divide the dataset into lupin / no lupin
lupin <- df %>% 
  filter(Lupin =="Yes") %>% 
  mutate(County= as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(lupin)

no.lupin <- df %>% 
  filter(Lupin =="No") %>% 
  mutate(County= as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(no.lupin)

# Calculate Hill diversities and extract them

# Hill diversities for plots with lupin
Hill.lupin <- renyi(lupin[,5:223], hill = TRUE)
shannon.lupin <- Hill.lupin[,4]
# Hill diversities for plots without lupin
Hill.nolupin <- renyi(no.lupin[,5:223], hill = TRUE)
shannon.nolupin <- Hill.nolupin[,4]
# Create database for test
db.1 <- bind_cols(lupin[,1:4], shannon.lupin)
db.1 <- rename(db.1, shannon = ...5)
db.2 <- bind_cols(no.lupin[,1:4], shannon.nolupin)
db.2 <- rename(db.2, shannon = ...5)

db <- bind_rows(db.1,db.2)
hist(db$shannon)
# Run analysis
model.shannon <- glmmTMB(shannon ~ Lupin + County + (1|Site), family= "Gamma", data=db) #county is a fixed effect and Site is random
summary(model.shannon)
mod_dharma1 <- model.shannon %>% simulateResiduals(n=1000)
plot(mod_dharma1)
plotResiduals(model.shannon, rank = TRUE, quantreg = FALSE)
visreg(model.shannon, "Lupin", scale="response", rug=FALSE, line = list(col="black"), 
       xlab= "Lupin present", ylab="Effective number of species (shannon) / quadrat") 

model.shannon2 <- glmmTMB(shannon ~ Lupin +  County , family= "Gamma", data=db)
summary(model.shannon2)
mod_dharma1 <- model.shannon2 %>% simulateResiduals(n=1000)
plot(mod_dharma1)
visreg(model.shannon, "Lupin", scale="response", rug=FALSE, line = list(col="black"), 
       xlab= "Lupin present", ylab="Effective number of species (shannon) / quadrat") 
