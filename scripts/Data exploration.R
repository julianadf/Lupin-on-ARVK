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
raw.data <- read.csv2(here("data", "Raw data.csv" ))

# Check the data

# 1. There were duplicate species in some sites (t.ex. Gulvial in Backa). Remove duplicates
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

# With lupin
Fig.1x <- ggplot(final.selection, aes(x=fct_infreq(Species), y = Occ)) + geom_col() + facet_wrap(vars(Lupin.NoLupin), nrow =2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank()) +
  scale_x_discrete()
Fig.1x

# without lupin
data.1.nolupin <- final.selection %>% 
  filter(Species != "Lupinus polyphyllus") %>% 
  rename(Lupin = Lupin.NoLupin) %>%
  mutate(Lupin = case_when(Lupin == "Lupin" ~ "Yes", TRUE ~ "No")) %>% 
  mutate(Lupin = as.factor(Lupin))

Fig.S1 <- ggplot(data.1.nolupin, aes(x=fct_infreq(Species), y = Occ)) + geom_col() +
  facet_wrap(vars(Lupin), nrow =2) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(), 
        panel.grid.major = element_blank()) +
  ylab(bquote('Total occurrences')) + font("ylab", size=11) +
  scale_x_discrete()+ scale_y_continuous(expand = c(0,0)) 
Fig.S1 


# 2. 
data.2 <- final.selection %>% 
  mutate(Site = as.factor(Site)) %>% 
  #mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin.NoLupin)) %>% 
  group_by(County, Site,  Lupin.NoLupin) %>%
  summarise(rich = n())
data.2

lupin <- data.2 %>% 
  filter(Lupin.NoLupin =="Lupin")
summary(lupin$rich)

no.lupin <- data.2 %>% 
  filter(Lupin.NoLupin =="NoLupin")
summary(no.lupin$rich)

# Turn it into wide
data.frame <- data.1 %>%
  pivot_wider(names_from = Species, values_from = Occ) %>% 
  replace(is.na(.), 0)
data.frame


# Using all data----

# Load and check data
raw.data <- read.csv2(here("data", "Raw data.csv" ))

# Check the data

# 1. There were duplicate species in some sites (t.ex. Gulvial in Backa). Remove duplicates
data.1 <- raw.data %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(Parnr)) %>% 
  mutate(Lupin.NoLupin = as.factor(Lupin.NoLupin)) %>% 
  group_by(County, Site, parnr, Lupin.NoLupin, Species) %>%
  summarise(Occ= sum(Occurrence))
data.1

# With lupin
Fig.1x <- ggplot(data.1, aes(x=fct_infreq(Species), y = Occ)) + geom_col() + facet_wrap(vars(Lupin.NoLupin), nrow =2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank()) +
  scale_x_discrete()
Fig.1x
# without lupin
data.1.nolupin <- data.1 %>% 
  filter(Species != "Lupinus polyphyllus")

Fig.S1 <- ggplot(data.1.nolupin, aes(x=fct_infreq(Species), y = Occ)) + geom_col() +
  facet_wrap(vars(Lupin.NoLupin), nrow =2) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank()) +
  ylab("Occurrences per 1 km^2") + font("ylab", size=11) +
  scale_x_discrete()+ scale_y_continuous(expand = c(0,0))
Fig.S1 
# alternatively: (doesnt work)
Fig.S1x <- ggplot(data.1.nolupin) +
  geom_col(aes(x=fct_infreq(Species), y = Occ, fill= Lupin.NoLupin), position = position_dodge(width = 0.5)) +
  theme_bw() + ylim(0,500) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank()) + 
  scale_x_discrete() + scale_y_continuous(expand = c(0,0))
Fig.S1x 


# 2. 
data.2 <- data.1 %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin.NoLupin = as.factor(Lupin.NoLupin)) %>% 
  group_by(County, Site, parnr, Lupin.NoLupin) %>%
  summarise(rich = n())
data.2

data.3 <- data.2 %>% 
  group_by(County, Site, Lupin.NoLupin) %>% 
  summarise(n = n(), mean.rich = mean(rich), sd=sd(rich), se=sd/sqrt(n))
data.3  

# Not possible to add error bars because some sites only have one pair
Fig.S2 <- ggplot(data.3) + geom_col(aes(x=Site, y=mean.rich, fill=Lupin.NoLupin), position=position_dodge()) +
  theme_bw() + 
  #geom_errorbar(aes(x=Site, y=mean.rich, ymin=mean.rich-se, ymax=mean.rich+se), width=.2, position=position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank()) +
  ylab("Mean number of species per quadrat") + scale_y_continuous(expand = c(0,0), limits= c(0,26)) 
Fig.S2  

  

# Turn it into wide
data.frame <- data.1 %>%
  pivot_wider(names_from = Species, values_from = Occ) %>% 
  replace(is.na(.), 0)
data.frame


