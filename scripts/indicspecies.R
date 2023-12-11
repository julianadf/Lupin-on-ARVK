# Species-habitat association analysis using the "indicspecies" package
rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(ggpubr)
library(indicspecies)

df <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))
abund <- df[,5:ncol(df)] %>% dplyr::select(-Lupinus.polyphyllus)
Lupin <- df$Lupin
Region <- df$County

associations <- multipatt(abund, Region, func = "r.g", control = how(nperm=9999))
summary(associations)

# Run analysis for each region separately
# Uppland

Uppland <- df %>% filter(County =="Uppsala")
abund <- Uppland[,5:ncol(Uppland)] %>% dplyr::select(-Lupinus.polyphyllus)
Lupin <- Uppland$Lupin

associations <- multipatt(abund, Lupin, func = "r.g", control = how(nperm=9999))
summary(associations)


# Värmland
Värmland <- df %>% filter(County =="Värmland")
abund <- Värmland[,5:ncol(Värmland)] %>% dplyr::select(-Lupinus.polyphyllus)
Lupin <- Värmland$Lupin

associations <- multipatt(abund, Lupin, func = "r.g", control = how(nperm=9999))
summary(associations)


