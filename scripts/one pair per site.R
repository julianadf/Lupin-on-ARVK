rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(permute)
library(geosphere)
library(iNEXT)

# NMDS ----
# Pool abundances per site (not average, because it doesnt really matter for NMDS)
raw.data <- read.csv2(here("data", "Raw data.csv" ))

data.1 <- raw.data %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(Parnr)) %>% 
  mutate(Lupin.NoLupin = as.factor(Lupin.NoLupin)) %>% 
  group_by(County, Site, parnr, Lupin.NoLupin, Species) %>%
  summarise(Occ= sum(Occurrence))
data.1

# Select random pair in sites with more than one
onepair.data <- data.1 %>%
  group_by(Site) %>% 
  #summarise()
  slice_sample()
onepair.data

#save the pairnumbers into a csv to use again later
write.csv2(onepair.data[,1:3], "selected pairs per site.csv")
# Load the selected pairs
selected.pairs <- read.csv2(here("data", "selected pairs per site.csv"))

selected.pairs <- selected.pairs %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr))
str(selected.pairs)

# Use the selected pairs to subset the entire dataset
final.selection <- semi_join(data.1, selected.pairs, by= c("Site", "parnr"))

# Turn into wide format
df <- final.selection %>% 
  group_by(County, Site,Lupin.NoLupin) %>% 
  spread(., key="Species", value = "Occ") %>% 
  replace(is.na(.), 0)
df 

plants <- df[,5:188]
Lupin.NoLupin <- df[,4]

NMDS <- metaMDS(plants, distance="bray", k=2, try = 100, trymax = 1000, autotransform=TRUE)
NMDS
stressplot(NMDS)

##pull points from MDS
NMDS1 <- NMDS$points[,1] ##also found using scores(NMDS)
NMDS2 <- NMDS$points[,2]
plant.plot<-cbind(plants, NMDS1, NMDS2, Lupin)

# quick plot ordination
p<-ggplot(plant.plot, aes(NMDS1, NMDS2, color=Lupin.NoLupin))+
  geom_point(position=position_jitter(.1), shape=3)+##separates overlapping points
  stat_ellipse(type='t',size =1)+ ##draws 95% confidence interval ellipses
  theme_minimal()
p

# PERMANOVA ----
# NUll hypothesis: Groups do not differ in spread or position in multivariate space.
plant.dist <- vegdist(plants, method='bray') #abundance based
env <- plants[,1:3]

# Does the community composition differ in plots with/without lupin?
# set the blocks (sites)
plant.div2 <- adonis2(plant.dist ~ Lupin.NoLupin, data=df , method="bray")
plant.div2
