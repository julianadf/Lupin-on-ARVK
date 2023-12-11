rm(list=ls())
library(here)
library(tidyverse)
library(emmeans)
library(vegan)
library(betapart)
library(lmerTest)
library(cowplot)
library(knitr)
library(pairwiseAdonis)
library(png)
library(magick)
library(ggpubr)
library(dendextend)
library(betapart)

# Load data
df <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))

# Restructure
df.res <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  dplyr::select(-Rubus.fruticosus) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region =as.factor(Region))  
str(df.res)

lupin.res <- df.res[, c(4:ncol(df.res))]
comm <- df.res[, c(5:ncol(df.res))]
lupin <- lupin.res[,1]

# Abundance-based multiple site-dissimilarities
bc.dist <- beta.pair.abund(comm, index.family = "bray")
bd.bc <- betadisper(bc.dist[[3]],lupin, type = "centroid")
plot(bd.bc, col= c("#008837", "#7b3294"), main = NULL, sub= NULL, lwd=2, cex=1.5)

boxplot(bd.bc)
anova(bd.bc)
# calculate distance between centroids
sqrt(dist(bd.bc$centroids[,bd.bc$eig>0])^2 - dist(bd.bc$centroids[,bd.bc$eig<0]^2))

#Uppland
df.upp <- df.res %>% 
  filter(Region == "East region")
df.upp

comm.upp <- df.upp[, c(5:ncol(df.upp))]
lupin.upp <- df.upp[,4]

# Abundance-based multiple site-dissimilarities
bc.dist.upp <- beta.pair.abund(comm.upp, index.family = "bray")
bd.bc.upp <- betadisper(bc.dist.upp[[3]],lupin.upp, type = "centroid")
plot(bd.bc.upp, col= c("#008837", "#7b3294"), main = NULL, sub= NULL, lwd=2, cex=1.5)

boxplot(bd.bc.upp)
anova(bd.bc.upp)


#Värmland
df.värm <- df.res %>% 
  filter(Region == "West region")
df.värm

comm.värm <- df.värm[, c(5:ncol(df.värm))]
lupin.värm <- df.värm[,4]

# Abundance-based multiple site-dissimilarities
bc.dist.värm <- beta.pair.abund(comm.värm, index.family = "bray")
bd.bc.värm <- betadisper(bc.dist.värm[[3]],lupin.värm, type = "centroid")
plot(bd.bc.värm, col= c("#008837", "#7b3294"), main = NULL, sub= NULL, lwd=2, cex=1.5)

boxplot(bd.bc.värm)
anova(bd.bc.värm)

# Regions
comm <- df.res[, c(5:ncol(df.res))]
region <- df.res[,1]

# Abundance-based multiple site-dissimilarities
bc.dist <- beta.pair.abund(comm, index.family = "bray")
bd.bc <- betadisper(bc.dist[[3]], region, type = "centroid")
plot(bd.bc, col= c("#d7191c", "#2c7bb6"), main = NULL, sub= NULL, lwd=2, cex=1.5)

boxplot(bd.bc)
anova(bd.bc)


# Presence absence ------

# calculate the Jaccard index and its partitions of turnover and nestedness
# Based on presence / absence data
presabs.lupin <- ifelse(lupin.res[,5:ncol(lupin.res)]>0,1,0)
lupin <- lupin.res[,1]
dist.pa <- beta.pair(presabs.lupin, index.family="jaccard")
bd <- betadisper(dist.pa[[3]],lupin)
plot(bd, main="sorensen") #The betadisper plot indicates that there is NO difference in species compositions (as the NMDS)

# boxplot to observe the distance of values of beta diversity of each treatment in relation to their
# centroids (basically, this indicates homogeneity in how communities of a given treatment differ 
# from each other). 
boxplot(bd)
# Perform an ANOVA to test if treatments are significantly different.
anova(bd) #almost, but treatments do not differ significantly in relation to how communities vary from each other.

beta.tot <- beta.multi(presabs.lupin)
gg.beta <- as.data.frame(beta.tot)
# beta.SOR = total, beta.SIM = turnover, beta.NES = nestedness


# Uppland
df.upp <- df.res %>% 
  filter(Region == "East region")
df.upp

lupin.upp <- df.upp[, c(4:ncol(df.upp))]

# Based on presence / absence data
presabs.lupinupp <- ifelse(lupin.upp[,5:ncol(lupin.upp)]>0,1,0)
upp.lupin <- lupin.upp[,1]

# calculate the Jaccard index and its partitions of turnover and nestedness
dist.upp <- beta.pair(presabs.lupinupp, index.family="jaccard")
bd.upp <- betadisper(dist.upp[[3]],upp.lupin)
plot(bd.upp)
boxplot(bd.upp)
anova(bd.upp) #treatments do not differ significantly in relation to how communities vary from each other.

beta.totUpp <- beta.multi(presabs.lupinupp)
gg.betaUpp <- as.data.frame(beta.totUpp)


# Värmland 
df.värm <- df.res %>% 
  filter(Region == "West region")
df.värm

lupin.värm <- df.värm[, c(4:ncol(df.värm))]

# Based on presence / absence data
presabs.lupinvärm <- ifelse(lupin.värm[,5:ncol(lupin.värm)]>0,1,0)
värm.lupin <- lupin.värm[,1]

# calculate the Jaccard index and its partitions of turnover and nestedness
dist.värm <- beta.pair(presabs.lupinvärm, index.family="jaccard")
bd.värm <- betadisper(dist.värm[[3]],upp.lupin)
plot(bd.värm)
boxplot(bd.värm)
anova(bd.värm) # treatments do not differ significantly in relation to how communities vary from each other.

beta.totvärm <- beta.multi(presabs.lupinvärm)
gg.betaVärm <- as.data.frame(beta.totvärm)

# plot the 3 betas (all, uppland värmland)
# join the 3:
gg.beta$which <- "Both regions"
gg.betaUpp$which <- "East region"
gg.betaVärm$which <- "West region"

gg.temp <- merge(gg.beta, gg.betaUpp, all=TRUE)
gg.betaLupin <- merge(gg.temp, gg.betaVärm, all=TRUE)

# turn to long 
gg.lupin <- gg.betaLupin %>% 
  gather(beta.type, value, beta.SIM:beta.SOR) %>% 
  mutate(beta.type = dplyr::recode(beta.type, "beta.SOR"= "Total")) %>% 
  mutate(beta.type = dplyr::recode(beta.type, "beta.SIM"= "Turnover")) %>% 
  mutate(beta.type = dplyr::recode(beta.type, "beta.SNE"= "Nestedness"))
gg.lupin

gg.beta <- gg.lupin %>% filter(beta.type != "Total") %>% ggplot(aes(x = which, y = value, fill = beta.type)) + 
  geom_bar(stat = 'identity') +
  scale_y_continuous("Total beta-diversity\nbetween plot types") +
  scale_x_discrete("") +
  scale_fill_manual("", values = c("#969696", "#cccccc")) +
  theme_minimal() 
  #theme(axis.text.x=element_blank(), text = element_text(size=16))
gg.beta


# Beta diversity between regions -----
df <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))

# Restructure
df.res <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  dplyr::select(-Rubus.fruticosus) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region =as.factor(Region))  
str(df.res)

region.res <- df.res[, c(1,5:ncol(df.res))]

# Based on presence / absence data
presabs.lupin <- ifelse(region.res[,5:ncol(region.res)]>0,1,0)
region <- region.res[,1]

# calculate the Jaccard index and its partitions of turnover and nestedness
dist <- beta.pair(presabs.lupin, index.family="jaccard")
bd.region <- betadisper(dist[[3]],region)
plot(bd.region) #The betadisper plot indicates that there is NO difference in species compositions (as the NMDS)

# boxplot to observe the distance of values of beta diversity of each treatment in relation to their
# centroids (basically, this indicates homogeneity in how communities of a given treatment differ 
# from each other). 
boxplot(bd.region)
# Perform an ANOVA to test if treatments are significantly different.
anova(bd.region) #Regions do differ significantly in relation to how communities vary from each other.



# Based on abundance data ----
dist2 <- bray.part(lupin.res[,5:ncol(lupin.res)])
bd2 <- betadisper(dist2[[3]],lupin)
plot(bd2)
boxplot(bd2)
anova(bd2)

abund.lupin <- lupin.res[,5:ncol(lupin.res)]
lupin <- lupin.res[,1]

dist <- bray.part(abund.lupin)
#list 1: bray.bal = balanced variation = turnover (species replacement)
#List 2: bray.gra = abundance gradient = nestedness 
#list 3: bray =  total abundance-based dissimilarity between sites
bd <- betadisper(dist[[3]],lupin)
plot(bd) 

turnover <- betadisper(dist[[1]],lupin)
boxplot(turnover)

 boxplot to observe the distance of values of beta diversity of each treatment in relation to their
# centroids (basically, this indicates homogeneity in how communities of a given treatment differ 
# from each other). 
boxplot(bd)
# Perform an ANOVA to test if treatments are significantly different.
anova(bd) #almost, but treatments do not differ significantly in relation to how communities vary from each other.

beta.tot <- beta.multi.abund(abund.lupin)
gg.beta <- as.data.frame(beta.tot)





# plot it:
tot.beta <- gg.beta %>% filter(Type != "Total") %>% ggplot(aes(x = lupin, y = beta, fill = lupin)) + 
  geom_bar(stat = 'identity') +
  scale_y_continuous("Total beta-diversity\nbetween habitat types") +
  scale_x_discrete("") +
  scale_fill_manual("", values = c("#c7eae5", "#01665e")) +
  theme_minimal() +
  theme(axis.text.x=element_blank(), text = element_text(size=16))
tot.beta

