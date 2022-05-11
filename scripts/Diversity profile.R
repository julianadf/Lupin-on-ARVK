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
library(iNEXT)

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

# turn into wide data
selected.pairs <- final.selection %>% 
  group_by(County, Site, Lupin.NoLupin) %>% 
  spread(., key="Species", value = "Occ") %>% 
  replace(is.na(.), 0) %>% 
  rename(Lupin = Lupin.NoLupin) %>% 
  mutate(Lupin = case_when(Lupin == "Lupin" ~ "Yes", TRUE ~ "No"))
selected.pairs

# Divide the dataset into lupin / no lupin
lupin <- selected.pairs %>% 
  filter(Lupin =="Yes") %>% 
  mutate(County= as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(lupin)

no.lupin <- selected.pairs %>% 
  filter(Lupin =="No") %>% 
  mutate(County= as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(no.lupin)

# Calculate Hill diversities and plot them

# Hill diversities for plots with lupin
Hill.lupin <- renyi(lupin[,5:188], hill = TRUE)
plot(Hill.lupin)
Hill.lupin$Lupin <- "Yes"

# Hill diversities for plots without lupin
Hill.nolupin <- renyi(no.lupin[,5:188], hill = TRUE)
plot(Hill.nolupin)
Hill.nolupin$Lupin <- "No"

# Join the two results
Hill.data <- bind_rows(Hill.lupin, Hill.nolupin)

# Turn into tidy data 
Hill <- Hill.data %>% 
  gather(., key="Alpha", value="Diversity", 1:11,na.rm=TRUE) %>%
  mutate(Lupin = as.factor(Lupin))
Hill

# Calculate a mean per landscape type and select alpha 0,1,2,inf
HillMean <- Hill.data %>% 
  gather(., key="Alpha", value="Diversity", 1:11,na.rm=TRUE) %>%
  mutate(Lupin = as.factor(Lupin)) %>%
  group_by(., Lupin, Alpha) %>%
  summarise(Diversity=mean(Diversity)) %>% 
  filter(Alpha %in% c("0", "1", "2", "Inf"))
HillMean

# Calculate standard error per landscape type and select alpha 0,1,2,inf
HillSD <- Hill.data %>% 
  gather(., key="Alpha", value="Diversity", 1:11,na.rm=TRUE) %>%
  mutate(Lupin = as.factor(Lupin)) %>%
  group_by(., Lupin, Alpha) %>%
  summarise(SD=sd(Diversity)/sqrt(length(Diversity))) %>% 
  filter(Alpha %in% c("0", "1", "2", "Inf"))
HillSD

HillMean$SD <- HillSD$SD

# Diversity plot for all alpha values
gg.hill <- ggplot(Hill, aes(x=Alpha, y=Diversity))+ geom_point(aes(colour = Lupin))  + 
  ggtitle("Hill diversity") + theme_bw() +
  scale_color_manual("Lupin", values = c("green", "purple"))
gg.hill

# Plot alpha values for species richness (0), Shannon (1), Simpson (2), Berger-Parker (inf) for mean category
ggMean <- ggplot(HillMean, aes(x=Alpha, y=Diversity, group = Lupin)) + 
  geom_line(linetype="dashed") + 
  geom_point(aes(colour = Lupin), size=4) + 
  geom_pointrange(aes(ymin=Diversity-SD, ymax=Diversity + SD, colour=Lupin)) +
  geom_errorbar(aes(ymin=Diversity - SD, ymax= Diversity + SD, width=.02, colour=Lupin)) +
  ggtitle("Diversity profile (mean diversities)") + 
  scale_color_manual("Lupin", values = c("green", "purple")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=26), axis.text.y.left = element_text(size=26),
        axis.text.x = element_text(angle = 90, vjust=0.5), 
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28)) 
ggMean

# Species accumulation curves
# Lupin
#get richness estimators (for each sample, cumulative)
pool.lupin <- poolaccum(lupin[,5:188])
#plot all: obs richness and  estimators
plot(pool.lupin)

#build the species accumulation curve & rarefaction curve (expected)
lupin.specaccum <- specaccum(lupin[,5:188], method = "rarefaction")
#plot the curve with some predefined settings
plot(lupin.specaccum,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

#build a expected curve (randomization for boxplot comparison)
lupin.specaccum.rand <- specaccum(lupin[,5:188], "random")
#plot both curves ("observed" vs "randomized")
plot(lupin.specaccum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(lupin.specaccum.rand, col="yellow", add=TRUE, pch="+")


# No Lupin
#get richness estimators (for each sample, cumulative)
pool.nolupin <- poolaccum(no.lupin[,5:188])
#plot all: obs richness and  estimators
plot(pool.nolupin)

#build the species accumulation curve & rarefaction curve (expected)
nolupin.specaccum <- specaccum(no.lupin[,5:188], method = "rarefaction")
#plot the curve with some predefined settings
plot(nolupin.specaccum,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

#build a expected curve (randomization for boxplot comparison)
nolupin.specaccum.rand <- specaccum(no.lupin[,5:188], "random")
#plot both curves ("observed" vs "randomized")
plot(nolupin.specaccum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(nolupin.specaccum.rand, col="yellow", add=TRUE, pch="+")


