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

extrapol.lupin <- specpool(selected.pairs[,5:188])

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
  gather(., key="q", value="Diversity", 1:11,na.rm=TRUE) %>%
  mutate(Lupin = as.factor(Lupin))
Hill

# Calculate a mean per plot type and select alpha 0,1,2,inf
HillMean <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11,na.rm=TRUE) %>%
  mutate(Lupin = as.factor(Lupin)) %>%
  group_by(., Lupin, q) %>% 
  summarise(Diversity=mean(Diversity)) %>%  
  filter(q %in% c("0", "1", "2", "Inf"))
HillMean

# Calculate standard error per plot type and select alpha 0,1,2,inf
HillSD <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11,na.rm=TRUE) %>%
  mutate(Lupin = as.factor(Lupin)) %>%
  group_by(., Lupin, q) %>%
  summarise(SD=sd(Diversity)/sqrt(length(Diversity))) %>% 
  filter(q %in% c("0", "1", "2", "Inf"))
HillSD

HillMean$SD <- HillSD$SD

# Diversity plot for all alpha values
gg.hill <- ggplot(Hill, aes(x=q, y=Diversity))+ geom_point(aes(colour = Lupin))  + 
  ggtitle("Hill diversity") + theme_bw() +
  xlab("Order q") +
  scale_color_manual("Lupin", values = c("#f1a340", "#998ec3"))
gg.hill

# Plot alpha values for species richness (0), Shannon (1), Simpson (2), Berger-Parker (inf) for mean category
ggMean <- ggplot(HillMean, aes(x=q, y=Diversity, group = Lupin)) + 
  geom_line(linetype="dashed") + 
  geom_point(aes(colour = Lupin), size=4) + 
  geom_pointrange(aes(ymin=Diversity-SD, ymax=Diversity + SD, colour=Lupin)) +
  geom_errorbar(aes(ymin=Diversity - SD, ymax= Diversity + SD, width=.02, colour=Lupin)) +
  #ggtitle("Diversity profile (mean diversities)") + 
  xlab("Order q") +
  scale_color_manual("Lupin", values = c("#f1a340", "#998ec3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=22), axis.text.y.left = element_text(size=22),
        axis.text.x = element_text(angle = 90, vjust=0.5), 
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28),
        legend.title = element_text(size=18), legend.text = element_text(size=18)) 

ggMean

# Species accumulation curves
# Lupin
#get richness estimators (for each sample, cumulative)
pool.lupin <- poolaccum(lupin[,5:188])
#plot all: obs richness and  estimators
plot(pool.lupin)

#build the species accumulation curve & rarefaction curve (expected)
lupin.specaccum <- specaccum(lupin[,5:188], method = "rarefaction")

#build a expected curve (randomization for boxplot comparison)
lupin.specaccum.rand <- specaccum(lupin[,5:188], "random")
extrapol.lupin <- specpool(lupin[,5:188])
#plot both curves ("observed" vs "randomized")
plot(lupin.specaccum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(lupin.specaccum.rand, col="yellow", add=TRUE, pch="+")

#renyiaccum
x <- renyiaccum(lupin[,5:188], scales = c(0, 0.5, 1, 2, 4, Inf), permutations = 100,
                raw = FALSE, collector = FALSE)
## S3 method for class 'renyiaccum'
plot(x, what = c("Collector", "mean", "Qnt 0.025", "Qnt 0.975"), type = "l") # sp. richness=0, looks good
## S3 method for class 'renyiaccum'
persp(x, theta = 220, col = heat.colors(100))


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

#renyiaccum
y <- renyiaccum(no.lupin[,5:188], scales = c(0, 0.5, 1, 2, 4, Inf), permutations = 100,
                raw = FALSE, collector = FALSE)
## S3 method for class 'renyiaccum'
plot(y, what = c("Collector", "mean", "Qnt 0.025", "Qnt 0.975"), type = "l") # sp. richness=0, looks good
## S3 method for class 'renyiaccum'
persp(y, theta = 220, col = heat.colors(100))


# Plot the accum curves together
lupin.specaccum.rand
nolupin.specaccum.rand
plot(nolupin.specaccum, ci.type="poly", col="#e66101", lwd=1, ci.lty=0, ci.col="#fdb863")
plot(lupin.specaccum, ci.type="poly", col="#5e3c99", lwd=1, ci.lty=0, ci.col="#b2abd2", add=TRUE)
legend(title= "Lupin", x="bottomright", legend=c("Yes", "No"), fill = c("#b2abd2", "#fdb863"), box.lty=0)

# Check frequency of lupins in lupin-plots ----
freq <- lupin %>% 
  mutate(lupin.freq = `Lupinus polyphyllus`/16)
freq

lupin.freq<- freq[,c(1:2,189)] %>% mutate(County = dplyr::recode(County, "Örebro"= "Värmland"))
hist(lupin.freq$lupin.freq)

uppsala <- lupin.freq %>% filter(County =="Uppsala")
värmland <- lupin.freq %>% filter(County !="Uppsala") %>% mutate(County = dplyr::recode(County, "Örebro"= "Värmland"))

t.test(uppsala$lupin.freq, värmland$lupin.freq)

gg.county <-ggplot(data=lupin.freq, aes(x=County, y=lupin.freq)) + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Medelfrekvens lupin / yta")
gg.county

