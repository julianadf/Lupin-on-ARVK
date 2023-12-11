rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(permute)
library(geosphere)
library(ggplot2)
library(png)
library(grid)
library(cowplot)
library(magick)

# Load the data frames
df <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))

# Create species accumulation curves ----
# Species accumulation plots for each treatment ----
lupin <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  filter(Lupin =="Yes") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(lupin)

no.lupin <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  filter(Lupin =="No") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(no.lupin)

# Lupin
#get richness estimators (for each sample, cumulative)
pool.lupin <- poolaccum(lupin[,5:ncol(lupin)])
#plot all: obs richness and  estimators
plot(pool.lupin)

#build the species accumulation curve & rarefaction curve (expected)
lupin.specaccum <- specaccum(lupin[,5:ncol(lupin)], method = "random")

#build a expected curve (randomization for boxplot comparison)
lupin.specaccum.rand <- specaccum(lupin[,5:ncol(lupin)], "random")
extrapol.lupin <- specpool(lupin[,5:ncol(lupin)])
#plot both curves ("observed" vs "randomized")
plot(lupin.specaccum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main="Lupin")
boxplot(lupin.specaccum.rand, col="yellow", add=TRUE, pch="+")

#renyiaccum
x <- renyiaccum(lupin[,5:ncol(lupin)], scales = c(0, 0.5, 1, 2, 4, Inf), permutations = 100,
                raw = FALSE, collector = FALSE)
## S3 method for class 'renyiaccum'
plot(x, what = c("Collector", "mean", "Qnt 0.025", "Qnt 0.975"), type = "l") # sp. richness=0, looks good
## S3 method for class 'renyiaccum'
persp(x, theta = 220, col = heat.colors(100))


# No Lupin
#get richness estimators (for each sample, cumulative)
pool.nolupin <- poolaccum(no.lupin[,5:ncol(no.lupin)])
#plot all: obs richness and  estimators
plot(pool.nolupin)

#build the species accumulation curve & rarefaction curve (expected)
nolupin.specaccum <- specaccum(no.lupin[,5:ncol(no.lupin)], method = "rarefaction")
#plot the curve with some predefined settings
plot(nolupin.specaccum,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

#build a expected curve (randomization for boxplot comparison)
nolupin.specaccum.rand <- specaccum(no.lupin[,5:ncol(no.lupin)], "random")
#plot both curves ("observed" vs "randomized")
plot(nolupin.specaccum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main= "No Lupin")
boxplot(nolupin.specaccum.rand, col="yellow", add=TRUE, pch="+")

#renyiaccum
y <- renyiaccum(no.lupin[,5:ncol(no.lupin)], scales = c(0, 0.5, 1, 2, 4, Inf), permutations = 100,
                raw = FALSE, collector = FALSE)
## S3 method for class 'renyiaccum'
plot(y, what = c("Collector", "mean", "Qnt 0.025", "Qnt 0.975"), type = "l") # sp. richness=0, looks good
## S3 method for class 'renyiaccum'
persp(y, theta = 220, col = heat.colors(100))


# Plot the accum curves together
lupin.specaccum.rand
nolupin.specaccum.rand

# Extract data to plot with ggplot
gg.lupin <- data.frame(Sites = lupin.specaccum.rand$sites, Richness= lupin.specaccum.rand$richness, 
                       SD= lupin.specaccum.rand$sd)
gg.lupin$Lupin <- "Yes"
gg.nolupin <- data.frame(Sites = nolupin.specaccum.rand$sites, Richness= nolupin.specaccum.rand$richness, 
                         SD= nolupin.specaccum.rand$sd)
gg.nolupin$Lupin <- "No"

# join all data into a long format
gg.all <- bind_rows(gg.lupin, gg.nolupin)

# Plot
gg.treatments <-
  ggplot() +
  #geom_point(data=gg.all, aes(x=Sites, y=Richness, fill= Habitat, group= Habitat)) +
  geom_line(data=gg.all, aes(x=Sites, y=Richness, color=Lupin, group=Lupin), size= 0.9) +
  geom_ribbon(data=gg.all ,aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD), fill=Lupin),alpha=0.2) +
  scale_color_manual("Lupin", values = c("#008837", "#7b3294")) +
  scale_fill_manual("Lupin", values = c("#a6dba0", "#c2a5cf")) +
  theme_minimal() +
  #facet_wrap("Lupin", ncol=3, nrow = 2) +
  theme(legend.position="none", text = element_text(size = 16))
gg.treatments


# Species accumulation curves per region
df <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))

# Divide the data
Uppland <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  filter(County =="Uppsala") %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(Uppland)

Värmland <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  filter(County =="Värmland") %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(Värmland)

# Uppland
#get richness estimators (for each sample, cumulative)
pool.upp <- poolaccum(Uppland[,5:ncol(Uppland)])
#plot all: obs richness and  estimators
plot(pool.upp)
#build the species accumulation curve & rarefaction curve (expected)
Upp.specaccum <- specaccum(Uppland[,5:ncol(Uppland)], method = "rarefaction")

#build a expected curve (randomization for boxplot comparison)
Upp.specaccum.rand <- specaccum(Uppland[,5:ncol(Uppland)], "random")
extrapol.Upp <- specpool(Uppland[,5:ncol(Uppland)])
#plot both curves ("observed" vs "randomized")
plot(Upp.specaccum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main="Uppland")
boxplot(Upp.specaccum.rand, col="yellow", add=TRUE, pch="+")

upp.lupin <- Uppland %>% filter (Lupin == "Yes")
upp.nolupin <- Uppland %>% filter (Lupin == "No")
spec.upp.lup <- specaccum(upp.lupin [,5:ncol(upp.lupin)], "random")
spec.upp.nolup <- specaccum(upp.nolupin [,5:ncol(upp.nolupin)], "random")

gg.lupin1 <- data.frame(Sites = spec.upp.lup$sites, Richness= spec.upp.lup$richness, 
                       SD= spec.upp.lup$sd)
gg.lupin1$Lupin <- "Yes"
gg.nolupin1 <- data.frame(Sites = spec.upp.nolup$sites, Richness= spec.upp.nolup$richness, 
                         SD= spec.upp.nolup$sd)
gg.nolupin1$Lupin <- "No"

# join all data into a long format
gg.all1 <- bind_rows(gg.lupin1, gg.nolupin1)

# Plot
gg.uppland <-
  ggplot() +
  #geom_point(data=gg.all, aes(x=Sites, y=Richness, fill= Habitat, group= Habitat)) +
  geom_line(data=gg.all1, aes(x=Sites, y=Richness, color=Lupin, group=Lupin), size= 0.9) +
  geom_ribbon(data=gg.all1, aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD), fill=Lupin),alpha=0.2) +
  scale_color_manual("Lupin", values = c("#008837", "#7b3294")) +
  scale_fill_manual("Lupin", values = c("#a6dba0", "#c2a5cf")) +
  theme_minimal() + 
  #facet_wrap("Lupin", ncol=3, nrow = 2) +
  theme(legend.position="none", text = element_text(size = 16))
gg.uppland 


# Värmland
#get richness estimators (for each sample, cumulative)
pool.värm <- poolaccum(Värmland[,5:ncol(Värmland)])
#plot all: obs richness and  estimators
plot(pool.värm)
#build the species accumulation curve & rarefaction curve (expected)
värm.specaccum <- specaccum(Värmland[,5:ncol(Värmland)], method = "rarefaction")

#build a expected curve (randomization for boxplot comparison)
värm.specaccum.rand <- specaccum(Värmland[,5:ncol(Värmland)], "random")
extrapol.värm <- specpool(Värmland[,5:ncol(Värmland)])
#plot both curves ("observed" vs "randomized")
plot(värm.specaccum, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", main="Värmland")
boxplot(värm.specaccum.rand, col="yellow", add=TRUE, pch="+")

värm.lupin <- Värmland %>% filter (Lupin == "Yes")
värm.nolupin <- Värmland %>% filter (Lupin == "No")
spec.värm.lup <- specaccum(värm.lupin [,5:ncol(värm.lupin)], "random")
spec.värm.nolup <- specaccum(värm.nolupin [,5:ncol(värm.nolupin)], "random")

gg.lupin2 <- data.frame(Sites = spec.värm.lup$sites, Richness= spec.värm.lup$richness, 
                        SD= spec.värm.lup$sd)
gg.lupin2$Lupin <- "Yes"
gg.nolupin2 <- data.frame(Sites = spec.värm.nolup$sites, Richness= spec.värm.nolup$richness, 
                          SD= spec.värm.nolup$sd)
gg.nolupin2$Lupin <- "No"

# join all data into a long format
gg.all2 <- bind_rows(gg.lupin2, gg.nolupin2)

# Plot
gg.värmland <-
  ggplot() +
  #geom_point(data=gg.all, aes(x=Sites, y=Richness, fill= Habitat, group= Habitat)) +
  geom_line(data=gg.all2, aes(x=Sites, y=Richness, color=Lupin, group=Lupin), size= 0.9) +
  geom_ribbon(data=gg.all2, aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD), fill=Lupin),alpha=0.2) +
  scale_color_manual("Lupin", values = c("#008837", "#7b3294")) +
  scale_fill_manual("Lupin", values = c("#a6dba0", "#c2a5cf")) +
  theme_minimal() + 
  #facet_wrap("Lupin", ncol=3, nrow = 2) +
  theme(legend.position=c(1.3,0.5), text = element_text(size = 16))
gg.värmland 

# Specaccum for all data
df.all <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(df.all)

specaccum.rand <- specaccum(df.all[,5:ncol(df.all)], "random")

# Extract data to plot with ggplot
gg.data <- data.frame(Sites = specaccum.rand$sites, Richness= specaccum.rand$richness, 
                       SD= specaccum.rand$sd)


# Plot
gg <-
  ggplot() +
  #geom_point(data=gg.all, aes(x=Sites, y=Richness, fill= Habitat, group= Habitat)) +
  geom_line(data=gg.data, aes(x=Sites, y=Richness), size= 0.9) +
  geom_ribbon(data=gg.data ,aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2) +
  # scale_color_manual("Lupin", values = c("#008837", "#7b3294")) +
  # scale_fill_manual("Lupin", values = c("#a6dba0", "#c2a5cf")) +
  theme_minimal() + xlab("Plots") +
  #facet_wrap("Lupin", ncol=3, nrow = 2) +
  theme(legend.position="none", text = element_text(size = 16))
gg

# Plot together ----

figureS2 <- ggarrange(gg.treatments, gg.uppland, gg.värmland,
                      labels = c("a", "b", "c"),
                      ncol = 2, nrow = 2)
figureS2

# ----

# Diversity profiles ----
# Divide the dataset into lupin / no lupin
lupin <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  filter(Lupin =="Yes") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(lupin)

no.lupin <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  filter(Lupin =="No") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(no.lupin)


# Calculate Hill diversities and plot them

# Hill diversities for plots with lupin
Hill.lupin <- renyi(lupin[,5:ncol(lupin)], hill = TRUE)
plot(Hill.lupin)
Hill.lupin$Lupin <- "Yes"

# Hill diversities for plots without lupin
Hill.nolupin <- renyi(no.lupin[,5:ncol(no.lupin)], hill = TRUE)
plot(Hill.nolupin)
Hill.nolupin$Lupin <- "No"

# Join the two results
Hill.data <- bind_rows(Hill.lupin, Hill.nolupin)

# Turn into tidy data 
Hill <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11, na.rm=TRUE) %>%
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
  scale_color_manual("Lupin", values = c("#008837", "#7b3294"))
gg.hill

# Plot alpha values for species richness (0), Shannon (1), Simpson (2), Berger-Parker (inf) for mean category
ggMean <- ggplot(HillMean, aes(x=q, y=Diversity, group = Lupin)) + 
  geom_line(linetype="dashed") + 
  geom_point(aes(colour = Lupin), size=4) + 
  geom_pointrange(aes(ymin=Diversity-SD, ymax=Diversity + SD, colour=Lupin)) +
  geom_errorbar(aes(ymin=Diversity - SD, ymax= Diversity + SD, width=.02, colour=Lupin)) +
  #ggtitle("Diversity profile (mean diversities)") + 
  #xlab("Order q") +
  scale_color_manual("Lupine", values = c("#008837", "#7b3294")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = c(0.9,0.9),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=22), axis.text.y.left = element_text(size=22),
        axis.title.x = element_blank(),
        axis.text.x = element_text( vjust=0.5), 
        #axis.title.x = element_text(size=28),
        axis.title.y = element_text(size=28), 
        legend.title = element_text(size=18), legend.text = element_text(size=18)) 

ggMean
# Diversity profile: region ----

# Divide the dataset into Uppland /Värmland
Uppland <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  filter(Region == "Uppsala")
str(Uppland)

Värmland <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  filter(Region == "Värmland")
str(Värmland)

# Calculate Hill diversities and plot them

# Hill diversities for plots with lupin
Hill.Upp <- renyi(Uppland[,5:ncol(Uppland)], hill = TRUE)
plot(Hill.Upp)
Hill.Upp$Region <- "Uppland"

# Hill diversities for plots without lupin
Hill.Värm <- renyi(Värmland[,5:ncol(Värmland)], hill = TRUE)
plot(Hill.Värm)
Hill.Värm$Region <- "Värmland"

# Join the two results
Hill.data <- bind_rows(Hill.Upp, Hill.Värm)

# Turn into tidy data 
Hill <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11, na.rm=TRUE) %>%
  mutate(Region = as.factor(Region))
Hill

# Calculate a mean per plot type and select alpha 0,1,2,inf
HillMean <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11,na.rm=TRUE) %>%
  mutate(Region = as.factor(Region)) %>%
  group_by(., Region, q) %>% 
  summarise(Diversity=mean(Diversity)) %>%  
  filter(q %in% c("0", "1", "2", "Inf"))
HillMean

# Calculate standard error per plot type and select alpha 0,1,2,inf
HillSD <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11,na.rm=TRUE) %>%
  mutate(Region = as.factor(Region)) %>%
  group_by(., Region, q) %>%
  summarise(SD=sd(Diversity)/sqrt(length(Diversity))) %>% 
  filter(q %in% c("0", "1", "2", "Inf"))
HillSD

HillMean$SD <- HillSD$SD

# Diversity plot for all alpha values
gg.hill <- ggplot(Hill, aes(x=q, y=Diversity))+ geom_point(aes(colour = Region))  + 
  ggtitle("Hill diversity") + theme_bw() +
  xlab("Order q") +
  scale_color_manual("Region", values = c("#d7191c", "#2c7bb6"))
gg.hill

# Plot alpha values for species richness (0), Shannon (1), Simpson (2), Berger-Parker (inf) for mean category
ggMean.reg <- ggplot(HillMean, aes(x=q, y=Diversity, group = Region)) + 
  geom_line(linetype="dashed") + 
  geom_point(aes(colour = Region), size=4) + 
  geom_pointrange(aes(ymin=Diversity-SD, ymax=Diversity + SD, colour=Region)) +
  geom_errorbar(aes(ymin=Diversity - SD, ymax= Diversity + SD, width=.02, colour=Region)) +
  #ggtitle("Diversity profile (mean diversities)") + 
  #xlab("Order q") +
  scale_color_manual("Region", values = c("#d7191c", "#2c7bb6"), labels=c("East", "West")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.9, 0.9),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=22), axis.text.y.left = element_text(size=22),
        axis.title = element_blank(),
        #axis.text.x = element_text(angle = 90, vjust=0.5), 
        #axis.title.y = element_text(size=28), axis.title.x = element_text(size=28),
        legend.title = element_text(size=18), legend.text = element_text(size=18)) 

ggMean.reg

# Diversity profiles for each region separately ----

# Uppland - east
lupin.upp <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  filter(Lupin =="Yes") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  filter(Region =="Uppsala")
str(lupin.upp)

no.lupin.upp <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  filter(Lupin =="No") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  filter(Region =="Uppsala")
str(no.lupin.upp)


# Calculate Hill diversities and plot them

# Hill diversities for plots with lupin
Hill.lupin <- renyi(lupin.upp[,5:ncol(lupin.upp)], hill = TRUE)
plot(Hill.lupin)
Hill.lupin$Lupin <- "Yes"

# Hill diversities for plots without lupin
Hill.nolupin <- renyi(no.lupin.upp[,5:ncol(no.lupin.upp)], hill = TRUE)
plot(Hill.nolupin)
Hill.nolupin$Lupin <- "No"

# Join the two results
Hill.data <- bind_rows(Hill.lupin, Hill.nolupin)

# Turn into tidy data 
Hill <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11, na.rm=TRUE) %>%
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
  scale_color_manual("Lupin", values = c("#008837", "#7b3294"))
gg.hill

# Plot alpha values for species richness (0), Shannon (1), Simpson (2), Berger-Parker (inf) for mean category
ggMean <- ggplot(HillMean, aes(x=q, y=Diversity, group = Lupin)) + 
  geom_line(linetype="dashed") + 
  geom_point(aes(colour = Lupin), size=4) + 
  geom_pointrange(aes(ymin=Diversity-SD, ymax=Diversity + SD, colour=Lupin)) +
  geom_errorbar(aes(ymin=Diversity - SD, ymax= Diversity + SD, width=.02, colour=Lupin)) +
  #ggtitle("Diversity profile (mean diversities)") + 
  xlab("Order q") +
  scale_color_manual("Lupine", values = c("#008837", "#7b3294")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.9,0.9),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=22), axis.text.y.left = element_text(size=22),
        #axis.title = element_blank(),
        axis.text.x = element_text(vjust=0.5), 
        axis.title.y = element_blank(), axis.title.x = element_text(size=28),
        legend.title = element_text(size=18), legend.text = element_text(size=18)) 

ggMean 

# Värmland
lupin.värm <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  filter(Lupin =="Yes") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  filter(Region =="Värmland")
str(lupin.värm)

no.lupin.värm <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  filter(Lupin =="No") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  filter(Region =="Värmland")
str(no.lupin.värm)


# Calculate Hill diversities and plot them

# Hill diversities for plots with lupin
Hill.lupin <- renyi(lupin.värm[,5:ncol(lupin.värm)], hill = TRUE)
plot(Hill.lupin)
Hill.lupin$Lupin <- "Yes"

# Hill diversities for plots without lupin
Hill.nolupin <- renyi(no.lupin.värm[,5:ncol(no.lupin.värm)], hill = TRUE)
plot(Hill.nolupin)
Hill.nolupin$Lupin <- "No"

# Join the two results
Hill.data <- bind_rows(Hill.lupin, Hill.nolupin)

# Turn into tidy data 
Hill <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11, na.rm=TRUE) %>%
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
  scale_color_manual("Lupin", values = c("#008837", "#7b3294"))
gg.hill

# Plot alpha values for species richness (0), Shannon (1), Simpson (2), Berger-Parker (inf) for mean category
ggMean <- ggplot(HillMean, aes(x=q, y=Diversity, group = Lupin)) + 
  geom_line(linetype="dashed") + 
  geom_point(aes(colour = Lupin), size=4) + 
  geom_pointrange(aes(ymin=Diversity-SD, ymax=Diversity + SD, colour=Lupin)) +
  geom_errorbar(aes(ymin=Diversity - SD, ymax= Diversity + SD, width=.02, colour=Lupin)) +
  #ggtitle("Diversity profile (mean diversities)") + 
  xlab("Order q") +
  scale_color_manual("Lupine", values = c("#008837", "#7b3294")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.9,0.9),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=22), axis.text.y.left = element_text(size=22),
        #axis.title = element_blank(),
        axis.text.x = element_text(vjust=0.5), 
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28),
        legend.title = element_text(size=18), legend.text = element_text(size=18)) 

ggMean 


# Unbalanced -----
df.unb <- read.csv2(here("data/Final Datasets", "df.unbalanced.csv"))
# Diversity profiles ----
# Divide the dataset into lupin / no lupin
lupin <- df.unb %>% 
  #dplyr::select(-Lupinus.polyphyllus) %>% #already took it away
  filter(Lupin =="Yes") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(lupin)

no.lupin <- df.unb %>% 
  #dplyr::select(-Lupinus.polyphyllus) %>% #already took it away
  filter(Lupin =="No") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(no.lupin)


# Calculate Hill diversities and plot them

# Hill diversities for plots with lupin
Hill.lupin <- renyi(lupin[,5:ncol(lupin)], hill = TRUE)
plot(Hill.lupin)
Hill.lupin$Lupin <- "Yes"

# Hill diversities for plots without lupin
Hill.nolupin <- renyi(no.lupin[,5:ncol(no.lupin)], hill = TRUE)
plot(Hill.nolupin)
Hill.nolupin$Lupin <- "No"

# Join the two results
Hill.data <- bind_rows(Hill.lupin, Hill.nolupin)

# Turn into tidy data 
Hill <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11, na.rm=TRUE) %>%
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
  scale_color_manual("Lupin", values = c("#008837", "#7b3294"))
gg.hill

# Plot alpha values for species richness (0), Shannon (1), Simpson (2), Berger-Parker (inf) for mean category
ggMean <- ggplot(HillMean, aes(x=q, y=Diversity, group = Lupin)) + 
  geom_line(linetype="dashed") + 
  geom_point(aes(colour = Lupin), size=4) + 
  geom_pointrange(aes(ymin=Diversity-SD, ymax=Diversity + SD, colour=Lupin)) +
  geom_errorbar(aes(ymin=Diversity - SD, ymax= Diversity + SD, width=.02, colour=Lupin)) +
  #ggtitle("Diversity profile (mean diversities)") + 
  xlab("Order q") +
  scale_color_manual("Lupin", values = c("#008837", "#7b3294")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=22), axis.text.y.left = element_text(size=22),
        axis.text.x = element_text(angle = 90, vjust=0.5), 
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28),
        legend.title = element_text(size=18), legend.text = element_text(size=18)) 

ggMean

# Diversity profile: region ----

# Divide the dataset into Uppland /Värmland
Uppland <- df.unb %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  filter(Region == "Uppsala")
str(Uppland)

Värmland <- df.unb %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  filter(Region == "Värmland")
str(Värmland)


# Calculate Hill diversities and plot them

# Hill diversities for plots with lupin
Hill.Upp <- renyi(Uppland[,5:ncol(Uppland)], hill = TRUE)
plot(Hill.Upp)
Hill.Upp$Region <- "Uppland"

# Hill diversities for plots without lupin
Hill.Värm <- renyi(Värmland[,5:ncol(Värmland)], hill = TRUE)
plot(Hill.Värm)
Hill.Värm$Region <- "Värmland"

# Join the two results
Hill.data <- bind_rows(Hill.Upp, Hill.Värm)

# Turn into tidy data 
Hill <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11, na.rm=TRUE) %>%
  mutate(Region = as.factor(Region))
Hill

# Calculate a mean per plot type and select alpha 0,1,2,inf
HillMean <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11,na.rm=TRUE) %>%
  mutate(Region = as.factor(Region)) %>%
  group_by(., Region, q) %>% 
  summarise(Diversity=mean(Diversity)) %>%  
  filter(q %in% c("0", "1", "2", "Inf"))
HillMean

# Calculate standard error per plot type and select alpha 0,1,2,inf
HillSD <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11,na.rm=TRUE) %>%
  mutate(Region = as.factor(Region)) %>%
  group_by(., Region, q) %>%
  summarise(SD=sd(Diversity)/sqrt(length(Diversity))) %>% 
  filter(q %in% c("0", "1", "2", "Inf"))
HillSD

HillMean$SD <- HillSD$SD

# Diversity plot for all alpha values
gg.hill <- ggplot(Hill, aes(x=q, y=Diversity))+ geom_point(aes(colour = Region))  + 
  ggtitle("Hill diversity") + theme_bw() +
  xlab("Order q") +
  scale_color_manual("Region", values = c("#d7191c", "#2c7bb6"))
gg.hill

# Plot alpha values for species richness (0), Shannon (1), Simpson (2), Berger-Parker (inf) for mean category
ggMean <- ggplot(HillMean, aes(x=q, y=Diversity, group = Region)) + 
  geom_line(linetype="dashed") + 
  geom_point(aes(colour = Region), size=4) + 
  geom_pointrange(aes(ymin=Diversity-SD, ymax=Diversity + SD, colour=Region)) +
  geom_errorbar(aes(ymin=Diversity - SD, ymax= Diversity + SD, width=.02, colour=Region)) +
  #ggtitle("Diversity profile (mean diversities)") + 
  xlab("Order q") +
  scale_color_manual("Region", values = c("#d7191c", "#2c7bb6")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=22), axis.text.y.left = element_text(size=22),
        axis.text.x = element_text(angle = 90, vjust=0.5), 
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28),
        legend.title = element_text(size=18), legend.text = element_text(size=18)) 

ggMean

# Diversity profiles for each region separately ----

# Uppland
lupin.upp <- df.unb %>% 
  filter(Lupin =="Yes") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  filter(Region =="Uppsala")
str(lupin.upp)

no.lupin.upp <- df.unb %>% 
  filter(Lupin =="No") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  filter(Region =="Uppsala")
str(no.lupin.upp)


# Calculate Hill diversities and plot them

# Hill diversities for plots with lupin
Hill.lupin <- renyi(lupin.upp[,5:ncol(lupin.upp)], hill = TRUE)
plot(Hill.lupin)
Hill.lupin$Lupin <- "Yes"

# Hill diversities for plots without lupin
Hill.nolupin <- renyi(no.lupin.upp[,5:ncol(no.lupin.upp)], hill = TRUE)
plot(Hill.nolupin)
Hill.nolupin$Lupin <- "No"

# Join the two results
Hill.data <- bind_rows(Hill.lupin, Hill.nolupin)

# Turn into tidy data 
Hill <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11, na.rm=TRUE) %>%
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
  scale_color_manual("Lupin", values = c("#008837", "#7b3294"))
gg.hill

# Plot alpha values for species richness (0), Shannon (1), Simpson (2), Berger-Parker (inf) for mean category
ggMean <- ggplot(HillMean, aes(x=q, y=Diversity, group = Lupin)) + 
  geom_line(linetype="dashed") + 
  geom_point(aes(colour = Lupin), size=4) + 
  geom_pointrange(aes(ymin=Diversity-SD, ymax=Diversity + SD, colour=Lupin)) +
  geom_errorbar(aes(ymin=Diversity - SD, ymax= Diversity + SD, width=.02, colour=Lupin)) +
  #ggtitle("Diversity profile (mean diversities)") + 
  xlab("Order q") +
  scale_color_manual("Lupin", values = c("#008837", "#7b3294")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=22), axis.text.y.left = element_text(size=22),
        axis.text.x = element_text(angle = 90, vjust=0.5), 
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28),
        legend.title = element_text(size=18), legend.text = element_text(size=18)) 

ggMean + ggtitle("Uppland")

# Värmland
lupin.värm <- df.unb %>% 
  filter(Lupin =="Yes") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  filter(Region =="Värmland")
str(lupin.värm)

no.lupin.värm <- df.unb %>% 
  filter(Lupin =="No") %>% 
  rename(Region = County) %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  filter(Region =="Värmland")
str(no.lupin.värm)


# Calculate Hill diversities and plot them

# Hill diversities for plots with lupin
Hill.lupin <- renyi(lupin.värm[,5:ncol(lupin.värm)], hill = TRUE)
plot(Hill.lupin)
Hill.lupin$Lupin <- "Yes"

# Hill diversities for plots without lupin
Hill.nolupin <- renyi(no.lupin.värm[,5:ncol(no.lupin.värm)], hill = TRUE)
plot(Hill.nolupin)
Hill.nolupin$Lupin <- "No"

# Join the two results
Hill.data <- bind_rows(Hill.lupin, Hill.nolupin)

# Turn into tidy data 
Hill <- Hill.data %>% 
  gather(., key="q", value="Diversity", 1:11, na.rm=TRUE) %>%
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
  scale_color_manual("Lupin", values = c("#008837", "#7b3294"))
gg.hill

# Plot alpha values for species richness (0), Shannon (1), Simpson (2), Berger-Parker (inf) for mean category
ggMean <- ggplot(HillMean, aes(x=q, y=Diversity, group = Lupin)) + 
  geom_line(linetype="dashed") + 
  geom_point(aes(colour = Lupin), size=4) + 
  geom_pointrange(aes(ymin=Diversity-SD, ymax=Diversity + SD, colour=Lupin)) +
  geom_errorbar(aes(ymin=Diversity - SD, ymax= Diversity + SD, width=.02, colour=Lupin)) +
  #ggtitle("Diversity profile (mean diversities)") + 
  xlab("Order q") +
  scale_color_manual("Lupin", values = c("#008837", "#7b3294")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=22), axis.text.y.left = element_text(size=22),
        axis.text.x = element_text(angle = 90, vjust=0.5), 
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28),
        legend.title = element_text(size=18), legend.text = element_text(size=18)) 

ggMean + ggtitle("Värmland")
