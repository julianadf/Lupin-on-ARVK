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
library(MuMIn)

# Load the created species richness data set (unbalanced, main text) ---- 
data.1 <- read.csv2(here("data/Final Datasets", "species.richness.noLP.csv" ))

# Restructure
data.res <- data.1 %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region = as.factor(Region)) %>% 
  mutate(obs = 1:n())
str(data.res)

# Plot data
boxplot <- ggplot(data.res, aes(x=Region, y=rich, color=Lupin)) + geom_boxplot()
boxplot

# Run model
model.rich <- glmmTMB(rich ~ Lupin * Region , family= nbinom1, data=data.res)
summary(model.rich)
car::Anova(model.rich, type= "III")
AICc(model.rich)
mod_dharma1 <- model.rich %>% simulateResiduals(n=1000)
plot(mod_dharma1)
testDispersion(mod_dharma1)
testDispersion(mod_dharma1, alternative = "less")
plot1 <- visreg(model.rich, "Lupin", scale="response", by= "Region", rug=FALSE, line = list(col="black"), 
       xlab= "Lupine present", ylab="Species richness / plot", gg=TRUE) + theme_bw() +
  theme(text = element_text(size = 16))
# to switch the order of the panels
plot1 + facet_grid(~fct_rev(Region))#, levels = c("West region", "East region"))

# Model without interaction
model.rich2 <- glmmTMB(rich ~ Lupin + Region + (1|Site), family= "poisson", data= data.res)
summary(model.rich2)
car::Anova(model.rich, type= "III")
AICc(model.rich2)

# Null model
model.rich3 <- glmmTMB(rich ~ 1 + (1|Site), family= "poisson", data=data.res)
AICc(model.rich3)

# Load the created species richness data set (balanced, supplementary material) ----
balanced.data <- read.csv2(here("data/Final Datasets", "balanced.dataset.csv" ))

boxplot.balanced <- ggplot(balanced.data, aes(x=County, y=rich, color=Lupin)) + geom_boxplot()
boxplot.balanced


# Restructure
balanced.data <- balanced.data %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  mutate(County=as.factor(County))
str(balanced.data)

# Run model
model.rich.balanced <- glmmTMB(rich ~ Lupin * County + (1|Site), family= "poisson", data=balanced.data)
summary(model.rich.balanced)
car::Anova(model.rich.balanced, type= "III")
AICc(model.rich.balanced)
mod_dharma1 <- model.rich.balanced %>% simulateResiduals(n=1000)
plot(mod_dharma1)
testDispersion(mod_dharma1)
testDispersion(mod_dharma1, alternative = "less")
visreg(model.rich.balanced, "Lupin", scale="response", by= "County", rug=FALSE, 
       line = list(col="black"), xlab= "Lupin present", ylab="Species richness / plot")

# Model without interaction
model.rich.balanced2 <- glmmTMB(rich ~ Lupin + County + (1|Site), family= "poisson", data=balanced.data)
AICc(model.rich.balanced2)

# Null model
model.rich.balanced3 <- glmmTMB(rich ~ 1 + (1|Site), family= "poisson", data=balanced.data)
AICc(model.rich.balanced3)


# Run model for Uppsala only (12 sites)
uppsala <- balanced.data %>% 
  filter(County == "Uppsala")

model.uppsala <- glmmTMB(rich ~ Lupin + (1|Site), family= "poisson", data=uppsala)
summary(model.uppsala)
car::Anova(model.uppsala, type= "III")
mod_dharma1 <- model.uppsala %>% simulateResiduals(n=1000)
plot(mod_dharma1)
testDispersion(mod_dharma1)
testDispersion(mod_dharma1, alternative = "less")
visreg(model.uppsala, "Lupin", scale="response", rug=FALSE, 
       line = list(col="black"), xlab= "Lupin present", ylab="Species richness / plot")

# Load neighbouring habitat file 
neighbour <- read.csv2(here("data/Final Datasets", "neighbor.habitat.csv" ))

# Join with the database
data.neigh <- left_join(data.res, neighbour, by=c("Region", "Site", "Lupin"))
data.neigh <- data.neigh %>% 
  mutate(Open = as.factor(Open)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Region = as.factor(Region)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(data.neigh)

# Run model
model.neigh <- glmmTMB(rich ~ Lupin * Open + (1|Site), family= "poisson", data=data.neigh)
summary(model.neigh)
car::vif(model.neigh, type="terms")
car::Anova(model.neigh, type= "III")
AICc(model.neigh)
mod_dharma1 <-model.neigh %>% simulateResiduals(n=1000)
plot(mod_dharma1)
testDispersion(mod_dharma1)
testDispersion(mod_dharma1, alternative = "less")
visreg(model.neigh, "Lupin", scale="response", by= "Open", rug=FALSE, 
       line = list(col="black"), xlab= "Lupin present", ylab="Species richness / plot") 
       


Tommy.plot <- ggplot(data.neigh, aes(x=fct_inorder(Site), y=rich, color= Open)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank())
  #ggtitle("New sites & no trees & low frequency LP in no.Lupin plots")
Tommy.plot


# Test to remove regions from the analysis
model.rich.test <- glmmTMB(rich ~ Lupin + (1|Site), family= "poisson", data=balanced.data)
summary(model.rich.balanced)
car::Anova(model.rich.balanced, type= "III")
AICc(model.rich.balanced)
mod_dharma1 <- model.rich.balanced %>% simulateResiduals(n=1000)
plot(mod_dharma1)
testDispersion(mod_dharma1)
testDispersion(mod_dharma1, alternative = "less")
visreg(model.rich.balanced, "Lupin", scale="response", by= "County", rug=FALSE, 
       line = list(col="black"), xlab= "Lupin present", ylab="Species richness / plot")

