rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(glmmTMB)
library(DHARMa)
library(ggplot2)
library(ggpubr)
library(png)
library(grid)
library(cowplot)
library(magick)
library(PerformanceAnalytics)

raw.data <- read.csv2(here("data", "Raw data.env.csv" ))

# Extract the columns I want
raw.data.1 <- raw.data[,c(2:3, 6, 8, 10, 11, 16, 33)]


data <- raw.data %>% 
  mutate(County = dplyr::recode(County, "Örebro"= "Värmland")) %>%
  mutate(County = as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(Parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  group_by(County, Site, parnr, Lupin, Veghgt.avg, Litt.avg) %>%
  summarise()
data

hist(data$Veghgt.avg)
hist(data$Litt.avg)

# Look for correlations between veg. height and litter
dat <- data[, c(5:6)]
chart.Correlation(dat, histogram=TRUE, pch=19)

dat.lupin <- lupin[, c(5:6)]
chart.Correlation(dat.lupin, histogram=TRUE, pch=19)

dat.nolupin <- no.lupin[, c(5:6)]
chart.Correlation(dat.nolupin, histogram=TRUE, pch=19)

# Vegetation height: paired t-test

t.test(lupin$Veghgt.avg, no.lupin$Veghgt.avg, paired = TRUE)

comp.lupin <- list(c("Yes", "No"))

gg.vegetation <-ggplot(data=data, aes( x=Lupin, y=Veghgt.avg, fill= Lupin)) + geom_violin() + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Mean vegetation height (cm)") + scale_fill_manual(values=c("green", "purple")) + theme_classic()
gg.vegetation2 <- gg.lupin + stat_compare_means(comparisons = comp.lupin, paired = TRUE, method = "t.test")
gg.vegetation2

vegetation <- ggplot(data) + geom_dotplot(aes(Veghgt.avg))


#plot data
plot <- ggplot(data) + geom_point(aes(Veghgt.avg, Litt.avg, color=Lupin)) +
  scale_color_manual(values = c("Yes" = "purple", "No" = "green")) +
  xlab("Mean vegetation height (cm)") + ylab("Mean litter depth (cm)")
plot

lupin <- data %>% 
  filter(Lupin == "Yes")

plot.lupin <- ggplot(lupin) + geom_point(aes(Veghgt.avg, Litt.avg, color=County)) +
  scale_color_manual(values = c("Uppsala" = "orange", "Värmland" = "blue")) +
  ggtitle("With Lupin") +
  xlab("Mean vegetation height (cm)") + ylab("Mean litter depth (cm)")
plot.lupin

no.lupin <- data %>% 
  filter(Lupin == "No")

plot.lupin <- ggplot(no.lupin) + geom_point(aes(Veghgt.avg, Litt.avg, color=County)) +
  scale_color_manual(values = c("Uppsala" = "orange", "Värmland" = "blue")) +
  ggtitle("Without Lupin") +
  xlab("Mean vegetation height (cm)") + ylab("Mean litter depth (cm)")
plot.lupin


# Litter depth: paired t-test

t.test(lupin$Litt.avg, no.lupin$Litt.avg, paired = TRUE)

comp.lupin <- list(c("Yes", "No"))

gg.litter <-ggplot(data=data, aes(x=Lupin, y=Litt.avg, fill= Lupin)) + geom_violin() + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Mean litter depth (cm)") + scale_fill_manual(values=c("green", "purple")) + theme_classic()
gg.litter2 <- gg.litter + stat_compare_means(comparisons = comp.lupin, paired = TRUE, method = "t.test")
gg.litter2


# Run analysis ----
model.veghgt <- glmmTMB(Veghgt.avg  ~ Lupin * County + (1|Site), family= "Gamma", data=data) 
summary(model.veghgt)
mod_dharma1 <- model.veghgt %>% simulateResiduals(n=1000)
plotResiduals(model.veghgt, rank = TRUE, quantreg = FALSE)
plot(mod_dharma1)
visreg(model.veghgt, "Lupin", scale="response", by= "County", rug=FALSE, line = list(col="black")) 

model.litt <- glmmTMB(log10(Litt.avgLOG)  ~ Lupin * County + (1|Site), family= "", data=data) 
summary(model.litt)
mod_dharma1 <- model.litt %>% simulateResiduals(n=1000)
plotResiduals(model.litt, rank = TRUE, quantreg = FALSE)
plot(mod_dharma1)
visreg(model.litt, "Lupin", scale="response", by= "County", rug=FALSE, line = list(col="black")) 

