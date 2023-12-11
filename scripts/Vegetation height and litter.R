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
library(rstatix)
library(lsr)

# Vegetation height ----
# Load data
tidy <- read.csv2(here("data/Final Datasets", "final.tidy.csv"))

# Turn into wide
tidy.res <- tidy %>% 
  filter(Species !=  "Lupinus polyphyllus") %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region =as.factor(Region))  
str(tidy.res)

df <- tidy.res %>% 
  group_by(Region, Site, Lupin) %>% 
  summarise(Veghgt.avg, Litt.avg) %>% 
  distinct()
df

# Histograms
gg.veg <- ggplot(df, aes(x=Veghgt.avg,  fill= Lupin)) +geom_histogram()
gg.veg + ggtitle("Vegetation height")  
  
hist(df$Veghgt.avg, main="Vegetation height")
hist(df$Litt.avg, main="Litter depth")

gg.litter <- ggplot(df, aes(x=Litt.avg,  fill= Lupin)) +geom_histogram()
gg.litter + ggtitle("Litter depth")

# Look for correlations between veg. height and litter
dat <- df[, c(4:5)]
chart.Correlation(dat, histogram=TRUE, pch=19)
cor.test(dat$Veghgt.avg, dat$Litt.avg, method=c("pearson"))


lupin <- df %>% 
  filter(Lupin =="Yes") %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(lupin)

no.lupin <- df %>% 
  filter(Lupin =="No") %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(no.lupin)

dat.lupin <- lupin[, c(4:5)]
chart.Correlation(dat.lupin, histogram=TRUE, pch=19)
cor.test(dat.lupin$Veghgt.avg, dat.lupin$Litt.avg, method=c("pearson"))

dat.nolupin <- no.lupin[, c(4:5)]
chart.Correlation(dat.nolupin, histogram=TRUE, pch=19)
cor.test(dat.nolupin$Veghgt.avg, dat.nolupin$Litt.avg, method=c("pearson"))

# Vegetation height: paired t-test

t.test(lupin$Veghgt.avg, no.lupin$Veghgt.avg, paired = TRUE)

comp.lupin <- list(c("Yes", "No"))

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

gg.vegetation <-ggplot(data=df, aes(x=Lupin, y=Veghgt.avg, fill= Lupin)) + geom_violin() + 
  geom_boxplot(width = 0.1, fill = "white") + xlab("Lupine") +
  ylab("Mean vegetation height (cm)") + scale_fill_manual(values=alpha(c("#008837", "#7b3294"),0.6)) + theme_classic() +
  theme(legend.position = "none")

gg.vegetation2 <- gg.vegetation + stat_compare_means(comparisons = comp.lupin, paired = TRUE, method = "t.test",
                                                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                        symbols = c("****", "***", "**", "*", "ns")))
gg.vegetation2 

# Compare between regions

Uppsala <- df %>% 
  filter(Region == "East region")

Värmland <- df %>% 
  filter(Region == "West region")

t.test(Uppsala$Veghgt.avg, Värmland$Veghgt.avg, paired = TRUE)
comp.region <- list(c("East region" ,"West region"))

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

gg.vegetation <-ggplot(data=df, aes(x=Region, y=Veghgt.avg, fill= Region)) + geom_violin() + 
  geom_boxplot(width = 0.1, fill = "white")+
  ylab("Mean vegetation height (cm)") + scale_fill_manual(values=alpha(c("#d7191c", "#2c7bb6"),0.7)) + theme_classic()

gg.vegetation2 <- gg.vegetation + stat_compare_means(comparisons = comp.region, paired = FALSE, method = "t.test",
                                                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                        symbols = c("****", "***", "**", "*", "ns")))
gg.vegetation2 + ggtitle("All data")   

# Check each region separately
# Check Uppland
Uppsala <- df %>% 
  filter(Region == "East region")

Lupin.upp <- Uppsala %>% 
  filter(Lupin =="Yes")
Lupin.upp

NoLupin.upp <- Uppsala %>% 
  filter(Lupin =="No")
NoLupin.upp

t.test(Lupin.upp$Veghgt.avg, NoLupin.upp$Veghgt.avg, paired = TRUE)

comp.lupin <- list(c("Yes", "No"))

gg.vegh <-ggplot(data=Uppsala, aes(x=Lupin, y=Veghgt.avg, fill= Lupin)) + geom_violin() + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Mean vegetation height (cm)") + scale_fill_manual(values=alpha(c("#008837", "#7b3294"),0.7)) + theme_classic()
gg.vegh2 <- gg.vegh + stat_compare_means(comparisons = comp.lupin, paired = TRUE, method = "t.test", 
                                             symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                symbols = c("****", "***", "**", "*", "ns")))
gg.vegh2 + ggtitle("Uppsala")


# Värmland
Värmland <- df %>% 
  filter(Region == "West region")

Lupin.värm <- Värmland %>% 
  filter(Lupin =="Yes")
Lupin.värm

NoLupin.värm <- Värmland %>% 
  filter(Lupin =="No")
NoLupin.värm

t.test(Lupin.värm$Veghgt.avg, NoLupin.värm$Veghgt.avg, paired = TRUE)

comp.lupin <- list(c("Yes", "No"))

gg.vegh <-ggplot(data=Värmland, aes(x=Lupin, y=Veghgt.avg, fill= Lupin)) + geom_violin() + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Mean vegetation height (cm)") + scale_fill_manual(values=alpha(c("#008837", "#7b3294"),0.7)) + theme_classic()
gg.vegh2 <- gg.vegh + stat_compare_means(comparisons = comp.lupin, paired = TRUE, method = "t.test", 
                                         symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                            symbols = c("****", "***", "**", "*", "ns")))
gg.vegh2 + ggtitle("Värmland")


# Litter depth ---- 
#plot data
plot <- ggplot(df) + geom_point(aes(Veghgt.avg, Litt.avg, color=Lupin)) +
  scale_color_manual(values = c("Yes" = "purple", "No" = "green")) +
  xlab("Mean vegetation height (cm)") + ylab("Mean litter depth (cm)")
plot

lupin <- df %>% 
  filter(Lupin == "Yes")

plot.lupin <- ggplot(lupin) + geom_point(aes(Veghgt.avg, Litt.avg, color=Region)) +
  scale_color_manual(values = c("East region" = "orange", "West region" = "blue")) +
  ggtitle("With Lupin") +
  xlab("Mean vegetation height (cm)") + ylab("Mean litter depth (cm)")
plot.lupin

no.lupin <- df %>% 
  filter(Lupin == "No")

plot.lupin <- ggplot(no.lupin) + geom_point(aes(Veghgt.avg, Litt.avg, color=Region)) +
  scale_color_manual(values = c("East region" = "orange", "West region" = "blue")) +
  ggtitle("Without Lupin") +
  xlab("Mean vegetation height (cm)") + ylab("Mean litter depth (cm)")
plot.lupin


# Litter depth: paired t-test

t.test(lupin$Litt.avg, no.lupin$Litt.avg, paired = TRUE)

comp.lupin <- list(c("Yes", "No"))

gg.litter <-ggplot(data=df, aes(x=Lupin, y=Litt.avg, fill= Lupin)) + geom_violin() + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Mean litter depth (cm)") + xlab("Lupine") +
  scale_fill_manual(values=alpha(c("#008837", "#7b3294"),0.6)) + theme_classic() +
  theme(legend.position = "none")
gg.litter2 <- gg.litter + stat_compare_means(comparisons = comp.lupin, paired = TRUE, method = "t.test", 
                                             symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                symbols = c("****", "***", "**", "*", "ns")))
gg.litter2 

# Differences between regions
comp.region <- list(c("East region" ,"West region"))

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))


t.test(lupin$Litt.avg, no.lupin$Litt.avg, paired = TRUE)
cohensD(lupin$Litt.avg, no.lupin$Litt.avg, method = "paired") # 0.1459 (small effect size)

gg.litter <-ggplot(data=df, aes(x=Region, y=Litt.avg, fill= Region)) + geom_violin() + 
  geom_boxplot(width = 0.1, fill = "white")+
  ylab("Mean litter depth (cm)") + scale_fill_manual(values=alpha(c("#d7191c", "#2c7bb6"),0.7)) + theme_classic() +
  theme(legend.position = "none")
gg.litter2 <- gg.litter + stat_compare_means(comparisons = comp.region, paired = FALSE, method = "t.test",
                                             symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                symbols = c("****", "***", "**", "*", "ns")))
gg.litter2   



# Check each region separately
Uppsala <- df %>% 
  filter(Region == "East region")

Värmland <- df %>% 
  filter(Region == "West region")

t.test(Uppsala$Litt.avg, Värmland$Litt.avg, paired = TRUE)


# Check Uppland
Uppsala <- df %>% 
  filter(Region == "East region")

Lupin.upp <- Uppsala %>% 
  filter(Lupin =="Yes")
Lupin.upp

NoLupin.upp <- Uppsala %>% 
  filter(Lupin =="No")
NoLupin.upp

t.test(Lupin.upp$Litt.avg, NoLupin.upp$Litt.avg, paired = TRUE)

comp.lupin <- list(c("Yes", "No"))

gg.litter <-ggplot(data=Uppsala, aes(x=Lupin, y=Litt.avg, fill= Lupin)) + geom_violin() + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Mean litter depth (cm)") + scale_fill_manual(values=alpha(c("#008837", "#7b3294"),0.7)) + theme_classic()
gg.litter2 <- gg.litter + stat_compare_means(comparisons = comp.lupin, paired = TRUE, method = "t.test", 
                                             symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                symbols = c("****", "***", "**", "*", "ns")))
gg.litter2 + ggtitle("Uppsala")

# Värmland
Värmland <- df %>% 
  filter(Region == "West region")

Lupin.värm <- Värmland %>% 
  filter(Lupin =="Yes")
Lupin.värm

NoLupin.värm <- Värmland %>% 
  filter(Lupin =="No")
NoLupin.värm

t.test(Lupin.värm$Litt.avg, NoLupin.värm$Litt.avg, paired = TRUE)

comp.lupin <- list(c("Yes", "No"))

gg.litter <-ggplot(data=Värmland, aes(x=Lupin, y=Litt.avg, fill= Lupin)) + geom_violin() + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Mean litter depth (cm)") + scale_fill_manual(values=alpha(c("#008837", "#7b3294"),0.7)) + theme_classic()
gg.litter2 <- gg.litter + stat_compare_means(comparisons = comp.lupin, paired = TRUE, method = "t.test", 
                                             symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                symbols = c("****", "***", "**", "*", "ns")))
gg.litter2 + ggtitle("Värmland")

