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
library(ordinal)
library(mclogit)
library(effects)

db <- read.csv2(here("data/Final Datasets", "final.tidy.traits.csv" ))

db.res <- db %>% 
  filter(Species != "Lupinus polyphyllus") %>% 
  filter(Species != "Rubus fruticosus") %>% 
  filter(Species != "Glebionis segetum") %>% 
  #filter(!is.na(Management.dependent)) %>% 
  mutate(County=as.factor(County)) %>% 
  mutate(Lupin=as.factor(Lupin)) %>%  
  #mutate(Management.dependent=as.factor(Management.dependent)) %>% 
  mutate(Competitive =as.factor(Competitive)) %>% 
  mutate(Nutriphilous=as.factor(Nutriphilous)) %>% 
  mutate(sum.comp.nutr=as.factor(sum.comp.nutr)) %>% 
  mutate(sum.reclassified=as.factor(sum.reclassified)) %>% 
  mutate(Seminat.grassland.spp=as.factor(Seminat.grassland.spp))
str(db.res)

traits <- db.res[, c(1,2,4,5,15,22,23,24,25,26)]


#Species list ----
list.occ <- traits %>% 
  group_by(Lupin, Species, sum.reclassified) %>% 
  summarise(Occurrences.site = n())
list.occ

list.abund <- traits %>% 
  group_by(Lupin, Species, sum.reclassified) %>% 
  summarise(abundance = sum(Occurrences))
list.abund

list.all <- left_join(list.occ, list.abund, by=c("Species", "Lupin"))

species.list <- traits %>% 
  group_by(County, Lupin, Species, sum.reclassified) %>% 
  summarise(count = n())
species.list

# Species list: region
list.occ.reg <- traits %>% 
  group_by(County, Lupin, Species, sum.reclassified) %>% 
  summarise(Occurrences.site = n())
list.occ.reg

list.abund.reg <- traits %>% 
  group_by(County, Lupin, Species, sum.reclassified) %>% 
  summarise(abundance = sum(Occurrences))
list.abund.reg

list.county <- left_join(list.occ.reg, list.abund.reg, by=c("Species", "Lupin", "County"))

species.list <- traits %>% 
  group_by(County, Lupin, Species, sum.reclassified) %>% 
  summarise(count = n())
species.list
 
# ------
gg.competition <- traits %>% 
  group_by(Lupin, Competitive) %>% 
  summarise(count = n())
gg.competition

pie1 <- ggplot(gg.competition, aes(x = "", y= count, fill=Competitive)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + geom_label(aes(label = count), position = position_stack(vjust = 0.5), show.legend = FALSE)
pie1 + ggtitle("Occurrences in each competitive category / Lupin treatment")

gg.nutriphilous <- traits %>% 
  group_by(Lupin, Nutriphilous) %>% 
  summarise(count = n())
gg.nutriphilous

pie2 <- ggplot(gg.nutriphilous, aes(x = "", y= count, fill=Nutriphilous)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + geom_label(aes(label = count), position = position_stack(vjust = 0.5), show.legend = FALSE)
pie2 + ggtitle("Occurrences in each nutriphilous category / Lupin treatment")

gg.sum <- traits %>% 
  group_by(Lupin, sum.comp.nutr) %>% 
  summarise(count = n())
gg.sum

pie3 <- ggplot(gg.sum, aes(x = "", y= count, fill=sum.comp.nutr)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + geom_label(aes(label = count), position = position_stack(vjust = 0.5), show.legend = FALSE)
pie3 + ggtitle("Occurrences in sum / Lupin treatment")

gg.sum2 <- traits %>% 
  group_by(Lupin, sum.reclassified) %>% 
  summarise(count = n())
gg.sum2

pie4 <- ggplot(gg.sum2, aes(x = "", y= count, fill=sum.reclassified)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + geom_label(aes(label = count), position = position_stack(vjust = 0.5), show.legend = FALSE)
pie4 + ggtitle("Reclassified sum / Lupin treatment")

gg.sum3 <- traits %>% 
  group_by(Lupin, sum.reclassified) %>% 
  summarise(sum = sum(Occurrences))
gg.sum3


# Species list per category in reclassified sum
species.list <- traits %>% 
  group_by(County, Lupin, Species, sum.reclassified) %>% 
  summarise(count = n())
species.list


# Per county
uppland <- traits %>% filter(County == "Uppsala")

gg.upp <- uppland %>% 
  group_by(Lupin, sum.reclassified) %>% 
  summarise(count = n()) %>% 
  mutate(perc = `count` / sum(`count`)) %>% 
  mutate(labels = scales::percent(perc))
gg.upp

pie.upp <- ggplot(gg.upp, aes(x = "", y= count, fill=sum.reclassified)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + geom_label(aes(label = count), position = position_stack(vjust = 0.5), show.legend = FALSE)
pie.upp + ggtitle("Uppland: Reclassified sum / Lupin treatment")


gg.upp2 <- list.county %>%
  filter(County == "Uppsala") %>% 
  group_by(County, Lupin, sum.reclassified.x) %>% 
  summarise(abund = sum(abundance))
gg.upp2


pie.upp2 <- ggplot(gg.upp2, aes(x = "", y= abund, fill=sum.reclassified.x)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + geom_label(aes(label = abund), position = position_stack(vjust = 0.5), show.legend = FALSE)
pie.upp2 + ggtitle("") + scale_fill_discrete(name = "")


pie <- ggplot(gg.upp2, aes(x = "", y= abund, fill=sum.reclassified.x)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + 
  geom_label(aes(label = abund), position = position_stack(vjust = 0.5), show.legend = FALSE) +
  scale_fill_manual(values=alpha(c("#abdda4", "#80cdc1", "#fee08b", "#f46d43"),0.7)) +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
pie #+ ggtitle("Uppland: Number of occurrences of each management dependent category  / Lupin treatment")

värmland <- traits %>% filter(County == "Värmland")

gg.värm <- värmland %>% 
  group_by(Lupin, sum.reclassified) %>% 
  summarise(count = n()) %>% 
  mutate(perc = `count` / sum(`count`)) %>% 
  mutate(labels = scales::percent(perc))
gg.värm

pie.värm <- ggplot(gg.värm, aes(x = "", y= count, fill=sum.reclassified)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + geom_label(aes(label = count), position = position_stack(vjust = 0.5), show.legend = FALSE)
pie.värm + ggtitle("Värmland: Reclassified sum / Lupin treatment")


pie <- ggplot(gg.värm2, aes(x = "", y= abund, fill=sum.reclassified.x)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + 
  geom_label(aes(label = abund), position = position_stack(vjust = 0.5), show.legend = FALSE) +
  scale_fill_manual(values=alpha(c("#abdda4", "#80cdc1", "#fee08b", "#f46d43"),0.7)) +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
pie #+ ggtitle("Värmland: Number of occurrences of each management dependent category  / Lupin treatment")


gg.värm2 <- list.county %>%
  filter(County == "Värmland") %>% 
  group_by(County, Lupin, sum.reclassified.x) %>% 
  summarise(abund = sum(abundance))
gg.värm2


pie.värm2 <- ggplot(gg.värm2, aes(x = "", y= abund, fill=sum.reclassified.x)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + geom_label(aes(label = abund), position = position_stack(vjust = 0.5), show.legend = FALSE)
pie.värm2 + ggtitle("")+ scale_fill_discrete(name = "")

# Pies with percentage
gg <- management %>% 
  group_by(Lupin, Management.dependent) %>% 
  summarise(count = n()) %>% 
  mutate(perc = `count` / sum(`count`)) 
  #mutate(labels = scales::percent(perc))
gg

pie <- ggplot(gg, aes(x = "", y= perc, fill=Management.dependent)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + 
  geom_label(aes(label = labels), position = position_stack(vjust = 0.5), show.legend = FALSE) +
  scale_fill_manual(values=alpha(c("#abdda4", "#80cdc1", "#fee08b", "#f46d43"),0.7)) +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
pie + labs(fill="Management \ndependency \ncategory")


# Pie all occurrences 
gg <- traits %>% 
  group_by(Lupin, sum.reclassified) %>% 
  summarise(abund = sum(Occurrences))
gg

pie <- ggplot(gg, aes(x = "", y= abund, fill=sum.reclassified)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + 
  geom_label(aes(label = abund), position = position_stack(vjust = 0.5), show.legend = FALSE) +
  scale_fill_manual(values=alpha(c("#abdda4", "#80cdc1", "#fee08b", "#f46d43"),0.7)) +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "top", 
        legend.text=element_text(size=rel(2)), legend.title=element_text(size=rel(2)))
pie + labs(fill="Competitiveness category")


# Model
data.mod <- traits %>% 
  group_by(County, Site, Lupin, sum.comp.nutr) %>% 
  summarise(
    no.species = n(),
    tot.occ = sum(Occurrences)
  )
data.mod

management.mod <- glmmTMB(no.species ~ Management.dependent * Lupin + County + (1|Site) + offset(log(tot.occ)),
                          family=poisson, data=data.mod)
summary(management.mod)
plot(management.mod)
AICc(model.rich.balanced)

visreg(management.mod, "Management.dependent", scale = "response", by="Lupin", rug=FALSE, 
       line = list(col="black"))

#, xlab= "Lupin present", ylab="Species richness / plot"
