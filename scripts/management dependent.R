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


db <- read.csv2(here("data/Final Datasets", "final.tidy.csv" ))

db.res <- db %>% 
  filter(Species != "Lupinus.polyphyllus") %>% 
  filter(!is.na(Management.dependent)) %>% 
  mutate(County=as.factor(County)) %>% 
  mutate(Lupin=as.factor(Lupin)) %>%  
  mutate(Management.dependent=as.factor(Management.dependent))
str(db.res)

management <- db.res[, c(1,2,4,5,14,15)]

gg <- management %>% 
  group_by(Lupin, Management.dependent) %>% 
  summarise(count = n())
gg

pie <- ggplot(gg, aes(x = "", y= count, fill=Management.dependent)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + geom_label(aes(label = count), position = position_stack(vjust = 0.5), show.legend = FALSE)
pie + ggtitle("Number of occurrences of each management dependent category  / Lupin treatment")

# Per county
uppland <- management %>% filter(County == "Uppsala")

gg.upp <- uppland %>% 
  group_by(Lupin, Management.dependent) %>% 
  summarise(count = n()) %>% 
  mutate(perc = `count` / sum(`count`)) %>% 
  mutate(labels = scales::percent(perc))
gg.upp

pie <- ggplot(gg.upp, aes(x = "", y= perc, fill=Management.dependent)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + 
  geom_label(aes(label = labels), position = position_stack(vjust = 0.5), show.legend = FALSE) +
  scale_fill_manual(values=alpha(c("#f46d43", "#fee08b", "#a6d96a", "#1a9850"),0.7)) +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
pie #+ ggtitle("Uppland: Number of occurrences of each management dependent category  / Lupin treatment")

värmland <- management %>% filter(County == "Värmland")

gg.värm <- värmland %>% 
  group_by(Lupin, Management.dependent) %>% 
  summarise(count = n()) %>% 
  mutate(perc = `count` / sum(`count`)) %>% 
  mutate(labels = scales::percent(perc))
gg.värm

pie <- ggplot(gg.värm, aes(x = "", y= perc, fill=Management.dependent)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + 
  geom_label(aes(label = labels), position = position_stack(vjust = 0.5), show.legend = FALSE) +
  scale_fill_manual(values=alpha(c("#f46d43", "#fee08b", "#a6d96a", "#1a9850"),0.7)) +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
pie #+ ggtitle("Värmland: Number of occurrences of each management dependent category  / Lupin treatment")

# Pies with percentage
gg <- management %>% 
  group_by(Lupin, Management.dependent) %>% 
  summarise(count = n()) %>% 
  mutate(perc = `count` / sum(`count`)) %>% 
  mutate(labels = scales::percent(perc))
gg

pie <- ggplot(gg, aes(x = "", y= perc, fill=Management.dependent)) + geom_col(color="black") + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + 
  geom_label(aes(label = labels), position = position_stack(vjust = 0.5), show.legend = FALSE) +
  scale_fill_manual(values=alpha(c("#f46d43", "#fee08b", "#a6d96a", "#1a9850"),0.7)) +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
pie + labs(fill="Management \ndependency \ncategory")


# Model
data.mod <- management %>% 
  group_by(County, Site, Lupin, Management.dependent) %>% 
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
