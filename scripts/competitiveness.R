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
library(emmeans)
library(bbmle)

db <- read.csv2(here("data/Final Datasets", "final.tidy.traits.csv" ))

db.res <- db %>% 
  filter(Species != "Lupinus polyphyllus") %>% 
  filter(Species != "Rubus fruticosus") %>% 
  filter(Species != "Glebionis segetum") %>% 
  #filter(sum.reclassified != 0) %>% 
  #filter(!is.na(Management.dependent)) %>% 
  mutate(County=as.factor(County)) %>% 
  mutate(Lupin=as.factor(Lupin)) %>%  
  #mutate(Management.dependent=as.factor(Management.dependent)) %>% 
  mutate(Competitive =as.factor(Competitive)) %>% 
  mutate(Nutriphilous=as.factor(Nutriphilous)) %>% 
  mutate(sum.comp.nutr=as.factor(sum.comp.nutr)) %>% 
  mutate(sum.reclassified=as.factor(sum.reclassified)) %>% 
  mutate(Seminat.grassland.spp=as.factor(Seminat.grassland.spp)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region"))
str(db.res)

traits <- db.res[, c(1,2,4,5,15,25)]

data.mod <- traits %>% 
  group_by(County, Site, Lupin, sum.reclassified) %>% 
  summarise(
    no.species = n(),
    tot.occ = sum(Occurrences)
  )
data.mod


model <- mblogit(sum.reclassified ~ County * Lupin, weights = tot.occ, random = ~1|Site, data=data.mod)
summary(model)
model2 <- mblogit(sum.reclassified ~ County + Lupin, weights = tot.occ, random = ~1|Site, data=data.mod)
summary(model2)
#Delta AIC
bbmle::AICtab(model,model2, weights=TRUE)
# AIC
L <- list(model,model2)
bbmle::AICtab(L)

x <- emmeans(model, ~sum.reclassified|County + Lupin)
plot(x)
# p-values:
y <- emmeans(model, pairwise~County + Lupin|sum.reclassified)
summary(y)
plot(y)
z <- emmeans(model, pairwise~Lupin|County + sum.reclassified)
summary(z)
plot(z)


dodge <- position_dodge(width=.8)
(ggplot(model)
  +facet_wrap(~County)
  +geom_bar(
    aes(fill=response,
        x=Infl,
        y=pred),
    stat='identity',position=dodge,width=.8)
  +geom_errorbar(
    aes(x=Infl,
        ymin=lower,
        ymax=upper,group=response),
    position=dodge,width=.4))


# Merge group zero and group 1 and analyse -----
traitsrereclass <- traits
traitsrereclass$sum.reclassified[traitsrereclass$sum.reclassified == 0] <- 1
traitsrereclass

data.modre <- traitsrereclass %>% 
  group_by(County, Site, Lupin, sum.reclassified) %>% 
  summarise(
    no.species = n(),
    tot.occ = sum(Occurrences)
  )
data.modre


model.2 <- mblogit(sum.reclassified ~  Lupin * County, weights = tot.occ, random = ~1|Site, data=data.modre)
summary(model.2)
AIC(model.2)

model2.2 <- mblogit(sum.reclassified ~ Lupin + County, weights = tot.occ, random = ~1|Site, data=data.modre)
summary(model2.2)
AIC(model2.2)

model2.3 <- mblogit(sum.reclassified ~ 1, weights = tot.occ, random = ~1|Site, data=data.modre)
summary(model2.3)
AIC(model2.3)

#Delta AIC
bbmle::AICtab(model.2, model2.2, weights=TRUE)
# AIC
L <- list(model.2,model2.2)
bbmle::AICtab(L)


x2 <- emmeans(model.2, ~sum.reclassified|County + Lupin)
plot(x2)
# p-values:
y2 <- emmeans(model.2, pairwise~County + Lupin|sum.reclassified)
summary(y2)
plot(y2, )
z2 <- emmeans(model.2, pairwise~Lupin|County + sum.reclassified)
summary(z2)
plot(z2)

w2 <- emmeans(model.2, pairwise~ Lupin|sum.reclassified + County )
summary(w2)

# Make nice plot
# Compare least-square means of each habitat type
em <- emmeans(model.2, pairwise~County + Lupin|sum.reclassified)
emms <- emmeans(model.2, "Lupin")
pairs(emms)
plot(emms, comparisons=TRUE)

em.lup <- as.data.frame(em)
# create database for ggplot
gg.em <- em.lup %>% 
  filter(contrast == ".") %>% 
  unite("treat", County:Lupin)
gg.em

ggplot(gg.em, aes(treat, prob)) +
  facet_wrap(vars(sum.reclassified)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.1) +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank())


em <- multcomp::cld(y2, Letters = letters)
em <- as.data.frame(em)
em.2 <- em %>% 
  unite("treat", County:Lupin, remove = FALSE, sep = " ") %>% 
  mutate(sum.reclassified = dplyr::recode(sum.reclassified, "1"= "Poor competitors")) %>% 
  mutate(sum.reclassified = dplyr::recode(sum.reclassified, "2"= "Moderate competitors")) %>% 
  mutate(sum.reclassified = dplyr::recode(sum.reclassified, "3"= "Strong competitors")) 
em.2

# Change the order of the rows 
em.3 <- em.2[order(factor(em.2$treat, levels = c("West region No", "West region Yes","East region No", "East region Yes"), ordered = T)),]

colours <- em.2$color[order(em.2$Lupin)]

Lupin <- ggplot(em.3, aes(treat, prob)) +
  facet_wrap(vars(sum.reclassified)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL)) + #, width = 0.1) +
  #coord_cartesian(ylim = c(0,21)) +
  geom_text(aes(y = 0.8, label = .group)) +
  #xlab("Habitat type") +
  #ylab("Number of butterfly species")+
  theme_bw() +
  scale_x_discrete(limits=unique(em.3$treat)) +
  ylab("P(occurrence)")+
    theme(panel.grid.major.y = element_blank(), 
        axis.text.x = element_text(angle = 90, color = colours), 
        axis.title.x = element_blank(), 
        text = element_text(size=16))
Lupin 



# For presentations
Lupin <- ggplot(em.2, aes(treat, prob, col=County)) +
  facet_wrap(vars(sum.reclassified)) +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), size= 1) + #, width = 0.1) +
  scale_color_manual(values= c("#d7191c", "#2c7bb6"), labels = c("Uppland", "Värmland"))+ 
  #coord_cartesian(ylim = c(0,21)) +
  geom_text(aes(y = 0.8, label = .group)) +
  #xlab("Habitat type") +
  #ylab("Number of butterfly species")+
  theme_bw() +
  ylab("P(occurrence)") +
  theme(panel.grid.major.y = element_blank(), 
        legend.title= element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), 
        text = element_text(size=20))
Lupin
