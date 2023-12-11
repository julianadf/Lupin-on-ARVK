rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(lattice)
library(ggpubr)
library(png)
library(cowplot)
library(ordinal)


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
  mutate(Seminat.grassland.spp=as.factor(Seminat.grassland.spp)) %>% 
  mutate(sum.reclassified = dplyr::recode(sum.reclassified, "0"= "Very poor competitor")) %>% 
  mutate(sum.reclassified = dplyr::recode(sum.reclassified, "1"= "Poor competitor")) %>% 
  mutate(sum.reclassified = dplyr::recode(sum.reclassified, "2"= "Moderate competitor")) %>% 
  mutate(sum.reclassified = dplyr::recode(sum.reclassified, "3"= "Strong competitor"))
str(db.res)

traits <- db.res[, c(1,2,4,5,15,22,23,24,25,26)]

# Pies 
# all
gg <- traits %>% 
  group_by(Lupin, sum.reclassified) %>% 
  summarise(abund = sum(Occurrences))
gg

pie.all <- ggplot(gg, aes(x = "", y= abund, fill=sum.reclassified)) + geom_col(color="black") +# + coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + 
  geom_label(aes(label = abund), position = position_stack(vjust = 0.5), show.legend = FALSE) +
  scale_fill_manual(values=alpha(c("#abdda4", "#80cdc1", "#fee08b", "#f46d43"),0.7)) +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "", 
        axis.text.y = element_text(size = 12), strip.text.x = element_text(size = 12),
        legend.text=element_text(size=rel(3)), legend.title=element_text(size=rel(3)))+
  labs(fill="")
pie.all 

legend <-  cowplot::get_legend(pie.all)
as_ggplot(legend) + labs(fill="")

# Per county
uppland <- traits %>% filter(County == "Uppsala")

gg.upp <- uppland %>% 
  filter(County == "Uppsala") %>% 
  group_by(County, Lupin, sum.reclassified) %>% 
  summarise(abund = sum(Occurrences))
gg.upp

pie.upp <- ggplot(gg.upp, aes(x = "", y= abund, fill=sum.reclassified)) + geom_col(color="black") + #coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + 
  geom_label(aes(label = abund), position = position_stack(vjust = 0.5), show.legend = FALSE) +
  scale_fill_manual(values=alpha(c("#abdda4", "#80cdc1", "#fee08b", "#f46d43"),0.7)) +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 12), strip.text.x = element_text(size = 12),
        legend.position = "none")
pie.upp 

värmland <- traits %>% filter(County == "Värmland")

gg.värm <- värmland %>% 
  group_by(Lupin, sum.reclassified) %>% 
  group_by(County, Lupin, sum.reclassified) %>% 
  summarise(abund = sum(Occurrences))
gg.värm

pie.värm <- ggplot(gg.värm, aes(x = "", y= abund, fill=sum.reclassified)) + geom_col(color="black") + #coord_polar(theta = "y") +
  facet_wrap(vars(Lupin), nrow =1) + 
  geom_label(aes(label = abund), position = position_stack(vjust = 0.5), show.legend = FALSE) +
  scale_fill_manual(values=alpha(c("#abdda4", "#80cdc1", "#fee08b", "#f46d43"),0.7)) +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 12), strip.text.x = element_text(size = 12),
        legend.position = "none")
pie.värm 

