rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)
library(BiodiversityR)

# Load data
df <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))

df.res <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  mutate(County=as.factor(County))  
str(df.res)

comm <- df.res[,5:ncol(df.res)]
env <- df.res[,c(1,4)]

# All data ----
rank.abun.No <- rankabundance(comm, y=env,  factor="Lupin", level = "No")
rank.abun.Yes <- rankabundance(comm, y=env,  factor="Lupin", level = "Yes")

rankabunplot(rank.abun.No, scale='abundance', addit=FALSE, specnames=c(1:5), col="green")
rankabunplot(rank.abun.Yes, scale='abundance', addit=TRUE, specnames=c(1:5), col="purple")


RA.data <- rankabuncomp(comm, y=env, factor='Lupin', scale='abundance', legend=TRUE, specnames=c(1:3))

BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12, family="Arial"),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())
            
plotgg1 <- ggplot(data=RA.data, aes(x = rank, y = abundance)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=1) +
  geom_point(aes(colour=Grouping, shape=Grouping), size=5, alpha=0.7) +
  geom_text_repel(data=subset(RA.data, labelit == TRUE), 
                  aes(colour=Grouping, label=species), 
                  angle=0, nudge_x=13, nudge_y=1, show.legend=FALSE) +
  BioR.theme +
  scale_color_manual("Lupin", values = c("#008837", "#7b3294"))  +
  labs(x = "rank", y = "abundance", colour = "Lupin", shape = "Lupin")

plotgg1

# Uppland ----
df.upp <- df.res %>% filter(County=="East region")
comm.upp <- df.upp[,5:ncol(df.res)]
env.upp <- df.upp[,c(1,4)]
rank.abun.No <- rankabundance(comm.upp, y=env.upp,  factor="Lupin", level = "No")
rank.abun.Yes <- rankabundance(comm.upp, y=env.upp,  factor="Lupin", level = "Yes")

rankabunplot(rank.abun.No, scale='abundance', addit=FALSE, specnames=c(1:5), col="green")
rankabunplot(rank.abun.Yes, scale='abundance', addit=TRUE, specnames=c(1:5), col="purple")


RA.data <- rankabuncomp(comm.upp, y=env.upp, factor='Lupin', scale='abundance', legend=TRUE, specnames=c(1:3))

BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12, family="Arial"),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

plotgg1 <- ggplot(data=RA.data, aes(x = rank, y = abundance)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=1) +
  geom_point(aes(colour=Grouping, shape=Grouping), size=5, alpha=0.7) +
  geom_text_repel(data=subset(RA.data, labelit == TRUE), 
                  aes(colour=Grouping, label=species), 
                  angle=0, nudge_x=13, nudge_y=1, show.legend=FALSE) +
  BioR.theme +
  scale_color_manual("Lupin", values = c("#008837", "#7b3294"))  +
  labs(x = "rank", y = "abundance", colour = "Lupin", shape = "Lupin")

plotgg1 + ggtitle("Uppland")

# Värmland ----
df.värm <- df.res %>% filter(County=="West region")
comm.värm <- df.värm[,5:ncol(df.res)]
env.värm <- df.värm[,c(1,4)]
rank.abun.No <- rankabundance(comm.värm, y=env.värm,  factor="Lupin", level = "No")
rank.abun.Yes <- rankabundance(comm.värm, y=env.värm,  factor="Lupin", level = "Yes")

rankabunplot(rank.abun.No, scale='abundance', addit=FALSE, specnames=c(1:5), col="green")
rankabunplot(rank.abun.Yes, scale='abundance', addit=TRUE, specnames=c(1:5), col="purple")


RA.data <- rankabuncomp(comm.värm, y=env.värm, factor='Lupin', scale='abundance', legend=TRUE, specnames=c(1:3))

BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12, family="Arial"),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

plotgg1 <- ggplot(data=RA.data, aes(x = rank, y = abundance)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=1) +
  geom_point(aes(colour=Grouping, shape=Grouping), size=5, alpha=0.7) +
  geom_text_repel(data=subset(RA.data, labelit == TRUE), 
                  aes(colour=Grouping, label=species), 
                  angle=0, nudge_x=13, nudge_y=1, show.legend=FALSE) +
  BioR.theme +
  scale_color_manual("Lupin", values = c("#008837", "#7b3294"))  +
  labs(x = "rank", y = "abundance", colour = "Lupin", shape = "Lupin")

plotgg1 + ggtitle("Värmland")
