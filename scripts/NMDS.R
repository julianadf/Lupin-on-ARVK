# Multivariate analysis on plant species at landscape level

library(tidyverse)
library(vegan)

rm(list=ls())
setwd("//storage.slu.se/Home$/jada0002/My Documents/GINFRA/R/Analysis habitat")

# Butterflies ----------------
bf.tidy <- read.csv2("Multivariate/data/bf.tidy.csv")
# Create a data frame
bf.df <- bf.tidy %>% 
  group_by(Landscape, Transect_type) %>% 
  spread(., key="Species", value = "Indiv")
bf.df

# Load dataframe
bf.data <- read.csv2("Multivariate/data/butterfly.dataframe.csv")
Habitat_type <- bf.df$Transect_type

NMDS.bf <- metaMDS(bf.data, k=2, try = 100, trymax = 1000)
NMDS.bf
stressplot(NMDS.bf) #Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions. 

# Create figure with ggplot
bf.scores <- as.data.frame(scores(NMDS.bf))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
bf.scores$landscape <- rownames(bf.scores)  # create a column of site names, from the rownames of data.scores
bf.scores$habitat <- Habitat_type  #  add the grp variable created earlier
head(bf.scores)  #look at the data

bfsp.scores <- as.data.frame(scores(NMDS.bf, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
bfsp.scores$species <- rownames(bfsp.scores)  # create a column of species, from the rownames of species.scores
head(bfsp.scores)  #look at the data

bf.BF <- bf.scores[bf.scores$category == "Between fields", ][chull(bf.scores[bf.scores$category == 
                                                                              "Between fields", c("NMDS1", "NMDS2")]), ]  # hull values 
bf.BR <- bf.scores[bf.scores$category == "Big road", ][chull(bf.scores[bf.scores$category == 
                                                                                  "Big road", c("NMDS1", "NMDS2")]), ]  # hull 
bf.P <- bf.scores[bf.scores$category == "Pasture", ][chull(bf.scores[bf.scores$category == 
                                                                              "Pasture", c("NMDS1", "NMDS2")]), ]  # hull values 
bf.PL <- bf.scores[bf.scores$category == "Powerline", ][chull(bf.scores[bf.scores$category == 
                                                                                  "Powerline", c("NMDS1", "NMDS2")]), ]  # hull values
bf.SR <- bf.scores[bf.scores$category == "Small road", ][chull(bf.scores[bf.scores$category == 
                                                                       "Small road", c("NMDS1", "NMDS2")]), ]  # hull values 

hull.bf <- rbind(bf.BF, bf.BR, bf.P, bf.PL, bf.SR)  #combine habitats
hull.bf

gg.bf<- ggplot() + 
  geom_polygon(data=hull.bf,aes(x=NMDS1,y=NMDS2,fill=habitat,group=habitat),alpha=0.25) + # add the convex hulls
  #geom_text(data=bfsp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=bf.scores,aes(x=NMDS1,y=NMDS2,shape=habitat,colour=habitat),size=3) + # add the point markers
  geom_text(data=bf.scores,aes(x=NMDS1,y=NMDS2,label=habitat),size=1,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("Between fields" = "yellow", "Big road" = "grey", "Pasture" = "green", "Powerline" = "blue", "Small road" = "red")) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
gg.bf


# Bumblebees ----------------
bb.df <- read.csv2("Multivariate/data/bumblebee.dataframe.csv")
category <- bb.df$Category
bb.df <- bb.df[,3:21]

NMDS.bb <- metaMDS(bb.df, k=2, try = 100, trymax = 1000)
NMDS.bb
stressplot(NMDS.bb) #Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions. 

# Create figure with ggplot
bb.scores <- as.data.frame(scores(NMDS.bb))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
bb.scores$landscape <- rownames(bb.scores)  # create a column of site names, from the rownames of data.scores
bb.scores$category <- category  #  add the grp variable created earlier
head(bb.scores)  #look at the data

bbsp.scores <- as.data.frame(scores(NMDS.bb, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
bbsp.scores$species <- rownames(bbsp.scores)  # create a column of species, from the rownames of species.scores
head(bbsp.scores)  #look at the data

bb.PL.LRD <- bb.scores[bb.scores$category == "PL.LRD", ][chull(bb.scores[bb.scores$category == 
                                                                              "PL.LRD", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
bb.NoPL.HRD <- bb.scores[bb.scores$category == "NoPL.HRD", ][chull(bb.scores[bb.scores$category == 
                                                                                  "NoPL.HRD", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
bb.PL.HRD <- bb.scores[bb.scores$category == "PL.HRD", ][chull(bb.scores[bb.scores$category == 
                                                                              "PL.HRD", c("NMDS1", "NMDS2")]), ]  # hull values for grp c
bb.NoPL.LRD <- bb.scores[bb.scores$category == "NoPL.LRD", ][chull(bb.scores[bb.scores$category == 
                                                                                  "NoPL.LRD", c("NMDS1", "NMDS2")]), ]  # hull values for grp d

hull.bb <- rbind(bb.PL.LRD, bb.PL.HRD, bb.NoPL.LRD, bb.NoPL.HRD)  #combine grp.a and grp.b
hull.bb

gg.bb<- ggplot() + 
  geom_polygon(data=hull.bb,aes(x=NMDS1,y=NMDS2,fill=category,group=category),alpha=0.30) + # add the convex hulls
  geom_text(data=bbsp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=bb.scores,aes(x=NMDS1,y=NMDS2,shape=category,colour=category),size=3) + # add the point markers
  geom_text(data=bb.scores,aes(x=NMDS1,y=NMDS2,label=landscape),size=1,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("PL.LRD" = "green", "NoPL.HRD" = "yellow", "PL.HRD" = "blue", "NoPL.LRD" = "red")) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
gg.bb

# Plants ---------------

plant.df <- read.csv2("Multivariate/data/plant.dataframe.csv")
category <- plant.df$Category
plant.df <- plant.df[,3:136]


NMDS.pp <- metaMDS(plant.df, k=2, try = 100, trymax = 1000)
NMDS.pp
stressplot(NMDS.pp) #Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions. 

# Create figure with ggplot
pp.scores <- as.data.frame(scores(NMDS.pp))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
pp.scores$landscape <- rownames(pp.scores)  # create a column of site names, from the rownames of data.scores
pp.scores$category <- category  #  add the grp variable created earlier
head(pp.scores)  #look at the data

ppsp.scores <- as.data.frame(scores(NMDS.pp, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
ppsp.scores$species <- rownames(ppsp.scores)  # create a column of species, from the rownames of species.scores
head(ppsp.scores)  #look at the data

pp.PL.LRD <- pp.scores[pp.scores$category == "PL.LRD", ][chull(pp.scores[pp.scores$category == 
                                                                              "PL.LRD", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
pp.NoPL.HRD <- pp.scores[pp.scores$category == "NoPL.HRD", ][chull(pp.scores[pp.scores$category == 
                                                                                  "NoPL.HRD", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
pp.PL.HRD <- pp.scores[pp.scores$category == "PL.HRD", ][chull(pp.scores[pp.scores$category == 
                                                                              "PL.HRD", c("NMDS1", "NMDS2")]), ]  # hull values for grp c
pp.NoPL.LRD <- pp.scores[pp.scores$category == "NoPL.LRD", ][chull(pp.scores[pp.scores$category == 
                                                                                  "NoPL.LRD", c("NMDS1", "NMDS2")]), ]  # hull values for grp d

hull.pp <- rbind(pp.PL.LRD, pp.PL.HRD, pp.NoPL.LRD, pp.NoPL.HRD)  #combine grp.a and grp.b
hull.pp

gg.pp<- ggplot() + 
  geom_polygon(data=hull.pp,aes(x=NMDS1,y=NMDS2,fill=category,group=category),alpha=0.30) + # add the convex hulls
  #geom_text(data=ppsp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=pp.scores,aes(x=NMDS1,y=NMDS2,shape=category,colour=category),size=3) + # add the point markers
  geom_text(data=pp.scores,aes(x=NMDS1,y=NMDS2,label=landscape),size=1,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("PL.LRD" = "green", "NoPL.HRD" = "yellow", "PL.HRD" = "blue", "NoPL.LRD" = "red")) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
gg.pp


# Indicator species ####
indicator.df <- read.csv2("Multivariate/data/indicator.dataframe.csv")
category <- indicator.df$Category
indicator.df <- indicator.df[,3:29]

NMDS.pi <- metaMDS(indicator.df, k=2, try = 100, trymax = 1000)
NMDS.pi
stressplot(NMDS.pi) #Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions. 

# Create figure with ggplot
pi.scores <- as.data.frame(scores(NMDS.pi))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
pi.scores$landscape <- rownames(pi.scores)  # create a column of site names, from the rownames of data.scores
pi.scores$category <- category  #  add the grp variable created earlier
head(pi.scores)  #look at the data

pisp.scores <- as.data.frame(scores(NMDS.pi, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
pisp.scores$species <- rownames(pisp.scores)  # create a column of species, from the rownames of species.scores
head(pisp.scores)  #look at the data

pi.PL.LRD <- pi.scores[pi.scores$category == "PL.LRD", ][chull(pi.scores[pi.scores$category == 
                                                                           "PL.LRD", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
pi.NoPL.HRD <- pi.scores[pi.scores$category == "NoPL.HRD", ][chull(pi.scores[pi.scores$category == 
                                                                               "NoPL.HRD", c("NMDS1", "NMDS2")]), ]  # hull values for grp B
pi.PL.HRD <- pi.scores[pi.scores$category == "PL.HRD", ][chull(pi.scores[pi.scores$category == 
                                                                           "PL.HRD", c("NMDS1", "NMDS2")]), ]  # hull values for grp c
pi.NoPL.LRD <- pi.scores[pi.scores$category == "NoPL.LRD", ][chull(pi.scores[pi.scores$category == 
                                                                               "NoPL.LRD", c("NMDS1", "NMDS2")]), ]  # hull values for grp d

hull.pi <- rbind(pi.PL.LRD, pi.PL.HRD, pi.NoPL.LRD, pi.NoPL.HRD)  #combine grp.a and grp.b
hull.pi

gg.pi<- ggplot() + 
  geom_polygon(data=hull.pi,aes(x=NMDS1,y=NMDS2,fill=category,group=category),alpha=0.30) + # add the convex hulls
  geom_text(data=pisp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=pi.scores,aes(x=NMDS1,y=NMDS2,shape=category,colour=category),size=3) + # add the point markers
  geom_text(data=pi.scores,aes(x=NMDS1,y=NMDS2,label=landscape),size=1,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("PL.LRD" = "green", "NoPL.HRD" = "yellow", "PL.HRD" = "blue", "NoPL.LRD" = "red")) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
gg.pi

