rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(permute)
library(geosphere)

# NMDS with subset of all data----
raw.data <- read.csv2(here("data", "Raw data.csv" ))

data.1 <- raw.data %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(Parnr)) %>% 
  mutate(Lupin.NoLupin = as.factor(Lupin.NoLupin)) %>% 
  group_by(County, Site, parnr, Lupin.NoLupin, Species) %>%
  summarise(Occ= sum(Occurrence))
data.1

# use the randomly pre-selected pairs to select data
# Load the selected pairs
selected.pairs <- read.csv2(here("data", "selected pairs per site.csv"))

selected.pairs <- selected.pairs %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr))
str(selected.pairs)

# Use the selected pairs to subset the entire dataset
final.selection <- semi_join(data.1, selected.pairs, by= c("Site", "parnr"))

data <- final.selection %>% 
  mutate(County = dplyr::recode(County, "Örebro"= "Värmland")) %>%
  rename(Lupin = Lupin.NoLupin) %>%
  mutate(Lupin = case_when(Lupin == "Lupin" ~ "Yes", TRUE ~ "No")) %>% 
  mutate(County = as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) 
str(data)

df <- data %>% 
  group_by(County, Site, Lupin) %>% 
  spread(., key="Species", value = "Occ") %>% 
  replace(is.na(.), 0)
df 

plants <- sqrt(df[,5:188])#square root transform to reduce the relative influence of the most frequent species, which otherwise will tend to dominate the dissimilarity matrix
lupinus <- as.data.frame(df[,4])
NMDS <- metaMDS(df[,5:188], distance="bray", k=2, try = 100, trymax = 1000)
NMDS
stressplot(NMDS) #Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions. 

##pull points from NMDS
NMDS1 <- NMDS$points[,1] ##also found using scores(NMDS)
NMDS2 <- NMDS$points[,2]
envfit(NMDS ~ Lupin, data=df)
envfit(NMDS ~ County, data=df)

plant.plot <-cbind(plants, NMDS1, NMDS2, lupinus)

# quick plot ordination
p <- ggplot(plant.plot, aes(NMDS1, NMDS2, color=Lupin)) +
  geom_point(position=position_jitter(.1), shape=3) + ##separates overlapping points
  stat_ellipse(type='t',size =1)+ ##draws 95% confidence interval ellipses
  theme_minimal()
p
# complicated plot with labels
# Create figure with ggplot
plant.scores <- as.data.frame(scores(NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
#plant.scores$Site <- rownames(plant.scores)   # create a column of site names, from the rownames of data.scores
plant.scores$Lupin <- lupinus  #  add the grp variable created earlier
head(plant.scores)  #look at the data

sp.scores <- as.data.frame(scores(NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
sp.scores$species <- rownames(sp.scores)  # create a column of species, from the rownames of species.scores
head(sp.scores)  #look at the data

Lupin <- plant.scores[plant.scores$Lupin == "Yes", ][chull(plant.scores[plant.scores$Lupin == 
                                                                            "Yes", c("NMDS1", "NMDS2")]), ]  # hull values 
NoLupin <- plant.scores[plant.scores$Lupin == "No", ][chull(plant.scores[plant.scores$Lupin == 
                                                                                "No", c("NMDS1", "NMDS2")]), ]  # hull 


hull.plants <- rbind(Lupin, NoLupin)  #combine plots
hull.plants

gg.lupin<- ggplot() + 
  geom_polygon(data=hull.plants,aes(x=NMDS1,y=NMDS2,fill=Lupin , group=Lupin),alpha=0.25) + # add the convex hulls
  geom_text(data=sp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=4) +  # add the species labels
  geom_point(data=plant.scores,aes(x=NMDS1,y=NMDS2,shape=Lupin,colour=Lupin),size=2) + # add the point markers
  #geom_text(data=plant.scores,aes(x=NMDS1,y=NMDS2,label=Lupin),size=1,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("Lupin" = "#fdae61", "NoLupin" = "#2c7bb6")) +
  scale_fill_manual(values=c("Lupin" = "#fdae61", "NoLupin" = "#2c7bb6")) +
  scale_shape_discrete("Lupin") +
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
gg.lupin

# Mantel test for spatial autocorrelation ----
# Load coordinate data
spatial.data <- read.csv2(here("data", "sites.spatial.csv"))

#remove duplicates
spatial <- spatial.data %>% 
  group_by(Site, lat, lng) %>% 
  summarise()
spatial # I have to manually remove duplicates for some sites

spatial.data <- read.csv2(here("data", "1coordset.csv"))

geo <- data.frame(as.numeric(spatial.data$lng), as.numeric(spatial.data$lat))

# calcluate distance matrix
d.geo = distm(geo, fun = distHaversine)
dist.geo = as.dist(d.geo)

#abundance data frame - bray curtis dissimilarity
plant.dist <- vegdist(plants, method='bray') #abundance based

#abundance vs geographic 
abund_geo  = mantel(plant.dist, dist.geo, method = "pearson", permutations = 9999, na.rm = TRUE)
abund_geo

# PERMANOVA ----
# NUll hypothesis: Groups do not differ in spread or position in multivariate space.
plants.transformed <- sqrt(df[,5:188])#square root transform# square root transform
plants <- cbind(df[,1:4], plants.transformed)
plant.dist <- vegdist(plants[,5:188], method='bray') #abundance based
env <- plants[,c(1:2, 4)]
perm <- how(within = Within(type = "free"),
            blocks = env$County, nperm = 999)

plant.div <-with(env, adonis2(plant.dist ~ env$Lupin, data= plants.transformed, method="bray", permutations = perm))
plant.div

plant.div <-with(env, adonis2(plant.dist ~ env$Lupin * env$County, data= plants.transformed, method="bray"))
plant.div

# NMDS with all data----
# Pool abundances per site 
raw.data <- read.csv2(here("data", "Raw data.csv" ))

data.1 <- raw.data %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(Parnr)) %>% 
  mutate(Lupin.NoLupin = as.factor(Lupin.NoLupin)) %>% 
  group_by(County, Site, parnr, Lupin.NoLupin, Species) %>%
  summarise(Occ= sum(Occurrence))
data.1

data.pooled <- data.1 %>% 
  group_by(County, Site, Lupin.NoLupin, Species) %>% 
  summarise(tot.abund = sum(Occ))
data.pooled

# Turn into wide format
df <- data.pooled %>% 
  group_by(County, Site,Lupin.NoLupin) %>% 
  spread(., key="Species", value = "tot.abund") %>% 
  replace(is.na(.), 0)
df # saved as dataframe.pooled.csv

# Load dataframe
data <- read.csv2(here("data","dataframe.pooled.csv"))
Lupin <- data[,3]
# transform the data to remove the influence of the most abundant species
plants <- sqrt(data[,4:222])#square root transform

NMDS <- metaMDS(plants, distance="bray", k=2, try = 100, trymax = 1000, autotransform=TRUE)
NMDS
stressplot(NMDS) #Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions. 

##pull points from MDS
NMDS1 <- NMDS$points[,1] ##also found using scores(NMDS)
NMDS2 <- NMDS$points[,2]
plant.plot<-cbind(plants, NMDS1, NMDS2, Lupin)

# quick plot ordination
p<-ggplot(plant.plot, aes(NMDS1, NMDS2, color=Lupin))+
  geom_point(position=position_jitter(.1), shape=3)+##separates overlapping points
  stat_ellipse(type='t',size =1)+ ##draws 95% confidence interval ellipses
  theme_minimal()
p

# complicated plot with labels
# Create figure with ggplot
plant.scores <- as.data.frame(scores(NMDS))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
plant.scores$Site <- rownames(plant.scores)   # create a column of site names, from the rownames of data.scores
plant.scores$Lupin <- Lupin  #  add the grp variable created earlier
head(plant.scores)  #look at the data

sp.scores <- as.data.frame(scores(NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
sp.scores$species <- rownames(sp.scores)  # create a column of species, from the rownames of species.scores
head(sp.scores)  #look at the data

Lupin <- plant.scores[plant.scores$Lupin == "Lupin", ][chull(plant.scores[plant.scores$Lupin == 
                                                                               "Lupin", c("NMDS1", "NMDS2")]), ]  # hull values 
NoLupin <- plant.scores[plant.scores$Lupin == "NoLupin", ][chull(plant.scores[plant.scores$Lupin == 
                                                                         "NoLupin", c("NMDS1", "NMDS2")]), ]  # hull 
 

hull.plants <- rbind(Lupin, NoLupin)  #combine plots
hull.plants

gg.lupin<- ggplot() + 
  geom_polygon(data=hull.plants,aes(x=NMDS1,y=NMDS2,fill=Lupin , group=Lupin),alpha=0.25) + # add the convex hulls
  geom_text(data=sp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=4) +  # add the species labels
  geom_point(data=plant.scores,aes(x=NMDS1,y=NMDS2,shape=Lupin,colour=Lupin),size=2) + # add the point markers
  #geom_text(data=plant.scores,aes(x=NMDS1,y=NMDS2,label=Lupin),size=1,vjust=0) +  # add the site labels
  scale_colour_manual(values=c("Lupin" = "#fdae61", "NoLupin" = "#2c7bb6")) +
  scale_fill_manual(values=c("Lupin" = "#fdae61", "NoLupin" = "#2c7bb6")) +
  scale_shape_discrete("Lupin") +
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
gg.lupin

# Mantel test for spatial autocorrelation ----
# Load coordinate data
spatial.data <- read.csv2(here("data", "sites.spatial.csv"))

#remove duplicates
spatial <- spatial.data %>% 
  group_by(Site, lat, lng) %>% 
  summarise()
spatial # I have to manually remove duplicates for some sites

spatial.data <- read.csv2(here("data", "1coordset.csv"))

geo <- data.frame(as.numeric(spatial.data$lng), as.numeric(spatial.data$lat))

# calcluate distance matrix
d.geo = distm(geo, fun = distHaversine)
dist.geo = as.dist(d.geo)

#abundance data frame - bray curtis dissimilarity
plant.dist <- vegdist(plants, method='bray') #abundance based

#abundance vs geographic 
abund_geo  = mantel(plant.dist, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo

# PERMANOVA ----
# NUll hypothesis: Groups do not differ in spread or position in multivariate space.
plants.transformed <- sqrt(data[,4:222])#square root transform# square root transform
plants <- cbind(data[,1:3], plants.transformed)
plant.dist <- vegdist(plants[,4:222], method='bray') #abundance based
env <- plants[,1:3]

# Does the community composition differ in plots with/without lupin?
# set the blocks (sites)
perm <- how(nperm = 999)
setStrata(perm) <- with(plants, Site)

plant.div <-with(env, adonis2(plant.dist ~ Lupin, data= plants.transformed, permutations = perm, method="bray"))
plant.div

plant.div2 <- adonis2(plant.dist ~ Lupin, data= plants.transformed, method="bray", strata = Site)
plant.div2

