rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(permute)
library(geosphere)
library(ncf)
library(glmmTMB)
library(effects)
library(betapart)

# Species richness ----
# Load the data
tidy.data <- read.csv2(here("data/Final Datasets", "final.tidy.csv"))

# Restructure
db <- tidy.data %>% 
  filter(Species !=  "Lupinus polyphyllus") %>% 
  mutate(Site = as.factor(Site)) %>%
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>%
  rename(Region = County) %>% 
  mutate(Region=as.factor(Region))  
str(db)

# Species richness
data.rich <- db %>% 
  group_by(Region, Site, Lupin, N, E, lat, lng) %>% 
  summarise(rich = n())
str(data.rich)

# Are sites that are near each other more similar in terms of species richness?
# read data set wit only one set of coords
data.rich <- read.csv2(here("data/One set of coords for spatial autocorrelation", "rich.coords.csv"))
# Run model
model.rich <- glmmTMB(rich ~ Lupin * Region + (1|Site), family= "poisson", data=data.rich)
summary(model.rich)
plot(allEffects(model.rich))
ncf.sprich <- spline.correlog(data.rich$N, data.rich$E, latlon = F, 
                              resid(model.rich, type  = "pearson"), resamp=1000, na.rm = TRUE)
plot1 <- plot(ncf.sprich, cex.lab=1.5, xlab="Distance (m)") #, main ="New data")


# Effective number of species ----
db.shannon <- read.csv2(here("data/Final Datasets", "db.shannon.csv"))

# Restructure
shannon <- db.shannon %>% 
  mutate(Region =as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>%
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) 
str(shannon)

# Are sites that are near each other more similar in terms of the effective number of species?
# Run model
model.shannon <- glmmTMB(shannon ~ Lupin * Region + (1|Site), family= "gaussian" (link="log"), data=db.shannon)
summary(model.shannon)
plot(allEffects(model.shannon))
ncf.effective <- spline.correlog(shannon$N, shannon$E, latlon = F, 
                                 resid(model.shannon , type  = "pearson"), resamp=1000, na.rm = TRUE)
plot(ncf.effective, cex.lab=1.5, xlab="Distance (m)")#, main ="New data")
 # Plot together ----
par(mfcol=c(1,2))
plot(ncf.sprich, cex.lab=1.5, xlab="Distance (m)", main = "a", adj=0)
plot(ncf.effective, cex.lab=1.5, xlab="Distance (m)", main = "b", adj=0)

# -----

# Species composition ----
rm(list=ls())
# Are sites that are near each other more similar in terms of species composition?
## load data
df <- read.csv2(here("data/Final Datasets", "final.dataframe.coords.2.csv"))

# compute Bray Curtis dissimilarity between sites
bray.distance <- bray.part(df[,-c(1:8)])$bray # selects dissimilarity matrix accounting for total abundance-based dissimilarity between sites
str(bray.distance)

# calculate distance matrix
geo <- data.frame(as.numeric(df$N), as.numeric(df$E))

colnames(geo)[1] = "N"
colnames(geo)[2] = "E"

# compute spline correlogram
ncf.bray <- spline.correlog(x=geo$N, y=geo$E, z= as.matrix(bray.distance), latlon = F, resamp=1000, na.rm = TRUE)
plot(ncf.bray, cex.lab=1.5, xlab="Distance (m)")#, main ="New data")


# Mantel test for spatial autocorrelation
df <- read.csv2(here("data/Final Datasets", "final.dataframe.coords.csv"))
comm <- df[,-c(1:8)]

geo <- data.frame(as.numeric(df$lng), as.numeric(df$lat))
d.geo <- distm(geo, fun = distHaversine)
dist.geo <- as.dist(d.geo)
#check that the distance is correct
max(dist.geo)

#abundance data frame - bray curtis dissimilarity
plant.dist <- vegdist(comm, method='bray') #abundance based

#abundance vs geographic 
abund_geo  = mantel(plant.dist, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo

abund_geo  = mantel(plant.dist, dist.geo, method = "pearson", permutations = 9999, na.rm = TRUE)
abund_geo

