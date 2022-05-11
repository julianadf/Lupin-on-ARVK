rm(list=ls())
library(here)
library(tidyverse)
library(sf)
library(spdep)
library(tmap)

#Load the shapefile
site.richness <- read_sf(here("data/Shapefiles", "Site.shp" ))
names(site.richness)
hist(site.richness$rich, main=NULL)#check for outliers or skewedness
boxplot(site.richness$rich, horizontal=TRUE) # no outliers

tm_shape(site.richness) + tm_bubbles(col="rich", style="quantile", n=8, palette="Greens") +
  tm_legend(outside=TRUE)
coords <- st_centroid(st_geometry(site.richness), of_largest_polygon=TRUE)
test <- knearneigh(site.richness, k=4)
nb <- knn2nb(knearneigh(site.richness))

plot(st_geometry(site.richness), border="grey")
plot(knn2nb(test), coords, add=TRUE)
