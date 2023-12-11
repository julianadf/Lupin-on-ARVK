rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(permute)
library(geosphere)
library(funrar)
library(ggplot2)
library(labdsv)

# Load the  dataframe ----
df <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))

# Restructure
df.res <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region =as.factor(Region))  
str(df.res)

# Square-root transformed abundances ----
tot.abundance <- df.res[,5:ncol(df.res)]
sqrt.abund <- sqrt(tot.abundance)
Lupin <- df.res$Lupin
Region <- df.res$Region

# With Bray-Curtis distances
NMDS <- metaMDS(tot.abundance, distance = "bray",  try = 100, trymax = 1000, maxit = 20000, k=3)
stressplot(NMDS, main="Bray-Curtis")
gof <- goodness(NMDS)
plot(NMDS, type="t", main="Goodness of fit")
points(NMDS, display = "sites", cex=gof*300)

# Create figure with ggplot
site.scores <- as.data.frame(scores(NMDS)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
site.scores$Lupin <- rownames(scores)  # create a column of site names, from the rownames of data.scores
site.scores$Lupin <- Lupin
head(site.scores)  #look at the data

sp.scores <- as.data.frame(scores(NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
sp.scores$species <- rownames(sp.scores)  # create a column of species, from the rownames of species.scores
head(sp.scores)  #look at the data

lupin <- site.scores[site.scores$Lupin == "Yes", ][chull(site.scores[site.scores$Lupin =="Yes", c("NMDS1", "NMDS2")]), ]  # hull values 
no.lupin <- site.scores[site.scores$Lupin == "No", ][chull(site.scores[site.scores$Lupin =="No", c("NMDS1", "NMDS2")]), ]  # hull 


hull <- rbind(lupin, no.lupin)  #combine treatments
hull

# main text ----
gg.1<- ggplot() + 
  geom_polygon(data=hull,aes(x=NMDS1,y=NMDS2,fill=Lupin),alpha=0.25) + # add the convex hulls
  #geom_text(data=sp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=3) + # add the species labels
  geom_point(data= site.scores,aes(x=NMDS1,y=NMDS2,shape=Lupin,colour=Lupin),size=2) +  # add the point markers
  #geom_text(data=site.scores,aes(x=NMDS1,y=NMDS2,label=site),size=2,vjust=0) +  # add the site labels
  scale_color_manual("Lupin", values = c("#008837", "#7b3294")) +
  scale_fill_manual("Lupin", values = c("#008837", "#7b3294")) +
  scale_shape_discrete("Lupin") +
  coord_equal() +
  ggtitle("a") +
  theme(legend.position = "bottom", legend.title = element_text(size=24), legend.text = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=26), axis.text.y.left = element_text(size=26),
        axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28), 
        plot.title = element_text(hjust = 0, size= 26))
gg.1 

# NMDS per region ---- 
Region <- df.res$Region

site.scores <- as.data.frame(scores(NMDS)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
site.scores$Region <- rownames(scores)  # create a column of site names, from the rownames of data.scores
site.scores$Region <- Region
head(site.scores)  #look at the data

sp.scores <- as.data.frame(scores(NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
sp.scores$species <- rownames(sp.scores)  # create a column of species, from the rownames of species.scores
head(sp.scores)  #look at the data

east <- site.scores[site.scores$Region == "East region", ][chull(site.scores[site.scores$Region =="East region", c("NMDS1", "NMDS2")]), ]  # hull values 
west <- site.scores[site.scores$Region == "West region", ][chull(site.scores[site.scores$Region =="West region", c("NMDS1", "NMDS2")]), ]  # hull 


hull <- rbind(east, west)  #combine treatments
hull

# main text
gg.2 <- ggplot() + 
  geom_polygon(data=hull,aes(x=NMDS1,y=NMDS2,fill=Region),alpha=0.25) + # add the convex hulls
  #geom_text(data=sp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=3) + # add the species labels
  geom_point(data= site.scores,aes(x=NMDS1,y=NMDS2,shape=Region,colour=Region),size=2) +  # add the point markers
  #geom_text(data=site.scores,aes(x=NMDS1,y=NMDS2,label=site),size=2,vjust=0) +  # add the site labels
  scale_color_manual("Region", values = c("#d7191c", "#2c7bb6")) +
  scale_fill_manual("Region", values = c("#d7191c", "#2c7bb6")) +
  scale_shape_discrete("Region") +
  coord_equal() +
  ggtitle("b")+
  theme(legend.position = "bottom", legend.title = element_text(size=24), legend.text = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=26), axis.text.y.left = element_text(size=26),
        axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28), 
        plot.title = element_text(hjust = 0, size= 26))
gg.2


# Load the unbalanced dataframe ----
df <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))

# Restructure
df.res <- df %>% 
  #dplyr::select(-Lupinus.polyphyllus) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region =as.factor(Region))  
str(df.res)

# Square-root transformed abundances ----
tot.abundance <- df.res[,5:ncol(df.res)]
sqrt.abund <- sqrt(tot.abundance)
Lupin <- df.res$Lupin
Region <- df.res$Region

# With square root abundances
NMDS <- metaMDS(sqrt.abund, try = 100, trymax = 1000, maxit = 20000, k=3)
stressplot(NMDS, main="sqrt(abundance)")
gof <- goodness(NMDS)
plot(NMDS, type="t", main="Goodness of fit")
points(NMDS, display = "sites", cex=gof*300)
# Stress = 0.2553801 

# With Bray-Curtis distances
NMDS <- metaMDS(tot.abundance, distance = "bray",  try = 100, trymax = 1000, maxit = 20000, k=4)
stressplot(NMDS, main="Bray-Curtis")
gof <- goodness(NMDS)
plot(NMDS, type="t", main="Goodness of fit")
points(NMDS, display = "sites", cex=gof*300)


# Check stress at different dimensions
# Code Lutz sent
n=10
stress <- vector(length = n)
for (i in 1:n){
  stress[i] <- metaMDS(tot.abundance, distance="bray", k=i)$stress
}
names(stress) <- paste0(1:n, "Dim")

par(mar = c(3.5,3.5,1,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 2)
barplot(stress, ylab = "stress")

# Spiderplot
plot(NMDS)
spider <- ordispider(NMDS, groups = Lupin)

# Create figure with ggplot
site.scores <- as.data.frame(scores(NMDS)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
site.scores$Lupin <- rownames(scores)  # create a column of site names, from the rownames of data.scores
site.scores$Lupin <- Lupin
head(site.scores)  #look at the data

sp.scores <- as.data.frame(scores(NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
sp.scores$species <- rownames(sp.scores)  # create a column of species, from the rownames of species.scores
head(sp.scores)  #look at the data

lupin <- site.scores[site.scores$Lupin == "Yes", ][chull(site.scores[site.scores$Lupin =="Yes", c("NMDS1", "NMDS2")]), ]  # hull values 
no.lupin <- site.scores[site.scores$Lupin == "No", ][chull(site.scores[site.scores$Lupin =="No", c("NMDS1", "NMDS2")]), ]  # hull 


hull <- rbind(lupin, no.lupin)  #combine treatments
hull

# main text ----
gg.1<- ggplot() + 
  geom_polygon(data=hull,aes(x=NMDS1,y=NMDS2,fill=Lupin),alpha=0.25) + # add the convex hulls
  #geom_text(data=sp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=3) + # add the species labels
  geom_point(data= site.scores,aes(x=NMDS1,y=NMDS2,shape=Lupin,colour=Lupin),size=2) +  # add the point markers
  #geom_text(data=site.scores,aes(x=NMDS1,y=NMDS2,label=site),size=2,vjust=0) +  # add the site labels
  scale_color_manual("Lupin", values = c("#008837", "#7b3294")) +
  scale_fill_manual("Lupin", values = c("#008837", "#7b3294")) +
  scale_shape_discrete("Lupin") +
  coord_equal() +
  ggtitle("a") +
  theme(legend.position = "bottom", legend.title = element_text(size=24), legend.text = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=26), axis.text.y.left = element_text(size=26),
        axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28), 
        plot.title = element_text(hjust = 0, size= 26))
gg.1 

# NMDS per region ---- 
Region <- df.res$Region

site.scores <- as.data.frame(scores(NMDS)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
site.scores$Region <- rownames(scores)  # create a column of site names, from the rownames of data.scores
site.scores$Region <- Region
head(site.scores)  #look at the data

sp.scores <- as.data.frame(scores(NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
sp.scores$species <- rownames(sp.scores)  # create a column of species, from the rownames of species.scores
head(sp.scores)  #look at the data

east <- site.scores[site.scores$Region == "East region", ][chull(site.scores[site.scores$Region =="East region", c("NMDS1", "NMDS2")]), ]  # hull values 
west <- site.scores[site.scores$Region == "West region", ][chull(site.scores[site.scores$Region =="West region", c("NMDS1", "NMDS2")]), ]  # hull 


hull <- rbind(east, west)  #combine treatments
hull

# main text
gg.2 <- ggplot() + 
  geom_polygon(data=hull,aes(x=NMDS1,y=NMDS2,fill=Region),alpha=0.25) + # add the convex hulls
  #geom_text(data=sp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=3) + # add the species labels
  geom_point(data= site.scores,aes(x=NMDS1,y=NMDS2,shape=Region,colour=Region),size=2) +  # add the point markers
  #geom_text(data=site.scores,aes(x=NMDS1,y=NMDS2,label=site),size=2,vjust=0) +  # add the site labels
  scale_color_manual("Region", values = c("#d7191c", "#2c7bb6")) +
  scale_fill_manual("Region", values = c("#d7191c", "#2c7bb6")) +
  scale_shape_discrete("Region") +
  coord_equal() +
  ggtitle("b")+
  theme(legend.position = "bottom", legend.title = element_text(size=24), legend.text = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=26), axis.text.y.left = element_text(size=26),
        axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28), 
        plot.title = element_text(hjust = 0, size= 26))
gg.2



# Load the dataframe and remove all species that have 3 occcurrences or less----
df <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))

# Restructure
df.res <- df %>% 
  #dplyr::select(-Lupinus.polyphyllus) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region =as.factor(Region))  
str(df.res)

# Square-root transformed abundances ----
tot.abundance <- df.res[,5:ncol(df.res)] # ncol=129

tot.abund.morethan3 <- 
  tot.abundance %>%
  select_if(~max(., na.rm = TRUE) >= 3)
tot.abund.morethan3 #ncol= 98

sqrt.abund <- sqrt(tot.abund.morethan3)
Lupin <- df.res$Lupin
Region <- df.res$Region

# With square root abundances
NMDS <- metaMDS(sqrt.abund, try = 100, trymax = 1000, maxit = 20000, k=4)
stressplot(NMDS, main="sqrt(abundance)")
gof <- goodness(NMDS)
plot(NMDS, type="t", main="Goodness of fit")
points(NMDS, display = "sites", cex=gof*300)
# Stress = 0.2553801 

# With Bray-Curtis distances
NMDS <- metaMDS(tot.abund.morethan3, distance = "bray",  try = 100, trymax = 1000, maxit = 20000, k=3)
stressplot(NMDS, main="Bray-Curtis")
gof <- goodness(NMDS)
plot(NMDS, type="t", main="Goodness of fit")
points(NMDS, display = "sites", cex=gof*300)


# Check stress at different dimensions
# Code Lutz sent
n=10
stress <- vector(length = n)
for (i in 1:n){
  stress[i] <- metaMDS(tot.abundance, distance="bray", k=i)$stress
}
names(stress) <- paste0(1:n, "Dim")

par(mar = c(3.5,3.5,1,1), mgp = c(2, 0.6, 0), cex = 0.8, las = 2)
barplot(stress, ylab = "stress")

# Spiderplot
plot(NMDS)
spider <- ordispider(NMDS, groups = Lupin)

# Create figure with ggplot
site.scores <- as.data.frame(scores(NMDS)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
site.scores$Lupin <- rownames(scores)  # create a column of site names, from the rownames of data.scores
site.scores$Lupin <- Lupin
head(site.scores)  #look at the data

sp.scores <- as.data.frame(scores(NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
sp.scores$species <- rownames(sp.scores)  # create a column of species, from the rownames of species.scores
head(sp.scores)  #look at the data

lupin <- site.scores[site.scores$Lupin == "Yes", ][chull(site.scores[site.scores$Lupin =="Yes", c("NMDS1", "NMDS2")]), ]  # hull values 
no.lupin <- site.scores[site.scores$Lupin == "No", ][chull(site.scores[site.scores$Lupin =="No", c("NMDS1", "NMDS2")]), ]  # hull 


hull <- rbind(lupin, no.lupin)  #combine treatments
hull

# main text ----
gg.1<- ggplot() + 
  geom_polygon(data=hull,aes(x=NMDS1,y=NMDS2,fill=Lupin),alpha=0.25) + # add the convex hulls
  #geom_text(data=sp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=3) + # add the species labels
  geom_point(data= site.scores,aes(x=NMDS1,y=NMDS2,shape=Lupin,colour=Lupin),size=2) +  # add the point markers
  #geom_text(data=site.scores,aes(x=NMDS1,y=NMDS2,label=site),size=2,vjust=0) +  # add the site labels
  scale_color_manual("Lupin", values = c("#008837", "#7b3294")) +
  scale_fill_manual("Lupin", values = c("#008837", "#7b3294")) +
  scale_shape_discrete("Lupin") +
  coord_equal() +
  ggtitle("a") +
  theme(legend.position = "bottom", legend.title = element_text(size=24), legend.text = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=26), axis.text.y.left = element_text(size=26),
        axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28), 
        plot.title = element_text(hjust = 0, size= 26))
gg.1 

# NMDS per region ---- 
Region <- df.res$Region

site.scores <- as.data.frame(scores(NMDS)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
site.scores$Region <- rownames(scores)  # create a column of site names, from the rownames of data.scores
site.scores$Region <- Region
head(site.scores)  #look at the data

sp.scores <- as.data.frame(scores(NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
sp.scores$species <- rownames(sp.scores)  # create a column of species, from the rownames of species.scores
head(sp.scores)  #look at the data

east <- site.scores[site.scores$Region == "East region", ][chull(site.scores[site.scores$Region =="East region", c("NMDS1", "NMDS2")]), ]  # hull values 
west <- site.scores[site.scores$Region == "West region", ][chull(site.scores[site.scores$Region =="West region", c("NMDS1", "NMDS2")]), ]  # hull 


hull <- rbind(east, west)  #combine treatments
hull

# main text
gg.2 <- ggplot() + 
  geom_polygon(data=hull,aes(x=NMDS1,y=NMDS2,fill=Region),alpha=0.25) + # add the convex hulls
  #geom_text(data=sp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=3) + # add the species labels
  geom_point(data= site.scores,aes(x=NMDS1,y=NMDS2,shape=Region,colour=Region),size=2) +  # add the point markers
  #geom_text(data=site.scores,aes(x=NMDS1,y=NMDS2,label=site),size=2,vjust=0) +  # add the site labels
  scale_color_manual("Region", values = c("#d7191c", "#2c7bb6")) +
  scale_fill_manual("Region", values = c("#d7191c", "#2c7bb6")) +
  scale_shape_discrete("Region") +
  coord_equal() +
  ggtitle("b")+
  theme(legend.position = "bottom", legend.title = element_text(size=24), legend.text = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=26), axis.text.y.left = element_text(size=26),
        axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28), 
        plot.title = element_text(hjust = 0, size= 26))
gg.2
