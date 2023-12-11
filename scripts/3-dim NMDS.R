rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(permute)
library(geosphere)
library(funrar)
library(ggplot2)
library(labdsv)
library(ggalt)
library(ggh4x)
library(patchwork)
library(permute)
library(vegan3d)
library(rgl)

# Load the  dataframe ----
df <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))

# Restructure
df.res <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  dplyr::select(-Rubus.fruticosus) %>% 
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
# remove all species with 3 occurrences or less
#tot.abund2 <- tot.abundance %>%  select(all_of(names(.)[1:3]), where(~ is.numeric(.) && 
                                                                       #sum(., na.rm = TRUE) > 3))
#sqrt.abund <- sqrt(tot.abundance)
Lupin <- df.res$Lupin
Region <- df.res$Region

# test
# data(tot.abundance, Lupin)
# d <- vegdist(tot.abundance)
# m <- metaMDS(d)
# cl <- hclust(d, "aver")
# orditree3d(m, cl, pch=16, col=cutree(cl, 3))
# 
# ordiplot3d(NMDS, display="sites", col=c("blue", "red"))
# ordirgl(NMDS, display = "sites")
# ordispider(NMDS, groups = "Lupin")

### NMDS ###
# ndms.plots <- c()
# nmds.stress <- c()

NMDS <- metaMDS(dist.bc, maxit = 20000, maxtry = 100, k= 3, trace = 0, distance = "bray")
# ordiplot(NMDS)
ordiplot (NMDS, type = "points", display = "sites")
ordispider(NMDS, groups = df.res$Lupin, col= c("#008837", "#7b3294"), label =T)
ordihull(NMDS, groups = df.res$Lupin, col= c("#008837", "#7b3294"), label =F)
specfit <- envfit(NMDS ~  Region, df.res)
plot(specfit)

dist <- beta.pair(dist.bc, index.family="sorensen")
bd <- betadisper(dist[[3]],lupin)
plot(bd, main="sorensen") #The betadisper plot indicates that there is NO difference in species compositions (as the NMDS)

# boxplot to observe the distance of values of beta diversity of each treatment in relation to their
# centroids (basically, this indicates homogeneity in how communities of a given treatment differ 
# from each other). 
boxplot(bd)
# Perform an ANOVA to test if treatments are significantly different.
anova(bd) #almost, but treatments do not differ significantly in relation to how communities vary from each other.

beta.tot <- beta.multi(presabs.lupin)
gg.beta <- as.data.frame(beta.tot)


# Quick and dirty plot -----
# NMDS1 <- NMDS$points[,1] ##also found using scores(NMDS)
# NMDS2 <- NMDS$points[,2]
# envfit(NMDS ~ Lupin, data=df.res)
# envfit(NMDS ~ Region, data=df.res)
# 
# plant.plot <-cbind(tot.abundance, NMDS1, NMDS2, Lupin)
# 
# quick.dirty <- ggplot(plant.plot, aes(NMDS1, NMDS2, color=Lupin)) +
#   geom_point(position=position_jitter(.1), shape=3) + ##separates overlapping points
#   stat_ellipse(type='t',size =1)+ ##draws 95% confidence interval ellipses
#   theme_minimal()
# quick.dirty

# ------


# prepare scores
site.scores <- as.data.frame(scores(NMDS)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
site.scores$Lupin <- rownames(scores)  # create a column of site names, from the rownames of data.scores
site.scores$Lupin <- Lupin
#site.scores$Region <- Region
head(site.scores)  #look at the data

sp.scores <- as.data.frame(scores(NMDS, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
sp.scores$species <- rownames(sp.scores)  # create a column of species, from the rownames of species.scores
head(sp.scores)  #look at the data

lupin <- site.scores[site.scores$Lupin == "Yes", ][chull(site.scores[site.scores$Lupin =="Yes", c("NMDS1", "NMDS2")]), ]  # hull values 
no.lupin <- site.scores[site.scores$Lupin == "No", ][chull(site.scores[site.scores$Lupin =="No", c("NMDS1", "NMDS2")]), ]  # hull 


hull <- rbind(lupin, no.lupin)  #combine treatments
hull

p1 <- ggplot(
  data = 
    bind_rows(
      site.scores %>% 
        select(NMDS1,NMDS2, Lupin) %>% 
        rename(x = NMDS2, y = NMDS1) %>% 
        mutate(Axis_x = "NMDS2", Axis_y = "NMDS1"),
      site.scores %>% 
        select(NMDS1,NMDS3, Lupin) %>% 
        rename(x = NMDS3, y = NMDS1) %>% 
        mutate(Axis_x = "NMDS3", Axis_y = "NMDS1"),
      site.scores %>% 
        select(NMDS3,NMDS2, Lupin) %>% 
        rename(x = NMDS3, y = NMDS2) %>% 
        mutate(Axis_x = "NMDS3", Axis_y = "NMDS2")
    ), # prepare data for showing all three dimensions
  aes(x=x,y=y,shape=Lupin,colour=Lupin)
) +
  geom_point(size=3) + # points
  geom_encircle(
    aes(fill = Lupin),
    alpha=0.25,
    s_shape = 1,
    expand = 0,
    spread = 0
  ) + # polygons showing habitat types
  facet_manual(
    Axis_x~Axis_y,  
    design = matrix(c(1,NA,2,3), ncol = 2), 
    respect = T, 
    axes = "all",
    strip = strip_split(
      position = c("bottom", "left"),
      text_y = list(element_text(),element_text(color = "white"))
    )
   ) +
   theme_classic()  +
    scale_y_continuous("") + 
    scale_x_continuous("") +
    scale_color_manual("Lupin", values = c("#008837", "#7b3294"), labels = c("No", "Yes")) +
    scale_fill_manual("Lupin", values = c("#008837", "#7b3294"), labels = c("No", "Yes")) +
    #scale_shape_discrete("Lupin", labels = c("Yes", "No")) +
   theme(legend.position = c(.2, .2),
         legend.title= element_blank(),
         panel.background = element_rect(colour = 'black'),
         strip.background = element_blank(),
         strip.placement = "outside") 
   #ggtitle("Stress = 0.1957066 , dimensions =3, Bray-Curtis") #, subtitle = paste("Stress =", round(nmds.stress[1,2],2)))
p1

# Region
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

p2 <- ggplot(
 data =
   bind_rows(
     site.scores %>%
       select(NMDS1,NMDS2, Region) %>%
       rename(x = NMDS2, y = NMDS1) %>%
       mutate(Axis_x = "NMDS2", Axis_y = "NMDS1"),
     site.scores %>%
       select(NMDS1,NMDS3, Region) %>%
       rename(x = NMDS3, y = NMDS1) %>%
       mutate(Axis_x = "NMDS3", Axis_y = "NMDS1"),
     site.scores %>%
       select(NMDS3,NMDS2, Region) %>%
       rename(x = NMDS3, y = NMDS2) %>%
       mutate(Axis_x = "NMDS3", Axis_y = "NMDS2")
   ), # prepare data for showing all three dimensions
 aes(x=x,y=y,shape=Region,colour=Region)
) +
 geom_point(size=3) + # points
 geom_encircle(
   aes(fill = Region),
   alpha=0.25,
   s_shape = 1,
   expand = 0,
   spread = 0
 ) + # polygons showing habitat types
 facet_manual(
   Axis_x~Axis_y,
   design = matrix(c(1,NA,2,3), ncol = 2),
   respect = T,
   axes = "all",
   strip = strip_split(
     position = c("bottom", "left"),
     text_y = list(element_text(),element_text(color = "white"))
   )
 ) +
 theme_classic()  +
 scale_y_continuous("") +
 scale_x_continuous("") +
 scale_color_manual("Region", values = c("#d7191c", "#2c7bb6")) +
                    #labels = c("West region", "East region")) +
 scale_fill_manual("Region", values = c("#d7191c", "#2c7bb6")) +
                   #labels = c("West region", "East region")) +
 #scale_shape_discrete("Region", labels = c("West region", "East region")) +
 theme(legend.position = c(.2, .2),
       legend.title= element_blank(),
       panel.background = element_rect(colour = 'black'),
       strip.background = element_blank(),
       strip.placement = "outside")
#ggtitle("Stress = 0.1957066 , dimensions =3, Bray-Curtis") #, subtitle = paste("Stress =", round(nmds.stress[1,2],2)))
p2


# PerMANOVA ----
tot.abundance <- df.res[,5:ncol(df.res)]
bray.curtis <- vegdist(tot.abundance, method="bray")
Site <- df.res$Site
# NUll hypothesis: Groups do not differ in spread or position in multivariate space.
perm <- adonis2(tot.abundance ~ Lupin + Region, data = df.res)
perm


perm <- how(within = Within(type = "free"),
            blocks = env$County, nperm = 999, strata= Site)

h <- how(plots = Plots(strata = env$Site, type = "free"), blocks = env$County,
         nperm = 999)

plant.div <-with(env, adonis2(plant.dist ~ env$Lupin,  data= plants.transformed, method="bray", permutations = h))
plant.div

plant.div2 <-with(env, adonis2(plant.dist ~ env$Lupin * env$County, data= plants.transformed, method="bray", permutations = h))
plant.div2


test <- adonis2(plant.dist ~ env$Lupin, data= plants.transformed, method="bray")
test      

test2 <- nested.npmanova(plant.dist ~ Lupin + Site, data = env, permutations = 999)
test2


# Each region separately -----
df <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))

# Uppland -----
df.upp <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  dplyr::select(-Rubus.fruticosus) %>% 
  filter(County == "Uppsala") %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  #mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region =as.factor(Region))  
str(df.upp)

# Square-root transformed abundances ----
tot.abundance <- df.upp[,5:ncol(df.upp)]
Lupin <- df.upp$Lupin
# remove all species with 3 occurrences or less
tot.abund2 <- tot.abundance %>%  select(all_of(names(.)[1:3]), where(~ is.numeric(.) && 
  sum(., na.rm = TRUE) > 3))


### NMDS ###
NMDS.upp <- metaMDS(tot.abundance, maxit = 20000, maxtry = 100, k= 3, trace = 0, distance = "bray")

# plot ----
site.scores <- as.data.frame(scores(NMDS.upp)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
site.scores$Lupin <- rownames(scores)  # create a column of site names, from the rownames of data.scores
site.scores$Lupin <- Lupin
head(site.scores)  #look at the data

sp.scores <- as.data.frame(scores(NMDS.upp, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
sp.scores$species <- rownames(sp.scores)  # create a column of species, from the rownames of species.scores
head(sp.scores)  #look at the data

lupin <- site.scores[site.scores$Lupin == "Yes", ][chull(site.scores[site.scores$Lupin =="Yes", c("NMDS1", "NMDS2")]), ]  # hull values 
no.lupin <- site.scores[site.scores$Lupin == "No", ][chull(site.scores[site.scores$Lupin =="No", c("NMDS1", "NMDS2")]), ]  # hull 


hull <- rbind(lupin, no.lupin)  #combine treatments
hull

# main text
gg<- ggplot() + 
  geom_polygon(data=hull,aes(x=NMDS1,y=NMDS2,fill=Lupin),alpha=0.25) + # add the convex hulls
  #geom_text(data=sp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=3) + # add the species labels
  geom_point(data= site.scores,aes(x=NMDS1,y=NMDS2,shape=Lupin,colour=Lupin),size=2) +  # add the point markers
  #geom_text(data=site.scores,aes(x=NMDS1,y=NMDS2,label=site),size=2,vjust=0) +  # add the site labels
  scale_color_manual("Lupin", values = c("#008837", "#7b3294")) +
  scale_fill_manual("Lupin", values = c("#008837", "#7b3294")) +
  scale_shape_discrete("Lupin") +
  coord_equal() +
  theme(legend.position = "right", legend.title = element_text(size=24), legend.text = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=26), axis.text.y.left = element_text(size=26),
        axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28))
gg + ggtitle("East region - Uppland. stress = 0.13, 3 dimensions, Bray-Curtis\npermanova not significant")

#permanova
tot.abund2 <- tot.abundance %>%  select(all_of(names(.)[1:3]), where(~ is.numeric(.) && 
                                                                       sum(., na.rm = TRUE) > 3))
Site <- df.upp$Site
# NUll hypothesis: Groups do not differ in spread or position in multivariate space.
perm <- adonis2(tot.abund2 ~ Lupin, data = df.upp)
perm

# -----

df <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))

# Restructure
df.värm <- df %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  filter(County == "Värmland") %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  #mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region =as.factor(Region))  
str(df.värm)

# Square-root transformed abundances ----
tot.abundance <- df.värm[,5:ncol(df.värm)]
#sqrt.abund <- sqrt(tot.abundance)
Lupin <- df.värm$Lupin

### NMDS ###
NMDS.värm <- metaMDS(tot.abundance, maxit = 20000, maxtry = 100, k= 3, trace = 0, distance = "bray")

# plot ----
site.scores <- as.data.frame(scores(NMDS.värm)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
site.scores$Lupin <- rownames(scores)  # create a column of site names, from the rownames of data.scores
site.scores$Lupin <- Lupin
head(site.scores)  #look at the data

sp.scores <- as.data.frame(scores(NMDS.värm, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
sp.scores$species <- rownames(sp.scores)  # create a column of species, from the rownames of species.scores
head(sp.scores)  #look at the data

lupin <- site.scores[site.scores$Lupin == "Yes", ][chull(site.scores[site.scores$Lupin =="Yes", c("NMDS1", "NMDS2")]), ]  # hull values 
no.lupin <- site.scores[site.scores$Lupin == "No", ][chull(site.scores[site.scores$Lupin =="No", c("NMDS1", "NMDS2")]), ]  # hull 

hull <- rbind(lupin, no.lupin)  #combine treatments
hull

# main text
gg<- ggplot() + 
  geom_polygon(data=hull,aes(x=NMDS1,y=NMDS2,fill=Lupin),alpha=0.25) + # add the convex hulls
  #geom_text(data=sp.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5, size=3) + # add the species labels
  geom_point(data= site.scores,aes(x=NMDS1,y=NMDS2,shape=Lupin,colour=Lupin),size=2) +  # add the point markers
  #geom_text(data=site.scores,aes(x=NMDS1,y=NMDS2,label=site),size=2,vjust=0) +  # add the site labels
  scale_color_manual("Lupin", values = c("#008837", "#7b3294")) +
  scale_fill_manual("Lupin", values = c("#008837", "#7b3294")) +
  scale_shape_discrete("Lupin") +
  coord_equal() +
  theme(legend.position = "right", legend.title = element_text(size=24), legend.text = element_text(size=20),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x.bottom = element_text(size=26), axis.text.y.left = element_text(size=26),
        axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.y = element_text(size=28), axis.title.x = element_text(size=28))
gg + ggtitle("West region - Värmland. stress = 0.16, 3 dimensions, Bray-Curtis\npermanova not significant")

#permanova
tot.abund2 <- tot.abundance %>%  select(all_of(names(.)[1:3]), where(~ is.numeric(.) && 
                                                                       sum(., na.rm = TRUE) != 0))
Site <- df.värm$Site
# NUll hypothesis: Groups do not differ in spread or position in multivariate space.
perm <- adonis2(tot.abund2 ~ Lupin, data = df.värm)
perm

