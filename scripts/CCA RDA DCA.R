rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(permute)
library(geosphere)
library(funrar)

df <- read.csv2(here("data/Final Datasets", "final.dataframe.names.csv"))

df.res <- df %>% 
  dplyr::select(-Lup.pol.) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  mutate(County=as.factor(County)) 
str(df.res)

comm <- df.res[, 5:ncol(df.res)]
env.1 <- read.csv2(here("data/Final Datasets", "env.csv"))
env.2 <- df[,c(1,2,4)]
env.3 <-merge(env.1, env.2, by=c("Site", "Lupin"))
env <- env.3[, c(1:5)]

env.geo <- env[,1:3]
env.veg <- env[,4:5]

# Need to transform explanatory vars
hist(env$Veghgt.avg)
hist(env$Litt.avg)
comm.hel <- decostand(comm, method="hellinger")

# DCA ----
cca.plants <- cca(comm)
summary(cca.plants)
summary(cca.plants, scaling = 1)
screeplot(cca.plants, bstick=TRUE)

par(mfrow = c(1,2))
# Scaling 1: sites are centroids of species
plot(cca.plants, scaling = 1, main = "CA plant abundances - biplot scaling 1", type = "text")
# Scaling 2: species are centroids of sites
plot(cca.plants, main = "CA plant abundances - biplot scaling 2",  type = "text")

# A posteriori curve fitting Lupin and Region
plot(cca.plants, main = "CA plant abundances - biplot scaling 2", sub="Veg. height (red), Litt.dep(green)")
spe.ca.env <- envfit(cca.plants ~ Veghgt.avg + Litt.avg, env)
plot(spe.ca.env)
ordisurf(cca.plants, env$Veghgt.avg, add =T)
ordisurf(cca.plants, env$Litt.avg, add =T, col="green")

vegemite(comm, cca.plants, scale = T)

# RDA ----
spe.rda <- rda(comm.hel ~ Lupin + Veghgt.avg + Litt.avg + Region + Site, data=env)
summary(spe.rda)

# Forward selection of variables:
fwd.sel <- ordiR2step(rda(comm.hel ~ 1, data = env), # lower model limit (simple!)
                      scope = formula(spe.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE) # change to TRUE to see the selection process!


# Check the new model with forward-selected variables
fwd.sel$call


# Write our new model
spe.rda.signif <- rda(comm.hel ~ Lupin, data = env)
# check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(spe.rda.signif)

anova.cca(spe.rda, step = 1000)

# Type 1 scaling
ordiplot(spe.rda.signif, scaling = 1, type = "text")
# Type 2 scaling
ordiplot(spe.rda.signif, scaling = 2, type = "text")


# Partial RDA ----
spe.partial.rda <- rda(comm.hel ~ Lupin + Veghgt.avg + Litt.avg + Region + Condition(Site), data=env)
summary(spe.partial.rda)

# Forward selection of variables:
fwd.sel <- ordiR2step(rda(comm.hel ~ 1, data = env), # lower model limit (simple!)
                      scope = formula(spe.rda), # upper model limit (the "full" model)
                      direction = "forward",
                      R2scope = TRUE, # can't surpass the "full" model's R2
                      pstep = 1000,
                      trace = FALSE) # change to TRUE to see the selection process!


# Check the new model with forward-selected variables
fwd.sel$call


# Write our new model
spe.rda.signif <- rda(comm.hel ~ Lupin + Region + Veghgt.avg + Litt.avg, data = env)
# check the adjusted R2 (corrected for the number of explanatory variables)
RsquareAdj(spe.rda.signif)

anova.cca(spe.rda.signif, step = 1000)

# Type 1 scaling
ordiplot(spe.rda.signif, scaling = 1, type = "text")
# Type 2 scaling
ordiplot(spe.rda.signif, scaling = 2, type = "text")

# Custom triplot code! ----

## extract % explained by the first 2 axes
perc <- round(100*(summary(spe.rda.signif)$cont$importance[2, 1:2]), 2)

## extract scores - these are coordinates in the RDA space
sc_si <- scores(spe.rda.signif, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(spe.rda.signif, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(spe.rda.signif, display="bp", choices=c(1, 2), scaling=1)

## Custom triplot, step by step

# Set up a blank plot with scaling, axes, and labels
plot(spe.rda.signif,
     scaling = 2, # set scaling type 
     type = "none", # this excludes the plotting of any points from the results
     frame = FALSE,
     # set axis limits
     xlim = c(-1,1), 
     ylim = c(-1,1),
     # label the plot (title, and axes)
     main = "Triplot RDA - scaling 2",
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)") 
)
# add points for site scores
points(sc_si, 
       pch = 21, # set shape (here, circle with a fill colour)
       col = "black", # outline colour
       bg = "steelblue", # fill colour
       cex = 1.2) # size
# add points for species scores
points(sc_sp, 
       pch = 22, # set shape (here, square with a fill colour)
       col = "black",
       bg = "#f2bd33", 
       cex = 1.2)
# add text labels for species abbreviations
text(sc_sp + c(0.03, 0.09), # adjust text coordinates to avoid overlap with points 
     labels = rownames(sc_sp), 
     col = "grey40", 
     font = 2, # bold
     cex = 0.6)
# add arrows for effects of the expanatory variables
arrows(0,0, # start them from (0,0)
       sc_bp[,1], sc_bp[,2], # end them at the score value
       col = "red", 
       lwd = 3)
# add text labels for arrows
text(x = sc_bp[,1] -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2] - 0.03, 
     labels = rownames(sc_bp), 
     col = "red", 
     cex = 1, 
     font = 2)

