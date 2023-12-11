rm(list=ls())
library(here)
library(tidyverse)
library(vegan)
library(permute)
library(geosphere)
library(ggplot2)
library(png)
library(grid)
library(cowplot)
library(magick)
library(glmmTMB)
library(DHARMa)
library(visreg)
library(MuMIn)
#library(diverse)

# Load data
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

db <- df.res[, c(5:ncol(df.res))]

# Calculate Shannon evenness index for plants in each site:

db.H <- diversity(db, index="shannon")
db.S <- db.H / log(specnumber(db))
rich <- specnumber(db)
res <- df.res[, c(1,2,4)]

db.even <- cbind(res, db.S)


# Test for correlations between species richness and evenness
cor.test(rich, db.S, method=c("pearson")) # not correlated

model.even <- glmmTMB(db.S ~ Lupin * Region + (1|Site), beta_family(), data=db.even)
summary(model.even)
car::Anova(model.even, type= "III")
AICc(model.even)
mod_dharma1 <- model.even %>% simulateResiduals(n=1000)
plot(mod_dharma1)
plotResiduals(model.even, rank = TRUE, quantreg = FALSE)
visreg(model.even, "Lupin", scale="response", by= "Region", rug=FALSE, line = list(col="black"), 
       xlab= "Lupin present", ylab="Evenness", gg=TRUE) + theme_bw()

# Uppland
upp <- df.res %>% 
  filter(Region == "East region")
str(upp)

upp.H <- diversity(upp[,5:ncol(upp)], index="shannon")
upp.S <- upp.H / log(specnumber(upp[, 5:ncol(upp)]))
res.upp <- upp[, c(1,2,4)]

rich.upp <- specnumber(upp[, 5:ncol(upp)])

db.upp <- cbind(res.upp, upp.S)


# Test for correlations between species richness and evenness
cor.test(rich.upp, upp.S, method=c("pearson")) # not correlated
ggplot(db.upp, aes(x = rich.upp, y=upp.S)) + geom_point() 

model.upp <- glmmTMB(upp.S ~ Lupin  + (1|Site), beta_family(), data=db.upp)
summary(model.upp )
car::Anova(model.upp , type= "III")
AICc(model.upp )
mod_dharma1 <- model.upp  %>% simulateResiduals(n=1000)
plot(mod_dharma1)
plotResiduals(model.upp , rank = TRUE, quantreg = FALSE)
visreg(model.upp, "Lupin", scale="response",  rug=FALSE, line = list(col="black"), 
       xlab= "Lupin present", ylab="Evenness", gg=TRUE) + theme_bw()



# Värmland
värm <- df.res %>% 
  filter(Region == "West region")
str(värm)

värm.H <- diversity(värm[,5:ncol(värm)], index="shannon")
värm.S <- värm.H / log(specnumber(värm[, 5:ncol(värm)]))
res.värm <- värm[, c(1,2,4)]

rich.värm <- specnumber(värm[, 5:ncol(värm)])

db.värm <- cbind(res.värm, värm.S)


# Test for correlations between species richness and evenness
cor.test(rich.värm, värm.S, method=c("pearson")) # not correlated
ggplot(db.värm, aes(x = rich.värm, y=värm.S)) + geom_point() 

model.värm <- glmmTMB(värm.S ~ Lupin  + (1|Site), beta_family() , data=db.värm)
summary(model.värm)
car::Anova(model.värm , type= "III")
AICc(model.värm)
mod_dharma1 <- model.värm  %>% simulateResiduals(n=1000)
plot(mod_dharma1)
plotResiduals(model.värm , rank = TRUE, quantreg = FALSE)
visreg(model.värm, "Lupin", scale="response",  rug=FALSE, line = list(col="black"), 
       xlab= "Lupin present", ylab="Evenness", gg=TRUE) + theme_bw()

# Completely ignore rare species: Berger-Parker
berger.parker <- diversity(df.res, type="berger-parker", category_row = TRUE)
