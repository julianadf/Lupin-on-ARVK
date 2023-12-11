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

# Load the dataframe ----
data <- read.csv2(here("data/Final Datasets", "final.dataframe.csv" ))

# Restructure
data.res <- data %>% 
  dplyr::select(-Lupinus.polyphyllus) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region = as.factor(Region)) 
str(data.res)

# Divide into the treatments
lupin <- data.res %>% 
  filter(Lupin =="Yes") %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(lupin)

no.lupin <- data.res %>% 
  filter(Lupin =="No") %>% 
  mutate(Region= as.factor(Region)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(no.lupin)

# Calculate Hill diversities and extract them

# Hill diversities for plots with lupin
Hill.lupin <- renyi(lupin[,5:ncol(lupin)], hill = TRUE)
shannon.lupin <- Hill.lupin[,4] # for D = 1

# Hill diversities for plots without lupin
Hill.nolupin <- renyi(no.lupin[,5:ncol(no.lupin)], hill = TRUE)
shannon.nolupin <- Hill.nolupin[,4] # for D = 1

# Create dataframe with the effective number of species (q=1)
db.1 <- bind_cols(lupin[,1:4], shannon.lupin)
db.1 <- rename(db.1, shannon = ...5)
db.2 <- bind_cols(no.lupin[,1:4], shannon.nolupin)
db.2 <- rename(db.2, shannon = ...5)

db <- bind_rows(db.1,db.2)
hist(db$shannon, main="New data", xlab="Hill Diversity, q=1 (Shannon)")

plot(density(db$shannon))#, xlim=c(-100,100), type="b")
plot(density(log(db$shannon)))
shapiro.test(db$shannon) #can assume normality
#write.csv2(db, "db.shannon.csv")

boxplot <- ggplot(db, aes(x=Region, y=shannon, color=Lupin)) + geom_boxplot()
boxplot


# Run analysis
model.shannon <- glmmTMB(shannon ~ Lupin * Region , family= "gaussian" (link="log"), data=db)
# About the warnings: This warning occurs when the optimizer visits a region of parameter space that
# is invalid. It is not a problem as long as the optimizer has left that region of parameter space
# upon convergence, which is indicated by an absence of the model convergence warnings described above.
summary(model.shannon)
car::Anova(model.shannon, type= "III")
AICc(model.shannon)
mod_dharma1 <- model.shannon %>% simulateResiduals(n=1000)
plot(mod_dharma1)
plotResiduals(model.shannon, rank = TRUE, quantreg = FALSE)
plot1 <- visreg(model.shannon, "Lupin", scale="response", by= "Region", rug=FALSE, line = list(col="black"), 
       xlab= "Lupine present", ylab="Effective number of species / plot", gg=TRUE) + theme_bw() +
  theme(text = element_text(size = 16))
# to switch the order of the panels
plot1 + facet_grid(~fct_rev(Region))#, levels = c("West region", "East region"))


# Without interaction
model.shannon2 <- glmmTMB(shannon ~ Lupin +  Region + (1|Site),family= "gaussian" (link="log"), data=db)
AICc(model.shannon2)

# Null model
model.shannon3 <- glmmTMB(shannon ~ 1 + (1|Site), family= "gaussian" (link="log"), data=db)
AICc(model.shannon3)

# Load the balanced dataframe ----
balanced <- read.csv2(here("data/Final Datasets", "balanced.df.final.csv"))

# Restructure
balanced.data <- balanced %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  mutate(County=as.factor(County))  
str(balanced.data)


# Divide into the treatments
lupin <- balanced.data %>% 
  filter(Lupin =="Yes") %>% 
  mutate(County= as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(lupin)

no.lupin <- balanced.data %>% 
  filter(Lupin =="No") %>% 
  mutate(County= as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin))
str(no.lupin)

# Calculate Hill diversities and extract them

# Hill diversities for plots with lupin
Hill.lupin <- renyi(lupin[,5:ncol(lupin)], hill = TRUE)
shannon.lupin <- Hill.lupin[,4] # for D = 1

# Hill diversities for plots without lupin
Hill.nolupin <- renyi(no.lupin[,5:ncol(no.lupin)], hill = TRUE)
shannon.nolupin <- Hill.nolupin[,4] # for D = 1

# Create dataframe with the effective number of species (q=1)
db.1 <- bind_cols(lupin[,1:4], shannon.lupin)
db.1 <- rename(db.1, shannon = ...5)
db.2 <- bind_cols(no.lupin[,1:4], shannon.nolupin)
db.2 <- rename(db.2, shannon = ...5)

db.balanced <- bind_rows(db.1,db.2)
hist(db.balanced$shannon, main="Balanced data", xlab="Hill Diversity, q=1 (Shannon)")
boxplot.balanced <- ggplot(db.balanced, aes(x=County, y=shannon, color=Lupin)) + geom_boxplot()
boxplot.balanced

# Run analysis
model.shannon <- glmmTMB(shannon ~ Lupin * County + (1|Site), family= "Gamma", data=db.balanced)
summary(model.shannon)
car::Anova(model.shannon, type= "III")
AICc(model.shannon)
mod_dharma1 <- model.shannon %>% simulateResiduals(n=1000)
plot(mod_dharma1)
plotResiduals(model.shannon, rank = TRUE, quantreg = FALSE)
visreg(model.shannon, "Lupin", scale="response", by= "County", rug=FALSE, line = list(col="black"), 
       xlab= "Lupin present", ylab="Effective number of species (shannon) / plot") 


model.shannon2 <- glmmTMB(shannon ~ Lupin +  County + (1|Site), family= "Gamma", data=db.balanced)
AICc(model.shannon2)

# Null model
model.shannon3 <- glmmTMB(shannon ~ 1 + (1|Site), family= "Gamma", data=db.balanced)
AICc(model.shannon3)
