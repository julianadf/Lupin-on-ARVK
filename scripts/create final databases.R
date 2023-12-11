rm(list=ls())
library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)

# Create database without trees and with the new sites
db <- read.csv2(here("data/Final Datasets", "complete.db.nograss.ENV.csv" ))

# Filter away trees and tidy
db.long <- db %>% 
  filter(!Träd==1) %>% 
  mutate(County = dplyr::recode(County, "Örebro"= "Värmland")) %>%
  mutate(County = as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  rename(pH = Soil.reaction..pH.) %>%
  rename(Nitrogen = Nitrogen..N.) %>% 
  rename(Sp.length = Medelhöjd..cm.) %>% 
  group_by(County, Site, parnr, Lupin, Species, Veghgt.avg, Litt.avg, 
           Light, Moisture, pH, Nitrogen, Family, Sp.length, Management.dependent) %>%
  distinct(Species, Occ, .keep_all = TRUE) %>% 
  summarise(Occurrences= sum(Occ)) 
db.long

# Load coordinate data
coords <- read.csv2(here("data/OLD sites", "sites.spatial.csv"))

coords <- coords %>% 
  mutate(Site=as.factor(Site)) %>% 
  mutate(parnr=as.factor(parnr)) %>% 
  mutate(Lupin=as.factor(Lupin)) %>% 
  mutate(N=as.numeric(N)) %>% 
  mutate(E=as.numeric(E)) %>% 
  mutate(lat=as.numeric(lat)) %>% 
  mutate(lng=as.numeric(lng))
str(coords)

db.tidy <- left_join(db.long, coords, by=c("Site", "Lupin", "parnr"))

# Turn into wide
species <- db.tidy %>% 
  mutate(County = dplyr::recode(County, "Örebro"= "Värmland")) %>%
  mutate(County = as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  group_by(County, Site, parnr, Lupin, Occurrences, Species) %>% 
  summarise()
str(species)

df <- species %>% 
  spread(., key="Species", value = "Occurrences") %>% 
  replace(is.na(.), 0)
df


# Load selected sites and filter
selected.sites <- read.csv2(here("data", "Urval Uppland_LowFreqLP.csv"))
selected.sites <- selected.sites %>%
  mutate(Site=as.factor(Site)) %>% 
  rename(County = Region) %>% 
  mutate(County = as.factor(County))

sites.tidy <- semi_join(db.tidy, selected.sites, by=c("Site"))

# Subset sites and pairs
# Load selected pairs and filter
pairs <- read.csv2(here("data", "selected pairs per site.csv"))
pairs <- pairs %>%  mutate(Site = as.factor(Site)) %>% mutate(parnr = as.factor(parnr))
final.tidy <- semi_join(sites.tidy, pairs, by=c("Site", "parnr"))

# Final tidy
#write.csv2(final.tidy, "final.tidy.csv")


# Turn into wide
df.final <- final.tidy %>% 
  mutate(County = dplyr::recode(County, "Örebro"= "Värmland")) %>%
  mutate(County = as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  group_by(County, Site, parnr, Lupin, Occurrences, Species) %>% 
  summarise()
str(df.final)

df <- df.final %>% 
  spread(., key="Species", value = "Occurrences") %>% 
  replace(is.na(.), 0)
df
# Final dataframe
#write.csv2(df, "final.dataframe.csv")

# Create species richness database
data <- final.tidy %>% 
  filter(Species != "Lupinus polyphyllus") %>% 
  mutate(County = as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  group_by(County, Site, Lupin) %>% 
  summarise(rich = n())
str(data)

write.csv2(data, "species.richness.noLP.csv")


# Create Tommy plot
data <- read.csv2(here("data/Final Datasets", "species.richness.noLP.csv"))

Tommy.plot <- ggplot(data, aes(x=fct_inorder(Site), y=rich, color= Lupin)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank()) +
  ggtitle("New sites & no trees & low frequency LP in no.Lupin plots")
Tommy.plot

# Figure S1
df.final <- read.csv2(here("data/Final Datasets", "final.tidy.csv"))

tidy.nolupin <- df.final %>% 
  filter(Species != "Lupinus polyphyllus") %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>%
  rename(Region = County) %>% 
  mutate(Lupin = as.factor(Lupin))

Fig.S1 <- ggplot(tidy.nolupin, aes(x=fct_infreq(Species), y = Occurrences, fill=Region)) + geom_col() +
  facet_wrap(vars(Lupin), nrow =2) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(), 
        panel.grid.major = element_blank()) +
  ylab(bquote('Total occurrences')) + font("ylab", size=11) +
  scale_x_discrete()+ scale_y_continuous(expand = c(0,0)) + ylim(0, 120) +
  scale_fill_manual("Region", values = c("#d7191c", "#2c7bb6"))

Fig.S1  + theme(legend.position = c(0.9,0.87), legend.title=element_blank())

unique(df.final$Site)

# % minskning av varje art / region
# per site
No.Lupin <- tidy.nolupin %>% 
  filter(Lupin == "No")
No.Lupin

Lupin.No <- No.Lupin[, c(1,2,5,15)]

Yes.Lupin <- tidy.nolupin %>% 
  filter(Lupin == "Yes")
Yes.Lupin

Lupin.Yes <- Yes.Lupin[, c(1,2,5,15)]

Lupin.all <- full_join(Lupin.No, Lupin.Yes, by=c("County", "Site", "Species")) %>% 
  rename(OccNo=Occurrences.x) %>% 
  rename(OccYes=Occurrences.y) %>% 
  replace(is.na(.), 0) %>% 
  mutate(Difference = OccNo - OccYes) %>%  
  mutate(percent = -(Difference / OccNo)*100) %>% 
  filter_all(all_vars(!is.infinite(.))) %>% 
  group_by(County, Species) %>% 
  summarise(mean.perc = mean(percent, na.rm=T))
Lupin.all


# plot
gg.diff <- ggplot(Lupin.all, aes(x=fct_infreq(Species), y=mean.perc, fill=County))+ geom_col() +
  facet_wrap(vars(County), nrow =2) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_blank(), 
        panel.grid.major = element_blank()) +
  ylab(bquote('Average percent change: without LP - with LP')) + font("ylab", size=11)
gg.diff


# Check frequency of lupins in lupin-plots ----
data <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))

lupin <- data %>% 
  filter(Lupin =="Yes") 
str(lupin)

no.lupin <- data %>% 
  filter(Lupin =="No")
str(no.lupin)

freq <- lupin %>% 
  mutate(lupin.freq =Lupinus.polyphyllus /16)
freq

lupin.freq <- freq[,c(1:2,ncol(freq))] 
hist(lupin.freq$lupin.freq)

uppsala <- lupin.freq %>% filter(County =="Uppsala")
värmland <- lupin.freq %>% filter(County !="Uppsala") 

t.test(uppsala$lupin.freq, värmland$lupin.freq)

gg.county <-ggplot(data=lupin.freq, aes(x=County, y=lupin.freq)) + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Medelfrekvens lupin / yta")
gg.county + ggtitle("plots with Lupin")

# In no-lupin plots
freq.no <- no.lupin %>% 
  mutate(lupin.freq = Lupinus.polyphyllus/16)
freq.no

nolupin.freq <- freq.no[,c(1:2,ncol(freq.no))] 
hist(nolupin.freq$lupin.freq)

uppsala.no <- nolupin.freq %>% filter(County =="Uppsala")
värmland.no <- nolupin.freq %>% filter(County !="Uppsala") 

t.test(uppsala.no$lupin.freq, värmland.no$lupin.freq)

gg.county <-ggplot(data=nolupin.freq, aes(x=County, y=lupin.freq)) + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Medelfrekvens lupin / yta")
gg.county + ggtitle("Plots without lupin")


# Unbalanced ----

# Create database without trees and with the new sites
unbalanced <- read.csv2(here("data/Final Datasets", "complete.db.nograss.ENV.csv" ))

# Filter away trees and tidy
db.unbal <- unbalanced %>% 
  filter(!Träd==1) %>% 
  mutate(County = dplyr::recode(County, "Örebro"= "Värmland")) %>%
  mutate(County = as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  rename(pH = Soil.reaction..pH.) %>%
  rename(Nitrogen = Nitrogen..N.) %>% 
  rename(Sp.length = Medelhöjd..cm.) %>% 
  group_by(County, Site, parnr, Lupin, Species, Veghgt.avg, Litt.avg, 
           Light, Moisture, pH, Nitrogen, Family, Sp.length, Management.dependent) %>%
  distinct(Species, Occ, .keep_all = TRUE) %>% 
  summarise(Occurrences= sum(Occ)) 
db.unbal

# Select pairs
selected.sites <- read.csv2(here("data", "selected pairs per site.csv"))

selected.sites <- selected.sites %>%
  mutate(Site=as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(County = dplyr::recode(County, "Örebro"= "Värmland")) %>%
  mutate(County = as.factor(County))
str(selected.sites)

unbalanced.tidy <- semi_join(db.unbal, selected.sites, by=c("Site", "parnr"))


# Check that there are 43 pairs
nobs(unique(unbalanced.tidy[,2])) # There are 43 sites
nobs(unique(unbalanced.tidy[,c(2:4)])) # and 86 pairs
nobs(unique(unbalanced.tidy[,5])) # and 145 species 


# Database for species richness
# Create species richness database
rich <- unbalanced.tidy %>% 
  filter(Species != "Lupinus polyphyllus") %>% 
  mutate(County = as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  group_by(County, Site, Lupin) %>% 
  summarise(rich = n())
str(rich)

# Create boxplot

boxplot <- ggplot(rich, aes(x=County, y=rich, color=Lupin)) + geom_boxplot()
boxplot + ggtitle("Unbalanced data: 31 sites in Uppland and 12 in Värmland")

model.rich <- glmmTMB(rich ~ Lupin * County + (1|Site), family= "poisson", data=rich)
summary(model.rich)
visreg(model.rich, "Lupin", scale="response", by= "County", rug=FALSE, 
       line = list(col="black"), xlab= "Lupin present", ylab="Species richness / plot") 

# Create final unbalanced dataframe

# Turn into wide
df.final <- unbalanced.tidy %>% 
  filter(Species != "Lupinus polyphyllus") %>%
  mutate(County = as.factor(County)) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  group_by(County, Site, parnr, Lupin, Occurrences, Species) %>% 
  summarise()
str(df.final)

df.unbalanced <- df.final %>% 
  spread(., key="Species", value = "Occurrences") %>% 
  replace(is.na(.), 0)
df.unbalanced

#write.csv2(df.unbalanced, "df.unbalanced.csv")

# Check frequency of lupins in lupin-plots 
data <- read.csv2(here("data/Final Datasets", "final.dataframe.csv"))

lupin <- data %>% 
  filter(Lupin =="Yes") 
str(lupin)

no.lupin <- data %>% 
  filter(Lupin =="No")
str(no.lupin)

freq <- lupin %>% 
  mutate(lupin.freq =Lupinus.polyphyllus /16)
freq

lupin.freq <- freq[,c(1:2,ncol(freq))] 
hist(lupin.freq$lupin.freq)

uppsala <- lupin.freq %>% filter(County =="Uppsala")
värmland <- lupin.freq %>% filter(County !="Uppsala") 

t.test(uppsala$lupin.freq, värmland$lupin.freq)

gg.county <-ggplot(data=lupin.freq, aes(x=County, y=lupin.freq)) + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Medelfrekvens lupin / yta")
gg.county + ggtitle("plots with Lupin")

# In no-lupin plots
freq.no <- no.lupin %>% 
  mutate(lupin.freq = Lupinus.polyphyllus/16)
freq.no

nolupin.freq <- freq.no[,c(1:2,ncol(freq.no))] 
hist(nolupin.freq$lupin.freq)

uppsala.no <- nolupin.freq %>% filter(County =="Uppsala")
värmland.no <- nolupin.freq %>% filter(County !="Uppsala") 

t.test(uppsala.no$lupin.freq, värmland.no$lupin.freq)

gg.county <-ggplot(data=nolupin.freq, aes(x=County, y=lupin.freq)) + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Medelfrekvens lupin / yta")
gg.county + ggtitle("Plots without lupin")

