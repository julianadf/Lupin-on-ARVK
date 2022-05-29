rm(list=ls())
library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)

# Load data
Tylers.lista <- read.csv2(here("data", "Tylers lista.csv" ))
raw.data <- read.csv2(here("data", "Raw data.env.csv" ))
species.list <- read.csv2(here("data", "Species list.csv" ))

# Extract what I need from the list
# rename columns to match
Tylers <- Tylers.lista %>% 
  rename(Swedish.name = Svenskt.namn) %>% 
  mutate(Swedish.name = sub("(.)", "\\U\\1", Swedish.name, perl=TRUE)) %>% 
  rename(Species = Scientific.name)
head(Tylers)

step1 <- left_join(raw.data, Tylers, by=c("Swedish.name"))
# identify the problematic species
step2 <- step1 %>%
  filter(is.na(Species.y)) %>% 
  group_by(Swedish.name, Species.x) %>% 
  summarise()
step2

# Change the names of those present in both databases so they match
raw.updated <-raw.data %>% 
  # mutate(Swedish.name = dplyr::recode(Swedish.name, "Björnbär"= "") %>% finns inte i Tylers lista
  # mutate(Swedish.name = dplyr::recode(Swedish.name, "Daggkåpa"= "")) %>% finns inte i Tylers lista
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Fläcknycklar"= "Jungfru Marie nycklar")) %>% 
  # mutate(Swedish.name = dplyr::recode(Swedish.name, "Fräken"= "")) %>% 
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Förgätmigej"= "Äkta förgätmigej")) %>% 
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Gråbinka"= "Vanlig gråbinka")) %>% 
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Gråfibbla"= "Vanlig gråfibbla")) %>% 
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Gråögontröst"= "Grå ögontröst")) %>% 
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Hundkex"= "Hundkäx")) %>% 
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Hårig fibbla"= "Hagfibblor")) %>%
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Johannesört"= "Äkta johannesört")) %>%
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Litenblåklocka"= "Liten blåklocka")) %>%
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Lupin"= "Blomsterlupin")) %>%
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Lönn"= "Skogslönn")) %>%
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Majsmörblomma"= "Majsmörblommor")) %>%
  #mutate(Swedish.name = dplyr::recode(Swedish.name, "Maskros"= "")) %>% finns inte i Tylers lista
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Mjölkört"= "Mjölke")) %>%
  # mutate(Swedish.name = dplyr::recode(Swedish.name, "Okänd klint"= "")) %>%
  # mutate(Swedish.name = dplyr::recode(Swedish.name, "Okänd ört"= "")) %>%
  # mutate(Swedish.name = dplyr::recode(Swedish.name, "Okänt gräs"= "")) %>%
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Rallarros"= "Mjölke")) %>% 
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Skogsnäva"= "Midsommarblomster")) %>%
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Stor getväppling"= "Getväppling")) %>%
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Storblåklocka"= "Stor blåklocka")) %>%
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Tussilago"= "Hästhov")) %>%
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Vanligsmörblomma"= "Smörblomma")) %>%
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Vildmorot"= "Morot")) %>%
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Vårbrodd"= "Sydvårbrodd")) %>%
  mutate(Swedish.name = dplyr::recode(Swedish.name, "Älgört"= "Älggräs"))
raw.updated

# Try joining again        
step3 <- left_join(raw.updated, Tylers, by=c("Swedish.name"))
# see if it worked
# identify the problematic species
step4 <- step3 %>%
  filter(is.na(Species.y)) %>% 
  group_by(Swedish.name, Species.x) %>% 
  summarise()
step4

# Extract the data I need and add the species heights
Tyler.database <- step3[,c(1:11, 16, 33:50)] %>% 
  rename(Species = Species.x)

step4 <- left_join(Tyler.database, species.list, by=c("Species"))

step5 <- step4 %>%
  filter(is.na(Swedish.name.y)) %>%  
  group_by(Swedish.name.x, Species) %>% 
  summarise()
step5 # looks good!

entire.database<- step4[,c(1:6,8, 9:13, 20:30, 36:37)]

# filter away pairs based on previous selection and then remove rows with NAs and trees
selected.pairs <- read.csv2(here("data", "selected pairs per site.csv"))

# Use the selected pairs to subset the entire dataset
final.selection <- semi_join(entire.database, selected.pairs, by= c("Site", "parnr"))

db <- final.selection %>% 
  filter(Träd==0) %>% 
  filter(!is.na(Medelhöjd..cm.)) %>% 
  filter(!is.na(Continentality)) %>% 
  rename(Swedish.name = Swedish.name.x) %>% 
  mutate(Site = as.factor(Site)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(Species = as.factor(Species))
db
  
# Calculate the community-weighted species height / traits
# with lupin
db.communityweighted <- db %>% 
  group_by(Site, parnr, Lupin) %>% 
  summarise(
    Weighted.height = weighted.mean(Medelhöjd..cm., Occurrence),
    Weighted.light = weighted.mean(Light, Occurrence),
    Weighted.moist = weighted.mean(Moisture, Occurrence),
    Weighted.soil = weighted.mean(Soil.reaction..pH., Occurrence),
    Weighted.nitrogen = weighted.mean(Nitrogen..N., Occurrence)
  )
db.communityweighted

# run paired t-test and make figures
lupin <- db.communityweighted %>% 
  filter(Lupin == "Yes")

no.lupin <- db.communityweighted %>% 
  filter(Lupin == "No")

t.test(lupin$Weighted.height, no.lupin$Weighted.height, paired = TRUE)

comp.lupin <- list(c("Yes", "No"))

# Community-weighted vegetation height
gg.communityheight <-ggplot(data= db.communityweighted, aes( x=Lupin, y=Weighted.height, fill= Lupin)) + geom_violin(trim=FALSE) + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Community-weighted vegetation height (cm)") + scale_fill_manual(values=c("#f1a340", "#998ec3")) + theme_classic()
gg.communityheight2 <- gg.communityheight + stat_compare_means(comparisons = comp.lupin, paired = TRUE, method = "t.test",
                                                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                                  symbols = c("****", "***", "**", "*", "ns")))
gg.communityheight2

vegetation <- ggplot(data) + geom_dotplot(aes(Veghgt.avg))

# other traits
tidy.db <- db.communityweighted %>% 
  gather(key = "trait", value = "value", Weighted.height:Weighted.nitrogen) %>% 
  filter(trait != "Weighted.height")
tidy.db

comp.traits <- list(c("Yes", "No"))

comp.traits <- list(c("Weighted.light"),c("Weighted.nitrogen"), c("Weighted.soil"),c("Weighted.moist"))

stat.test <- tidy.db %>%
  group_by(trait) %>%
  pairwise_t_test(
    value ~ Lupin, paired = TRUE, 
    p.adjust.method = "bonferroni") %>% 
    select(-df, -statistic, -p) # Remove details
stat.test

gg.communitytraits <-ggplot(data= tidy.db, aes(x=trait, y=value, color=Lupin)) + 
  geom_violin(trim=FALSE, aes(fill = Lupin)) +
  geom_boxplot(aes(color = Lupin), position=position_dodge(width=0.9), width=0.2,  show.legend = FALSE) + 
  scale_color_manual(values=c("black", "black"))+
  scale_fill_manual(values=c("#f1a340", "#998ec3")) + theme_classic() +
  scale_x_discrete(labels=c(c("Light", "Moisture", "Soil reaction", "Nitrogen"))) +
  theme(axis.title.x = element_blank()) 
  #ylab("Ecological value index")
gg.communitytraits

stat.test <- stat.test %>% add_xy_position(x = "trait")

gg.communitytraits + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08, hide.ns = TRUE, tip.length = 0
)   

