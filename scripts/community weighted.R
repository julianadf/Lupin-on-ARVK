rm(list=ls())
library(here)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(RColorBrewer)

# Load data
tidy <- read.csv2(here("data/Final Datasets", "final.tidy.csv"))

# Turn into wide
tidy.res <- tidy %>% 
  filter(!(Species %in% c("Lupinus polyphyllus", "Rubus fruticosus"))) %>%
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region =as.factor(Region))
str(tidy.res)

df <- tidy.res %>% 
  group_by(Region, Site, Lupin) %>% 
  summarise(Veghgt.avg, Litt.avg) %>% 
  distinct()
df

# Check that there are 43 pairs
nobs(unique(tidy.res[,2])) # There are 24 sites
nobs(unique(tidy.res[,c(2:4)])) # and 48 pairs
nobs(unique(tidy.res[,5])) # and 127 species without lupin & björnbär

  # Calculate the community-weighted species height / traits ----
# Remove NAs
data <- tidy.res %>%
  filter(!is.na(Light)) 
data
  

# Community weighted means - all data
db.communityweighted <- data %>% 
  filter(!is.na(Sp.length)) %>% 
  group_by(Region, Site, parnr, Lupin) %>% 
  summarise(
     Weighted.height = weighted.mean(Sp.length, Occurrences),
     Weighted.light = weighted.mean(Light, Occurrences),
     Weighted.moist = weighted.mean(Moisture, Occurrences),
     Weighted.soil = weighted.mean(pH, Occurrences),
     Weighted.nitrogen = weighted.mean(Nitrogen, Occurrences)
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
gg.communityheight <-ggplot(data= db.communityweighted, aes(x=Lupin, y=Weighted.height, fill= Lupin)) + geom_violin(trim=FALSE) + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Community-weighted vegetation height (cm)") + xlab("Lupine") +
  scale_fill_manual(values=alpha(c("#008837", "#7b3294"),0.6)) + theme_classic() +
  theme(legend.position = "none")
gg.communityheight2 <- gg.communityheight + stat_compare_means(comparisons = comp.lupin, paired = TRUE, method = "t.test",
                                                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                                  symbols = c("****", "***", "**", "*", "ns")))
gg.communityheight2

vegetation <- ggplot(data) + geom_dotplot(aes(Veghgt.avg))

# other traits
tidy.db <- db.communityweighted %>% 
  rename(Lupine = Lupin) %>% 
  gather(key = "trait", value = "value", Weighted.height:Weighted.nitrogen) %>% 
  filter(trait != "Weighted.height")
tidy.db

comp.traits <- list(c("Yes", "No"))

comp.traits <- list(c("Weighted.light"),c("Weighted.nitrogen"), c("Weighted.soil"),c("Weighted.moist"))

stat.test <- tidy.db %>%
  group_by(trait) %>%
  pairwise_t_test(
    value ~ Lupine, paired = TRUE, 
    p.adjust.method = "bonferroni")  
#select(-df, -statistic, -p) # Remove details
stat.test

gg.communitytraits <-ggplot(data= tidy.db, aes(x=trait, y=value, color=Lupine)) + 
  geom_violin(trim=FALSE, aes(fill = Lupine)) +
  geom_boxplot(aes(color = Lupine), position=position_dodge(width=0.9), width=0.2,  show.legend = FALSE) + 
  scale_color_manual(values=c("black", "black"))+
  scale_fill_manual( values=alpha(c("#008837", "#7b3294"),0.6)) + 
  theme_classic() +
  scale_x_discrete(labels=c(c("Light", "Moisture",  "Nitrogen", "Soil reaction"))) +
  theme(axis.title.x = element_blank(), text = element_text(size=16)) +
  ylab("Ecological value index")
gg.communitytraits 

stat.test <- stat.test %>% add_xy_position(x = "trait")

gg.communitytraits + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08, hide.ns = FALSE, tip.length = 0
)   


# Look at weighted veg. height between regions
# run paired t-test and make figures
db.communityweighted <- db.communityweighted %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  mutate(County = dplyr::recode(County, "Örebro"= "West region")) 
db.communityweighted

Uppsala <- db.communityweighted %>% 
  filter(Region == "East region")

Värmland <- db.communityweighted %>% 
  filter(Region == "West region")

t.test(Uppsala$Weighted.height, Värmland$Weighted.height, paired = TRUE)

comp.region <- list(c("East region", "West region"))

# Community-weighted vegetation height
gg.communityheight <-ggplot(data= db.communityweighted, aes(x=Region, y=Weighted.height, fill= Region)) + geom_violin(trim=FALSE) + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Community-weighted vegetation height (cm)") + scale_fill_manual(values=alpha(c("#d7191c", "#2c7bb6"),0.7)) + theme_classic()
gg.communityheight2 <- gg.communityheight + stat_compare_means(comparisons = comp.region, paired = FALSE, method = "t.test",
                                                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                                  symbols = c("****", "***", "**", "*", "ns")))
gg.communityheight2

# Uppland ----
tidy <- read.csv2(here("data/Final Datasets", "final.tidy.csv"))

tidy.res <- tidy %>%
  filter(!(Species %in% c("Lupinus polyphyllus", "Rubus fruticosus"))) %>%
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region =as.factor(Region))  
str(tidy.res)

uppsala <- tidy.res %>% 
  filter(Region == "East region")
uppsala

# Remove NAs
data <- uppsala %>%
  filter(!is.na(Light)) %>% 
  filter(!is.na(Sp.length))
data

# with lupin
db.communityweighted <- data %>% 
  group_by(Site, parnr, Lupin) %>% 
  summarise(
    Weighted.height = weighted.mean(Sp.length, Occurrences),
    Weighted.light = weighted.mean(Light, Occurrences),
    Weighted.moist = weighted.mean(Moisture, Occurrences),
    Weighted.soil = weighted.mean(pH, Occurrences),
    Weighted.nitrogen = weighted.mean(Nitrogen, Occurrences)
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
gg.communityheight <-ggplot(data= db.communityweighted, aes(x=Lupin, y=Weighted.height, fill= Lupin)) + geom_violin(trim=FALSE) + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Community-weighted vegetation height (cm)") + scale_fill_manual(values=alpha(c("#008837", "#7b3294"),0.7)) + theme_classic()
gg.communityheight2 <- gg.communityheight + stat_compare_means(comparisons = comp.lupin, paired = TRUE, method = "t.test",
                                                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                                  symbols = c("****", "***", "**", "*", "ns")))
gg.communityheight2 + ggtitle("Uppsala")

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
    p.adjust.method = "bonferroni")  
#select(-df, -statistic, -p) # Remove details
stat.test

gg.communitytraits <-ggplot(data= tidy.db, aes(x=trait, y=value, color=Lupin)) + 
  geom_violin(trim=FALSE, aes(fill = Lupin)) +
  geom_boxplot(aes(color = Lupin), position=position_dodge(width=0.9), width=0.2,  show.legend = FALSE) + 
  scale_color_manual(values=c("black", "black"))+
  scale_fill_manual(values=alpha(c("#008837", "#7b3294"),0.7), name= "Lupine") + theme_classic() +
  scale_x_discrete(labels=c(c("Light", "Moisture", "Nitrogen", "Soil reaction"))) +
  theme(axis.title.x = element_blank()) +
  ylab("Ecological value index")
gg.communitytraits #+ ggtitle("Uppsala")

stat.test <- stat.test %>% add_xy_position(x = "trait")

gg.communitytraits + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08, hide.ns = FALSE, tip.length = 0
)   

# Värmland ----
tidy <- read.csv2(here("data/Final Datasets", "final.tidy.csv"))

tidy.res <- tidy %>% 
  filter(!(Species %in% c("Lupinus polyphyllus", "Rubus fruticosus"))) %>%
  mutate(Site = as.factor(Site)) %>% 
  mutate(parnr = as.factor(parnr)) %>% 
  mutate(Lupin = as.factor(Lupin)) %>% 
  mutate(County = dplyr::recode(County, "Uppsala"= "East region")) %>% 
  mutate(County = dplyr::recode(County, "Värmland"= "West region")) %>% 
  rename(Region = County) %>% 
  mutate(Region =as.factor(Region))  
str(tidy.res)

Värmland <- tidy.res %>% 
  filter(Region == "West region")
Värmland

# Remove NAs
data <- Värmland %>%
  filter(!is.na(Light)) %>% 
  filter(!is.na(Sp.length))
data

# with lupin
db.communityweighted <- data %>% 
  group_by(Site, parnr, Lupin) %>% 
  summarise(
    Weighted.height = weighted.mean(Sp.length, Occurrences),
    Weighted.light = weighted.mean(Light, Occurrences),
    Weighted.moist = weighted.mean(Moisture, Occurrences),
    Weighted.soil = weighted.mean(pH, Occurrences),
    Weighted.nitrogen = weighted.mean(Nitrogen, Occurrences)
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
gg.communityheight <-ggplot(data= db.communityweighted, aes(x=Lupin, y=Weighted.height, fill= Lupin)) + geom_violin(trim=FALSE) + geom_boxplot(width = 0.1, fill = "white")+
  ylab("Community-weighted vegetation height (cm)") + scale_fill_manual(values=alpha(c("#008837", "#7b3294"),0.7)) + theme_classic()
gg.communityheight2 <- gg.communityheight + stat_compare_means(comparisons = comp.lupin, paired = TRUE, method = "t.test",
                                                               symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                                                  symbols = c("****", "***", "**", "*", "ns")))
gg.communityheight2 + ggtitle("Värmland")

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
    p.adjust.method = "bonferroni")  
#select(-df, -statistic, -p) # Remove details
stat.test

gg.communitytraits <-ggplot(data= tidy.db, aes(x=trait, y=value, color=Lupin)) + 
  geom_violin(trim=FALSE, aes(fill = Lupin)) +
  geom_boxplot(aes(color = Lupin), position=position_dodge(width=0.9), width=0.2,  show.legend = FALSE) + 
  scale_color_manual(values=c("black", "black"))+
  scale_fill_manual(values=alpha(c("#008837", "#7b3294"),0.7), name="Lupine") + theme_classic() +
  scale_x_discrete(labels=c(c("Light", "Moisture", "Nitrogen", "Soil reaction"))) +
  theme(axis.title.x = element_blank()) +
  ylab("Ecological value index")
gg.communitytraits

stat.test <- stat.test %>% add_xy_position(x = "trait")

gg.communitytraits + stat_pvalue_manual(
  stat.test, label = "p.adj.signif", 
  step.increase = 0.08, hide.ns = FALSE, tip.length = 0
)   


