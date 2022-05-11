rm(list=ls())
library(here)
library(tidyverse)

# Load data
Tylers.lista <- read.csv2(here("data", "Tylers lista.csv" ))
raw.data <- read.csv2(here("data", "Raw data.env.csv" ))

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
Tyler.database <- step3[,c(1:11, 16, 33:50)]

