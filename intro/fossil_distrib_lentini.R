library(tidyverse)


raw_datafile <-tibble(fossils = read_lines("~/projects/kakapo-genomics/data/lentini2017_kakapo_fossil_sites.txt", skip = 1))

# Tidy the data into a useable table
fossil_sites <- raw_datafile %>% 
  separate(fossils,
         into = c(1:7, "Site description", "Grid East", "Grid North", "Latitude", "Longitude", "Altitude", "Reference"), 
         sep = " ", 
         fill = "left") %>% 
  unite(col = "site", c(1:7, `Site description`), sep = " ", remove = TRUE, na.rm = TRUE) %>% 
  mutate(across(c(`Grid East`, `Grid North`, Latitude, Longitude, Altitude, Reference), as.numeric)) %>% 
  select(site, Latitude, Longitude) %>% 
  rename(latitude = Latitude, longitude = Longitude) %>% 
  mutate(type = "fossil")


# Historical sites mentioned in Lentini et al. 2018. The lat long was obtained from me looking up the places on NZTOPO, very very approximate
historical_sites_lentini <- tibble(site = c("Esperance River", "Dusky Sound", "Sinbad Gully", "Scollay's Flat", "Pegasus Creek", "Basin Creek"),
                           latitude = c(-44.73,-45.68, -44.64, -47.15, -47.14, -47.16),
                           longitude = c(167.98, 166.99, 167.83, 167.78, 167.70, 167.63)) %>% 
  mutate(type = "Historical sighting")


all_sites <- bind_rows(fossil_sites, historical_sites_lentini) %>% 
  filter(type != "introduction") %>% 
  mutate(type = case_when(type == "fossil" ~ "Subfossil remains",
                          type == "subfossil" ~ "Subfossil remains",
                          TRUE ~ type))
  
# write_delim(all_sites, "~/projects/kakapo-genomics/data/kakapo_sites_lentini_tidied.csv", delim = ",")

