library(tidyverse)

characters <- read_delim("~/projects/kakapo-genomics/data/phenotypic_data/Characters/kakapo_characters.csv", delim = ",") %>% 
  select(birdName, houseID, sex, alive, demise, dateHatched, dateLaid, discoveryDate)# %>% 
  # filter(alive != TRUE)

timeline_file <- read_delim("~/projects/kakapo-genomics/data/bird_info_timeline.csv", delim = ",") 


# Remove no-data indivs, or birds that died at 2 years old or younger
timeline_res <- timeline_file %>% 
  filter(is.na(hatch_year) == FALSE) %>% 
  mutate(age = death_year - hatch_year) %>% 
  filter(death_year != 2100) %>%
  filter(age > 1) %>% # was previously 2
  group_by(age) %>%
  summarise(n = n()) %>%
  mutate(previ = lag(cumsum(n),default = 0)) %>% 
  mutate(dead_tot = n /130 ) %>%
  mutate(dead_tot_previ = n/(130-previ)) %>% 
  mutate(huh = cumsum(n)/130) %>%   # These two are equivalent
  mutate(csum = cumsum(dead_tot))


ggplot(timeline_res, aes(x = age, y = csum))+
  geom_point() +
  labs(x = "Age", y = "Csum")+
  xlim(0,30)+
  ylim(0,1)+ 
  theme_classic()

ggplot(timeline_res, aes(x = age, y = dead_tot_previ))+
  geom_point() +
  labs(x = "Age", y = "Dead/Total recorded")+
  xlim(0,30)+
  ylim(0,0.25)+ 
  theme_classic()

# Don't think this is actually informative without adjustment
# Prop of total pop that died at that age
# ggplot(timeline_res, aes(x = age, y = dead_tot))+
#   geom_point() +
#   labs(x = "Age", y = "Dead/Total recorded")+
#   xlim(0,60)+
#   ylim(0,0.05)+ 
#   theme_classic()


####################
# Manually updated #
####################

timeline_updated <- timeline_file %>% 
  filter(is.na(hatch_year) == FALSE) %>% 
  mutate(death_year = case_when(ID == "Ihi" ~ 2021,
                                 ID == "Mila" ~ 2021,
                                 ID == "Millie" ~ 2021,
                                TRUE ~ death_year)) %>% 
  mutate(age = death_year - hatch_year) 

tc <- timeline_updated %>% 
  filter(age > 1) %>% # was previously 2
  filter(death_year != 2100) %>% 
  group_by(age) %>%
  summarise(n = n()) %>% 
  mutate(dead_tot = n /120 ) %>% 
  mutate(csum = cumsum(dead_tot))
  
ggplot(tc, aes(x = age, y = csum))+
  geom_point() +
  labs(x = "Age", y = "Proportion of pop deceased")+
  xlim(0,30)+
  ylim(0,1)+ 
  theme_classic()

###########################
# Minimum age of founders #
###########################

founders <- read_delim("~/projects/kakapo-genomics/data/cohort_by_gen.txt", delim = ",") %>% 
  filter(gen == 1) %>% 
  pull(ID)

names_of_interest <- c("John-Girl", "Sarah", "Sandra", "Fuchsia", "Gunner", "Lee", "Sass", "Richard_Henry", "Whiskas", "Waynebo", "Barnard", "Ben", "Lionel", "Smoko", "Jimmy", "Piripi", "Arab", "Gumboots", "Merty", "Felix")

res <- characters %>% 
  filter(birdName %in% c(founders, names_of_interest)) %>% 
  mutate(lay_year = str_extract(dateLaid, "(?<=/)[:digit:]{4}(?=\\s)")) %>% 
  mutate(hatch_year = str_extract(dateHatched, "(?<=/)[:digit:]{4}(?=\\s)")) %>% 
  mutate(death_year = str_extract(demise, "(?<=/)[:digit:]{4}(?=\\s)")) %>% 
  mutate(discovery_year = str_extract(discoveryDate, "(?<=/)[:digit:]{4}(?=\\s)"))


bird_info_rships <- read_delim("~/projects/kakapo-genomics/data/bird_info_rships_nov20.csv", delim = ",") %>% 
  pivot_longer(cols = c(mother, father), values_to = "parent") %>% 
  filter(parent %in% c(founders, names_of_interest)) %>% 
  group_by(parent) %>% 
  summarise(oldest_chick_hatched = min(hatch_year))

founder_approximate <- left_join(res, bird_info_rships, by = c(birdName = "parent")) %>% 
  mutate(min_age_based_offspring = case_when(sex == "Male" ~ (oldest_chick_hatched -7),
                                             sex == "Female" ~ (oldest_chick_hatched - 5)),
         best_guess = case_when(is.na(hatch_year) == FALSE ~ hatch_year,
                                is.na(hatch_year) == TRUE ~ pmin(discovery_year, min_age_based_offspring,na.rm = TRUE))) %>% 
  select(birdName, best_guess, death_year) %>% 
  mutate(age = case_when(is.na(death_year) ~ 2021 - as.numeric(best_guess),
                           TRUE ~ as.numeric(death_year) - as.numeric(best_guess))) %>% 
  mutate(alive = case_when(is.na(death_year) == TRUE ~ "yes",
                           is.na(death_year) == FALSE ~ "no"),
         death_year = as.numeric(death_year), 
         best_guess = as.numeric(best_guess)) %>% 
  rename(ID = birdName,
         hatch_year = best_guess)

fa <- founder_approximate %>% 
  group_by(age, alive) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  filter(alive == "no") %>% 
  mutate(dead_tot = n/43) %>% 
  mutate(csum = cumsum(dead_tot)) %>% 
  mutate(previ = lag(cumsum(n),default = 0)) %>% 
  mutate(dead_tot_previ = n/(43-previ))

ggplot(fa, aes(x = age, y = csum))+
  geom_point() +
  labs(x = "Age", y = "Csum")+
  xlim(0,30)+
  ylim(0,1)+ 
  theme_classic()

# Don't think this is actually informative without adjustment
# ggplot(fa, aes(x = age, y = dead_tot))+
#   geom_point() +
#   labs(x = "Age", y = "Prop")+
#   xlim(0,30)+
#   ylim(0,1)+ 
#   theme_classic()

ggplot(fa, aes(x = age, y = dead_tot_previ))+
  geom_point() +
  labs(x = "Age", y = "Dead/Total recorded")+
  xlim(0,30)+
  ylim(0,0.25)+ 
  theme_classic()

##############################
# Founders + updated younger #
##############################
all <- bind_rows(timeline_updated, founder_approximate) %>% 
  mutate(group = case_when(is.na(alive) == TRUE ~ "chick",
                           is.na(alive) == FALSE ~ "founder")) %>% 
  filter(death_year < 2100) %>% 
  filter(age > 1) %>%  # was previously 2
  group_by(age, group) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  arrange(age) %>% 
  mutate(dead_tot = n/180) %>% 
  mutate(csum = cumsum(dead_tot)) %>% 
  mutate(previ = lag(cumsum(n), default = 0)) %>% 
  mutate(dead_tot_previ = n/(180-previ))

ggplot(all, aes(x = age, y = csum, colour = group))+
  geom_point() +
  scale_colour_manual(values = c("black", "lightgrey"))+ 
  labs(x = "Age", y = "Cumulative proportion dead")+
  xlim(0,50)+
  ylim(0,1)+ 
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))

# Don't think this is actually informative without adjustment
# ggplot(all, aes(x = age, y = dead_tot, colour = group))+
#   geom_point() +
#   labs(x = "Age", y = "Prop")+
#   xlim(0,50)+
#   ylim(0,0.25)+ 
#   theme_classic()

ggplot(all, aes(x = age, y = dead_tot_previ, colour = group))+
  geom_point() +
  scale_colour_manual(values = c("black", "lightgrey"))+ 
  labs(x = "Age", y = "Died/Total that reached at least that age")+
  xlim(0,50)+
  ylim(0,0.1)+
  theme_classic()+
  theme(legend.position = "none",
                        axis.title = element_text(size = 18),
                        panel.background = element_rect(fill = NA, 
                                                        colour = "black", 
                                                        size = 0.7 ),
                        axis.text = element_text(size = 14))
  
# all %>% 
#   select(age, csum, dead_tot_previ) %>% 
#   rename(cumulative = csum, rate = dead_tot_previ) %>% 
#   write_delim("~/Documents/phd_python/deaths_data.csv", delim = ",")

# ########################################
# # Old approach to minimum founder ages #
# ########################################
# # Minimum age of founders that died.  First 7 are from characters.
# timeline_file %>% 
#   mutate(hatch_year = case_when(ID == "Lionel" ~ 1978,
#                                 ID == "Maggie" ~ 1978,
#                                 ID == "Barnard" ~ 1981,
#                                 ID == "Luke" ~ 1981,
#                                 ID == "Zephyr" ~ 1981,
#                                 ID == "Heather" ~ 1981,
#                                 ID == "Stumpy" ~ 1991,
#                                 ID == "Gunter" ~ 1981, # Age 6 or older
#                                 ID == "Pegasus" ~ 1984,
#                                 ID == "Ken" ~ 1988,
#                                 ID == "Bill" ~ 1982,
#                                 ID == "Richard_Henry" ~ 1975
#                                 
#   ))
# 

# 
# test <- characters %>% 
#   mutate(lay_year = str_extract(dateLaid, "(?<=/)[:digit:]{4}(?=\\s)")) %>% 
#   mutate(hatch_year = str_extract(dateHatched, "(?<=/)[:digit:]{4}(?=\\s)")) %>% 
#   mutate(death_year = str_extract(demise, "(?<=/)[:digit:]{4}(?=\\s)")) %>%
#   filter(hatch_year < 1995)




         