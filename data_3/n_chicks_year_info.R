library(tidyverse)

# bird_info <- read_csv('~/projects/kakapo-genomics/data/bird_info_rships.csv')
bird_infopre2019 <- read_csv("~/projects/kakapo-genomics/data/bird_info_rships_nov20.csv")
bird_info2019 <- read_delim("~/projects/kakapo-genomics/data/bird_info_rships_2019cohort.csv",
                            delim = "\t") %>% 
  mutate(hatch_year = 2019) %>% 
  filter(str_detect(ID, "19") == FALSE)

bird_info <- bind_rows(bird_infopre2019, bird_info2019)

female_n_chicks_season <- bird_info %>% 
  group_by(mother, hatch_year) %>% 
  summarise(n_chicks_year = n()) %>% 
  filter(is.na(mother) == FALSE)



male_n_chicks_season <- bird_info %>% 
  group_by(father, hatch_year) %>% 
  summarise(n_chicks_year = n()) %>% 
  filter(is.na(father) == FALSE)


mean(female_n_chicks_season$n_chicks_year)
sd(female_n_chicks_season$n_chicks_year)


mean(male_n_chicks_season$n_chicks_year)
sd(male_n_chicks_season$n_chicks_year)


# Easier to represent as % of breeders that have each number of chicks
table(female_n_chicks_season$n_chicks_year)
nrow(female_n_chicks_season)


table(male_n_chicks_season$n_chicks_year)
nrow(male_n_chicks_season)


multimating_males <- bird_info %>% 
  group_by(father, mother, hatch_year) %>% 
  summarise(n_chicks_in_pair = n()) %>% 
  ungroup() %>% 
  group_by(father, hatch_year) %>% 
  summarise(n_pairings = n()) %>% 
  filter(is.na(hatch_year) == FALSE)


table(multimating_males$n_pairings)
