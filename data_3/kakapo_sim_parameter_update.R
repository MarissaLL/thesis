library(lubridate)
library(readxl)
library(tidyverse)
library(xtable)

#### Function ####

find_time_interval <- function(name, df, mode){
  if (mode == "bird"){filter_col = sym("bird")}
  if (mode == "island"){filter_col = sym("Island")} 
  intervals <- df %>% 
    select(!!filter_col, year) %>% 
    distinct() %>% 
    filter(!!filter_col == name) %>%
    select(year) %>%
    pull() %>% 
    diff() %>% 
    as.numeric()
  return(intervals)
}

##### Data #####

kc <- read_delim("~/Downloads/kakapo_characters.csv", delim = ",")


kf <- read_xlsx("~/projects/kakapo-genomics/data/phenotypic_data/Fertility/kakapo_fertility.xlsx")

#### Breeding interval on the three main islands ####

mini_kf <- kf %>% 
  select(bird, Island, year, Fledged)

island_breeding <- mini_kf %>% 
  group_by(Island, year) %>% 
  summarise(breeding = n()) %>% 
  ungroup()

print(xtable(island_breeding, digits = c(0,0,0,0)), include.rownames = F)

islands <- c("Anchor", "Hauturu", "Whenua Hou")

island_br_interval <- lapply(islands, 
                             find_time_interval,
                             df = island_breeding,
                             mode = "island") %>% 
  unlist()


table(island_br_interval)

# Remove an interval of 19 years from Hauturu-o-toi as kakapo may not have been on the island at the time.
isl_br_excl19 <- island_br_interval[island_br_interval != 19]
mean(isl_br_excl19)
median(isl_br_excl19)
sd(isl_br_excl19)

#### Interval between breeding for females

females_years_bred <- mini_kf %>% 
  group_by(bird, year) %>% 
  summarise(times_bred = n()) %>% 
  group_by(bird) %>% 
  summarise(years_bred = n())

# Filter for birds that have bred at least twice so the interval between their breeding can be taken
females_bred2plus <- females_years_bred %>% 
  filter(years_bred >= 2) %>% 
  select(bird) %>% 
  pull()



f2plus <- mini_kf %>% 
  filter(bird %in% females_bred2plus == TRUE) %>% 
  distinct()




#### Regardless of if breeding was successful
intervals <- lapply(females_bred2plus, 
                    find_time_interval, 
                    df = mini_kf,
                    mode = "bird") %>% 
  unlist()


table(intervals)
mean(intervals)
median(intervals)
sd(intervals)

hist(intervals)

##### Only counting successful breeding (fledged > 0)

fledged_only <- mini_kf %>% 
  filter(Fledged > 0) 

f2chickplus <- fledged_only %>% 
  group_by(bird, year) %>% 
  summarise(times_bred = n()) %>% 
  group_by(bird) %>% 
  summarise(years_bred = n()) %>% 
  filter(years_bred >= 2) %>% 
  select(bird) %>% 
  pull()


intervals_fldg <- lapply(f2chickplus,
                         find_time_interval,
                         df = fledged_only, 
                         mode = "bird") %>% 
  unlist()

table(intervals_fldg)


a <- enframe(table(intervals))
b <- enframe(table(intervals_fldg))


c <- full_join(a, b, by = c(name= "name")) %>% 
  rename(interval_length = name, all_breeding = value.x, fledged_only = value.y) %>% 
  pivot_longer(cols = c(all_breeding, fledged_only), names_to = "which_breeding", values_to = "n") %>% 
  mutate(interval_length =as.numeric(interval_length))



ggplot(c, aes(x = interval_length, y = n ))+
  facet_wrap(~ which_breeding,nrow = 2)+ 
  geom_histogram(stat = "identity") + 
  labs(x = "Interval between breeding in female kakapo", y = "Occurrences") +
  theme_classic()+
  theme(legend.position = "none",
      axis.title = element_text(size = 14),
      panel.background = element_rect(fill = NA, 
                                      colour = "black", 
                                      size = 0.7 ),
      axis.text = element_text(size = 14))


mean(intervals)
mean(intervals_fldg)

sd(intervals)
sd(intervals_fldg)

median(intervals)
median(intervals_fldg)

#### Age at death

# Should really distinguish laid from discovered
# Strip hm data before importing date because it is inconsistent 
mort_data <- kc %>% 
  select(birdName, alive, dateLaid, discoveryDate, demise, sex) %>% 
  filter(alive == FALSE) %>% 
  # mutate(first_info = dateLaid) %>% # Trying with just known hatchdates
  mutate(first_info = case_when(is.na(dateLaid) == TRUE ~ discoveryDate,
                                TRUE ~ dateLaid)) %>%
  mutate(info_from =  case_when(is.na(dateLaid) == TRUE ~ "discovered",
                                TRUE ~ "laid")) %>%
  filter(is.na(first_info) == FALSE & is.na(demise) == FALSE) %>% 
  select(birdName, first_info, demise, info_from, sex) %>%
  mutate(first_info = dmy(str_remove(first_info, 
                                  "[:digit:]{1,2}:[:digit:]{2}$"))) %>% 
  mutate(demise = dmy(str_remove(demise, 
                             "[:digit:]{1,2}:[:digit:]{2}$"))) %>% 
  mutate(lifespan_days = interval(first_info, demise) %/% days(1)) %>% 
  mutate(lifespan_years = interval(first_info, demise) %/% years(1))


ggplot(data = subset(mort_data, 
                     lifespan_years > 0), 
       aes(x = lifespan_years, fill = info_from)) +
  geom_histogram(binwidth = 1) +
  theme_classic()+
  theme(legend.position = "top",
        axis.title = element_text(size = 14),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))

# Non-founders
table(mort_data$lifespan_years[mort_data$info_from == "laid"])

mort_adult <- mort_data %>% 
  filter(info_from == "laid") %>% 
  filter(lifespan_years > 0)

mean(mort_adult$lifespan_years)
sd(mort_adult$lifespan_years)


# Founders
table(mort_data$lifespan_years[mort_data$info_from == "discovered"])

mort_adult <- mort_data %>% 
  filter(info_from == "discovered") %>% 
  filter(lifespan_years > 0)

mean(mort_adult$lifespan_years)
sd(mort_adult$lifespan_years)

#Combined
table(mort_data$lifespan_years)

mort_adult <- mort_data %>% 
  filter(lifespan_years > 0)

mean(mort_adult$lifespan_years)
sd(mort_adult$lifespan_years)



## Lifetime offspring 
full_life <- mort_data %>% 
  mutate(demise_year = str_extract(demise, "^[:digit:]{4}"),
         hatch_year = str_extract(first_info, "^[:digit:]{4}")) %>% 
  filter(demise_year > 1995) %>% 
  filter(hatch_year != demise_year) 


full_life_m <- full_life %>%
  filter(sex == "Male") %>% 
  pull(birdName)

full_life_f <- full_life %>%
  filter(sex == "Female") %>% 
  pull(birdName)

bird_info <- read_delim("data/bird_info_rships_nov20.csv", delim = ",")

fn <- function(name, bird_info){
  aa <- bird_info %>% 
    filter(mother == name | father == name )
  tibble(name = name,
         offspring = nrow(aa))
  
}

lifetime_offspring_m <- lapply(full_life_m, fn, bird_info) %>% 
  reduce(bind_rows)

mean(lifetime_offspring_m$offspring)
sd(lifetime_offspring_m$offspring)

lifetime_offspring_f <- lapply(full_life_f, fn, bird_info) %>% 
  reduce(bind_rows)

mean(lifetime_offspring_f$offspring)
sd(lifetime_offspring_f$offspring)




bird_info %>% 
  filter(mother %in% full_life | father %in% full_life)



