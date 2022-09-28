library(tidyverse)



pop_size <- read_delim("~/Downloads/approx_kakapo_pop_size.csv", delim = ",") %>% 
  mutate(sim = "real")
  
ggplot(pop_size, aes(x = year, y = total_pop))+
  geom_line()

sims <- list.files(path = "/media/drive_6tb/projects/kakapo-genomics/output/40_growth_sims_off_mahuika/prop_f_test/prop_f_37o/", 
                   pattern = "pop_list_", full.names = TRUE) %>%
  setNames(., sub("\\.csv$", "", basename(.))) %>% 
  map(read_delim, delim = ",", trim_ws = TRUE)

summarise_sim <- function(num, data){
  name <- names(data[num])
  bind_rows(data[num]) %>% 
    group_by(year) %>%
    summarise(total_pop = n()) %>%
    mutate(year = year-24,
           sim = name)
}

sim_res <- lapply(1:length(sims), summarise_sim, data = sims) %>% 
  reduce(bind_rows) %>% 
  mutate(sim = str_extract(sim, "[:digit:]+"))


both <- bind_rows(pop_size, sim_res) %>% 
  filter(year < 2025) %>% 
  mutate(state = case_when(sim == "real" ~ "real",
                           TRUE ~ "sim"))

ggplot(both, aes(x = year, y = total_pop, group = sim, colour = state))+
  geom_line()+
  scale_colour_manual(values =  c("black","grey","#452736", "#dbcbbaff")) +
  labs(x = "Year",  y = "Total population") +
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)))

# cut off translocation effect
cutoff <- sim_res %>% 
  filter(total_pop >= 290) %>% 
  group_by(sim) %>% 
  summarise(cut_here = min(year))

cleaned <- left_join(sim_res, cutoff, by = "sim") %>% 
  filter(year < cut_here+1) %>% 
  bind_rows(pop_size) %>% 
  mutate(capped_pop = case_when(total_pop > 290 ~ 290,
                                TRUE ~ total_pop)) %>% 
  filter(year < 2025) %>% 
  mutate(state = case_when(sim == "real" ~ "real",
                           TRUE ~ "sim"))
  
ggplot(cleaned, aes(x = year, y = capped_pop, group = sim, colour = state))+
  geom_line()+
  scale_colour_manual(values =  c("black","grey","#452736", "#dbcbbaff")) +
  labs(x = "Year",  y = "Total population") +
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)))
  
## OLD

# real_pop_size <- read_delim("/media/drive_6tb/projects/kakapo-genomics/data/bird_info_timeline.csv", delim = ",") %>% 
#   mutate(hatch_year = case_when(is.na(hatch_year) == TRUE ~ 1980,
#                                 TRUE ~ hatch_year))
# 
# 
# count_alive_in_year <- function(year, dataset){
#   dataset %>% 
#     filter(hatch_year <= year & death_year > year) %>% 
#     nrow()
# }


# a <- real_pop_size %>% 
#   filter(hatch_year <= 1995 & death_year > 1995)
# 
# count_alive_in_year(1995, real_pop_size)
# 
# 
# years <-c(1980:2018)
# res <- lapply(years, count_alive_in_year, dataset = real_pop_size) %>% 
#   unlist()
# 
# pop <- tibble(year = years, pop_size = res)
# 
# ggplot(pop, aes(x = year, y = pop_size))+
#   geom_line()


# there seems to be an issue with the above


# # This isn't right either
# k_char <- read_delim("~/projects/kakapo-genomics/data/phenotypic_data/Characters/kakapo_characters.csv", delim = ",") %>% 
#   select(birdName, demise, dateHatched) %>% 
#   mutate(hatch = str_extract(dateHatched, "(?<=/)[:digit:]{4}(?=\\s)"),
#          death = str_extract(demise, "(?<=/)[:digit:]{4}(?=\\s)")) %>% 
#   select(birdName, hatch, death) %>% 
#   #filter((is.na(hatch) & is.na(death)) == FALSE) %>% 
#   filter(hatch < 1995 | is.na(hatch)) %>% 
#   filter(death > 1995 | is.na(death))


# # From Andrew Digby's plotly, need to figure out how to recreate though
# pop_trend_origin <- read_delim("~/Downloads/pop-trend-origin.csv", delim = ",", skip = 1, 
#                                col_names = c("year", "fiordland", "stewart", "transferred", "hatched")) %>% 
#    add_row(year = c(2018,2019,2020,2021), fiordland = c(0,0,0,0), stewart = c(0,0,0,0), transferred = c(0,0,0,0), hatched = c(0,0,210,201)) %>% 
#   mutate(total_pop = fiordland+stewart+transferred+hatched) %>% 
#   mutate(sim = "real")