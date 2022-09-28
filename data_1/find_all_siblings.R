library(tidyverse)



ped <- read_delim("~/projects/kakapo-genomics/data/bird_info_rships.csv", delim = ",")

### Parent-offspring ###
offspring_mother <- ped %>% 
  select(ID, mother) %>% 
  rename(indiv1 = ID, indiv2 = mother) %>% 
  filter(is.na(indiv2) == FALSE) %>% 
  mutate(rship = "offspring_mother")
  
offspring_father <- ped %>% 
  select(ID, father) %>% 
  rename(indiv1 = ID, indiv2 = father) %>% 
  filter(is.na(indiv2) == FALSE) %>% 
  mutate(rship = "offspring_father")


father_offspring <- offspring_father %>% 
  rename(indiv1 = indiv2, indiv2 = indiv1)

mother_offspring <- offspring_mother %>% 
  rename(indiv1 = indiv2, indiv2 = indiv1)


parent_offspring <- bind_rows(offspring_mother, offspring_father, father_offspring, mother_offspring)


### Full siblings ###
full_sib_groups <- ped %>% 
  group_by(mother, father) %>% 
  summarise(num_sibs = n()) %>% 
  filter(is.na(mother) == FALSE, is.na(father) == FALSE, num_sibs > 1 ) %>% 
  rownames_to_column(var = "group_num") %>% 
  mutate(group_num = as.numeric(group_num)) %>% 
  right_join(ped) %>% 
  ungroup() %>% 
  filter(is.na(group_num) == FALSE)

get_sibs <- function(sib_group, grouped_sibs){
  sibs <- grouped_sibs %>% 
    filter(group_num == sib_group) %>% 
    select(ID) %>% 
    pull()
  res <- crossing(sibs, sibs, .name_repair = "unique") %>% 
    rename(indiv1 = `sibs...1`, indiv2 = `sibs...2`) %>% 
    filter(indiv1 != indiv2)
  return(res)
}


full_sibs <- lapply(1:max(full_sib_groups$group_num), get_sibs, grouped_sibs = full_sib_groups) %>% 
  reduce(bind_rows) %>% 
  mutate(rship = "full_sibs")


### Half sibs ###

all_names <- ped %>% select(ID) %>% pull()

get_half_sibs <- function(bird_name, df){
  name_row <- df %>% 
    filter(ID == bird_name)
  mother_name <- name_row$mother
  father_name <- name_row$father
  father_hs <- df %>% 
    filter(father == father_name & mother != mother_name) %>% 
    select(ID) %>% 
    pull()
  mother_hs <- df %>% 
    filter(mother == mother_name & father != father_name) %>% 
    select(ID) %>% 
    pull()
  all <- c(father_hs, mother_hs) %>% 
    enframe(name = NULL, value = "indiv2") %>% 
    mutate(indiv1 = bird_name, 
           rship = "half_sibs")
  return(all)
}


all_hs <- lapply(all_names, get_half_sibs, df = ped) %>% 
  reduce(bind_rows)


###### All relationships ####

all_rships <- bind_rows(parent_offspring, full_sibs, all_hs)


# Then also moved this file to data directory
write_delim(all_rships, "~/projects/kakapo-genomics/data/pairwise_relationships_old.txt", delim = "\t")



  