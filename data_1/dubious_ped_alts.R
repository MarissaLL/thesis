#!/usr/bin/env Rscript

library(tidyverse)

#############
# FUNCTIONS #
#############

# Find the likely candidates for parents for the dubious trios
find_candidates <- function(indiv,  pedigree){
  chick_hatch_year <- pedigree %>% 
    filter(ID == indiv) %>% 
    select(hatch_year) %>% 
    pull()
  indiv_info <- pedigree %>% 
    filter(ID != indiv ) %>% 
    select(ID, hatch_year, sex) %>% 
    filter(hatch_year <= chick_hatch_year-2 | is.na(hatch_year) == TRUE) %>% 
    mutate(chick = indiv)
  return(indiv_info)
}

# Find the likely candidates for parents for the dubious trios where the chick hatch year is provided explicitly (not in the pedigree)
find_candidates_np <- function(indiv,  chick_hatch_year){
  indiv_info <- pedigree_file %>% 
    filter(ID != indiv ) %>% 
    select(ID, hatch_year, sex) %>% 
    filter(hatch_year <= chick_hatch_year-2 | is.na(hatch_year) == TRUE) %>% 
    mutate(chick = indiv) 
  return(indiv_info)
}

# Find all possible combinations of the candidate parents with each uncertain parentage chick
find_combinations <- function(indiv, indiv_candidates) {
  indiv_mothers <- indiv_candidates %>% 
    filter(chick == indiv) %>% 
    filter(sex == "F") %>% 
    select(ID)
  
  indiv_fathers <- indiv_candidates %>% 
    filter(chick == indiv) %>% 
    filter(sex == "M") %>% 
    select(ID)
  
  all_combos_nested <- crossing(mother = indiv_mothers, father = indiv_fathers) 
  all_combos <- all_combos_nested %>%
    mutate(chick = indiv,
           mother = pull(all_combos_nested$mother),
           father = pull(all_combos_nested$father))
  return(all_combos)
}


###########
# GLOBALS #
###########

files <- commandArgs(trailingOnly = TRUE)


pedigree_file <- read_csv(files[1]) 

options <- scan(files[2], what = "character", sep = "\t")

candidate_trios_output <- files[3]

mode <- files[4]

if (mode == "chick_in_pedigree") {
  chicks_to_check <- unlist(files[5:length(files)]) %>% 
    str_remove("\\[") %>% 
    str_remove("\\]") %>% 
    str_remove(",")
  message(paste0("Found these chicks to check: ", chicks_to_check))
  message("Expecting to find them in the pedigree")
  message(chicks_to_check[1])
}

if (mode == "chick_not_in_pedigree") {
  chicks_to_check <- files[5:(((length(files) - 4)/2)+4)]
  chick_hatch_years <- as.numeric(files[(((length(files) - 4)/2)+5):length(files)])
  message("Found these chicks to check: ", paste(chicks_to_check, collapse = ", "))
  message("Not expecting to find them in the pedigree")
  message("Hatch years are: ", paste(chick_hatch_years, collapse = ", "))
}



## Dev
# pedigree_file <- read_csv("data/bird_info_rships_nov20.csv")
# # pedigree_file <- read_csv("old_files/formatted_bird_list_manual.csv")
# # pedigree_file <- read_csv("data/konini.csv")
# # pedigree_file <- read_csv("data/bird_info_rships.csv")
# # #
# options <- scan("all_seq_birdnames.txt", what = "character", sep = "\t")
# candidate_trios_output <- "output/17_processing_np_chicks/np_chicks_full.txt"
# chicks_to_check <- c("Kuihi1B", "Kuihi3B", "Queenie3A")
# chicks_to_check <- c("Egilsay")
#
# candidate_trios_output <- "output/04_test_dubious_trios/candidate_trios_test_2020.txt"



########
# MAIN #
########


pedigree_data <-  pedigree_file %>% 
  mutate(ID = str_replace_all(ID, " ", "_"),
         mother = str_replace_all(mother, " ", "_"),
         father = str_replace_all(father, " ", "_")) %>% 
  filter(ID %in% options)
message("Checkpoint")
if (mode == "chick_in_pedigree"){
  candidates <- lapply(chicks_to_check, find_candidates, pedigree = pedigree_file) %>% 
    reduce(bind_rows)
}


if (mode == "chick_not_in_pedigree"){
  candidates <- mapply(find_candidates_np, chicks_to_check, chick_hatch_years, SIMPLIFY = FALSE) %>% 
    reduce(bind_rows)
}

message(candidates)
candidate_trios_test <- lapply(chicks_to_check, find_combinations, indiv_candidates = candidates) %>% 
  reduce(bind_rows)
message(candidate_trios_test)

candidate_trios_test <-  candidate_trios_test %>% 
  mutate(family_id = 1:nrow(candidate_trios_test))

# Unsure what the purpose of this was but it seems to break things so getting rid of for now
# candidate_trios_test <- candidate_trios_test[1:5,]

write_delim(candidate_trios_test, candidate_trios_output)


# Assign sex information based on true sex for offspring, and mother/father status for parents (not necessarily the real sex)
gathered_trios <- candidate_trios_test %>% 
  gather(type, individual, -family_id) 

pedigree_data <- pedigree_file %>% 
  mutate(ID = str_replace_all(ID, " ", "_"))

all_info <- inner_join(gathered_trios, pedigree_data, by = c(individual = "ID")) %>% 
  select(family_id, individual, sex, type) %>%
  mutate(sex2 = case_when(type == "mother" ~ 2,
                          type == "father" ~ 1,
                          type == "chick" & sex == "F" ~ 2,
                          type == "chick" & sex == "M" ~ 1)) %>% 
  select(-sex) %>% 
  rename(sex = sex2)
message("gets to here")



# Generate and output individuals to keep info files
keep_info <- all_info %>%
  group_by(family_id) %>%
  distinct(individual,.keep_all = TRUE) %>%
  select(individual) %>% 
  group_split()

# keep_info_no_FID <-  lapply(keep_info, FUN = function(x) subset(x, select=2)) 
names(keep_info) <- seq(1:length(keep_info))

lapply(names(keep_info),
       function(x){write_delim(keep_info[[x]],
                               path = paste0("family_", x ,"_keep_info.txt"),
                               col_names = FALSE)})
message("generates keep info")
# Generate and output id info files
id_info <- all_info %>%
  group_by(family_id) %>%
  distinct(individual,.keep_all = TRUE) %>%
  mutate(individual1 = individual, individual2 = individual) %>%
  mutate(plink2_FID = individual1) %>% 
  select(plink2_FID, individual1, family_id, individual2) %>%
  group_split()

# id_info <- lapply(id_info, FUN = function(x) subset(x, select=-1))


names(id_info) <- seq(1:length(id_info))


lapply(names(id_info),
       function(x){write_delim(id_info[[x]],
                               path = paste0("family_", x ,"_id_info.txt"),
                               col_names = FALSE)})

message("generates id info files")

# Generate and output parent info files
parent_info <- left_join(all_info, candidate_trios_test, by = c("individual"="chick" , "family_id" = "family_id")) %>%
  # rename(mother = ID, father = ID1) %>%
  select(family_id, individual, mother, father) %>%
  group_by(family_id, individual, mother, father) %>%
  distinct() %>%
  ungroup() %>%
  group_by(family_id) %>%
  group_split()

# parent_info <- lapply(parent_info, FUN = function(x) subset(x, select=-1))
names(parent_info) <- seq(1:length(parent_info))

lapply(names(parent_info),
       function(x){write_delim(parent_info[[x]],
                               path = paste0("family_", x ,"_parent_info.txt"),
                               col_names = FALSE)})


# Generate and output sex info files
sex_info <- all_info %>%
  group_by(family_id) %>%
  distinct(individual,.keep_all = TRUE) %>%
  select(family_id, individual, sex) %>%
  group_split()

# sex_info <- lapply(sex_info, FUN = function(x) subset(x, select=-1))

names(sex_info) <- seq(1:length(sex_info))


lapply(names(sex_info),
       function(x){write_delim(sex_info[[x]],
                               path = paste0("family_", x ,"_sex_info.txt"),
                               col_names = FALSE )})


message("generates sex info files. Script end")

################################
# Used to RUN m_err_dubious.sh # Now just let this script run through check_uncertain_parentage.nf
################################

#####################################
# PLOT with plot_dubious_ped_alts.R #
#####################################

