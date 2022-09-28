#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)
options(scipen = 999) # Prevent numbers converting to scientific notation - causes matching problems for locus IDs

###########
# GLOBALS #
###########

# files <- commandArgs(trailingOnly = TRUE)
# 
# all_trios_file <- files[1]
# trios_run_file <- files[2]
# me_vcf_path <- "."
# family_num <- files[3]



# all_trios_file <- "/cifs/biocldap/student_users/marissalelec/projects/kakapo-test/output/07_offspring_m_err/all_true_trio_parents.txt"
all_trios_file <- "/media/drive_2tb/projects/kakapo-genomics/07_offspring_m_err/all_true_trio_parents.txt"
trios_run_file <- "/media/drive_2tb/projects/kakapo-genomics/07_offspring_m_err/all_run_trios_part_num.txt"
me_vcf_path <- "/media/drive_2tb/projects/kakapo-genomics/07_offspring_m_err/"


family_num <-  16

########
# MAIN #
########

#################################
# Section 1: Family information #
#################################

# Identify and assign an ID to multisibling families from the full dataset of trios
# Actually only going to run for families 1 to 16, which are the families with at least 3 full siblings
multisib_fams <- read_delim(all_trios_file, delim = "\t", trim_ws = TRUE, col_names = FALSE) %>%
  filter(is.na(X3) == FALSE) %>% 
  mutate(parents = str_c(X3, X4, sep = ",")) %>% 
  group_by(parents) %>% 
  count() %>% 
  filter(n > 1) %>% 
  mutate(father = str_extract(parents, "[:graph:]+(?=,)"),
         mother = str_extract(parents, "(?<=,)[:graph:]+")) %>% 
  ungroup() %>% 
  select(father, mother, n) %>% 
  arrange(desc(as.numeric(n))) %>% 
  mutate(multisib_fam_id = paste0("family_", seq(1:length(father))))


# Read in the file with the mendel run number of each trio (run through 07_offspring_me.nf)
trios_run <- read_delim(trios_run_file, delim = '\t')


# Put together all the info about which IDs all the trios in a multisib family were run under
multisib_trio_runs <- left_join(trios_run, multisib_fams, by = c("mother", "father"))  


# CHECK. There should be no multisib families in here
not_multisib <- multisib_trio_runs %>% 
  filter(is.na(multisib_fam_id) == TRUE)

# Extract details for one family (at a time). 
multisib_single_family <- multisib_trio_runs %>% 
  filter(is.na(multisib_fam_id) == FALSE) %>% 
  filter(multisib_fam_id == paste0("family_", family_num))

relevant_file_nums <- multisib_single_family$random_group 

father <- multisib_single_family$father %>% 
  unique()

mother <- multisib_single_family$mother %>% 
  unique()

all_of_the_family <- c(mother, father, multisib_single_family$ID)

# End of section, remove excess stuff from environment
rm(trios_run, multisib_trio_runs, not_multisib)


##########################################
# Section 2: Mendelian error information #
##########################################

# Read in mendelian error info for the one family
all_me_in_multisib_family <- lapply(multisib_single_family$random_group, function(group_id) {
  read_delim(file = paste0(me_vcf_path, "all_details_", group_id, ".mendel"), 
             delim = ' ', trim_ws = TRUE, skip = 1, col_names = FALSE) %>% 
    select(-X1, -X3, -X7, -X9) %>% 
    rename_at(vars(-contains("X4")),
              funs(str_c(., "_file_", group_id))) 
  }) 


# Convert list into tibble and give columns more helpful names
all_me_flat <- all_me_in_multisib_family %>% 
  reduce(full_join, by = "X4") %>% 
  rename_at(vars(contains("X2")),
            funs(str_replace(., "X2", "chick_name"))) %>% 
  rename(locus = X4) %>% 
  rename_at(vars(contains("X5")),
            funs(str_replace(., "X5", "m_err_code"))) %>% 
  rename_at(vars(contains("X6")),
            funs(str_replace(., "X6", "parent1"))) %>% 
  rename_at(vars(contains("X8")),
            funs(str_replace(., "X8", "parent2"))) %>% 
  rename_at(vars(contains("X10")),
            funs(str_replace(., "X10", "chick"))) 
  

# Extract a vector of all the M.E. sites 
me_sites <- all_me_flat %>% 
  select(locus) %>% 
  unique() %>% 
  pull()

# Make a file of loci of interest (with M.E.) for grep to use (for more efficient read in of other data)
me_locus_file <- all_me_flat %>%
  select(locus) %>% 
  mutate(CHROM = str_extract(locus, "[:graph:]+(?=:[:digit:]+)")) %>% 
  mutate(POS = str_extract(locus, "(?<=:)[:digit:]+")) %>% 
  select(-locus)

write_delim(me_locus_file, paste0("me_locus_file_", family_num ,".txt"), delim = '\t')


# Provided the above checks look fine, next step is to tidy the M.E. output to reduce the number of duplicate columns
# Possibly consider just getting the relevant genotypes for every individual here so that they aren't missing in any
# Alternately, may not need genotypes at all, just the error code?

all_m_err_cols <- paste0("m_err_code_file_", sample(relevant_file_nums, length(relevant_file_nums), replace = FALSE))


me_attribution <- all_me_flat %>% 
  select(locus, all_of(all_m_err_cols)) %>% 
  mutate_at(vars(contains("m_err_code_file_")),
            funs(str_replace_all(., "2", "unattributable"))) %>% 
  mutate_at(vars(contains("m_err_code_file_")),
            funs(str_replace_all(., "5", "unattributable"))) %>%
    mutate_at(vars(contains("m_err_code_file_")),
            funs(str_replace_all(., "6", "father_chick"))) %>% 
  mutate_at(vars(contains("m_err_code_file_")),
            funs(str_replace_all(., "7", "mother_chick"))) %>% 
  replace(is.na(.), "no_error")

# Summarise the number of each attribution per site
me_attribution_wide <- me_attribution %>% 
  gather(var, value, -locus) %>% 
  spread(locus, value) 


unattributable_me_list <- lapply(me_sites, function(locus_name) {
  sum (str_count(as.matrix(me_attribution_wide[locus_name]), fixed("unattributable")))}) 

father_chick_me_list <- lapply(me_sites, function(locus_name) {
  sum (str_count(as.matrix(me_attribution_wide[locus_name]), fixed("father_chick")))}) 

mother_chick_me_list <- lapply(me_sites, function(locus_name) {
  sum (str_count(as.matrix(me_attribution_wide[locus_name]), fixed("mother_chick")))}) 

no_error_locus_list <- lapply(me_sites, function(locus_name) {
  sum (str_count(as.matrix(me_attribution_wide[locus_name]), fixed("no_error")))}) 

unattributable_me <- unlist(unattributable_me_list)
father_chick_me <-  unlist(father_chick_me_list)
mother_chick_me <-  unlist(mother_chick_me_list)
no_error_locus <-  unlist(no_error_locus_list)

summarised_attributions <- list(me_sites,unattributable_me,father_chick_me, mother_chick_me, no_error_locus) %>% 
  bind_cols() %>% 
  rename(locus = `...1`, unattributable = `...2`, father_chick = `...3`, mother_chick = `...4`, no_error = `...5`)


# Join back with the dataset of attribution codes
joined_me_codes_summary <- full_join(me_attribution, summarised_attributions, by = c(locus = "locus"))

# Rename cols with info about individual
me_attribution_named <- joined_me_codes_summary #%>% 
  # rename_at(vars(contains("file")), 
  #           funs(str_replace(., "_file_[:digit:]+", "" )))


# End of section, remove excess stuff from environment
rm(me_attribution_wide, unattributable_me_list,father_chick_me_list, mother_chick_me_list, no_error_locus_list,
   unattributable_me, father_chick_me, mother_chick_me, no_error_locus, me_attribution, summarised_attributions)



##############################
# Section 3/4: Variants & GQ #
##############################


# This pre-processing line no longer works????
# Import vcf lines for the mendelian error variants
multisib_family_variants_gq <- lapply(multisib_single_family$random_group, function(group_id) {
  vcf_to_import <- paste0(me_vcf_path,"part_", group_id, "_keep_indivs.vcf")
  
  preprocessing_command <- sprintf("grep --file %s %s", 
                                   paste0("me_locus_file_", family_num, ".txt"), 
                                   paste(vcf_to_import, collapse = " "))
  
  vcf <- data.table::fread(cmd = preprocessing_command) %>%
  select(-ID, -QUAL, -FILTER, -INFO, -FORMAT) %>%
    rename(CHROM = `#CHROM`) %>%
    rename_at(vars(-matches("CHROM|POS|REF|ALT")),
              funs(str_c(., "_indiv", sep = ""))) %>%
  mutate(locus = str_c(CHROM, POS, sep = ":")) %>%
  select(-CHROM, -POS) %>%
  rownames_to_column()
  return(vcf)})

# Flatten, remove duplicated columns (the parents have been read in multiple times), remove duplicated rows
# Why does filter locus %in% me sites remove any loci????????????????????????
flat_vcf <- multisib_family_variants_gq %>%
  reduce(full_join, by = c("locus", "rowname"), suffix = c(".x", ".y")) %>% 
  select(-ends_with(".x"), -ends_with(".y")) #%>% 
  # rename_at(vars(ends_with(".x")),
  #           funs(str_remove(.,"\\.x"))) 

# This is necessary in the cases where reduce decides to not add ".x" or ".y" to some duplicated columns
# flat_vcf <- flat_vcf %>% 
#   subset(select=which(!duplicated(names(.))))




separate_var_gq <- function(tibble, indiv_name) {
  separate(data = tibble, 
           col = paste0(indiv_name, "_indiv"), 
           into = c(paste0(indiv_name, "_var"), paste0(indiv_name, "_gq")), 
           sep = ":",  #"(?<=[:graph:]{3}):(?=[:graph:]+)"
           extra = "drop") %>% 
    select_at(vars(matches("_gq|_var|locus|rowname")))
}

flat_var_gq <- lapply(all_of_the_family, function(individual){
  separate_var_gq(flat_vcf, individual)
}) %>% 
  reduce(full_join, by = c("locus", "rowname"))

separate_alleles <- function(tibble, indiv_name) {  
  separate(data = tibble, col = paste0(indiv_name, "_var"), 
           into = c(paste0(indiv_name, "_allele1"), paste0(indiv_name, "_allele2")), 
           sep = "[/\\|]") %>% 
    select_at(vars(matches("_allele[12]|_gq|locus|rowname")))
}

separated_alleles <- lapply(all_of_the_family, function(individual){
  separate_alleles(flat_var_gq, individual)
}) %>% 
  reduce(full_join, by = c("locus", "rowname"), suffix = c(".x", ".y"))  %>% 
  select(-ends_with(".x"), -ends_with(".y"))  #%>% 
  # rename_at(vars(ends_with(".x")),
            # funs(str_remove(.,"\\.x")))

# This is necessary in the cases where reduce decides to not add ".x" or ".y" to some duplicated columns
separated_alleles <- separated_alleles %>% 
  subset(select=which(!duplicated(names(.))))


sep_alleles_info <- full_join(flat_vcf, separated_alleles, by = c("locus", "rowname")) %>%
  select_at(vars(-ends_with("_indiv")))



# Find the highest alt allele. Could use this to determine how many columns to make.
# Currently just using this to 
### FIND A BETTER WAY ####
### Apparently this doesn't even work ####
allele_cols_names <-  sep_alleles_info %>% 
  select_at(vars(matches("_allele[12]"))) %>% 
  names()

locus_max <-   sep_alleles_info %>%  mutate(most = pmax(!!!rlang::syms(allele_cols_names)))
most_alts <- as.numeric(max(locus_max$most))


# Select one line for each locus based on having the least missingness
added_counting_vcf <- sep_alleles_info %>%
  unite(for_counting, matches("[^locus|^rowname|^REF|^ALT]" ), sep = "") %>% 
  mutate(num_dots = str_count(for_counting, "\\.")) 

#################### NEED TO THINK MORE ABOUT WHY THERE ARE EXTRA LOCI ^^ HERE THAT SHOULDN'T BE ######  
# Is grep treating the two values I give it as an or?

numeric_allele_rep <- full_join(sep_alleles_info, added_counting_vcf, by = c("locus", "rowname", "REF", "ALT")) %>% 
  filter(locus %in% me_sites == TRUE) %>% 
  group_by(locus) %>% 
  top_n(n = -1, wt = num_dots) %>% 
  select(-for_counting, -rowname, -num_dots) 



## Fill the column with the actual alleles that individual has.
# Some of the unused alts will be dropped. Need to double check this is working properly
alt_allele_cols <- numeric_allele_rep %>% 
  separate(col = ALT, 
           into = paste0("alt_", 1:5),
           sep = ",",
           fill = "right",
           extra = "drop") %>% 
  rownames_to_column() %>% 
  group_by(locus) 

gqs_to_rejoin <- alt_allele_cols %>% 
  select_at(vars(matches("_gq|rowname")))



# Function to convert numeric allele representation into alphabetical representation
alpha_alleles <- function(df, name) {
  mutate(df, !!name := case_when(!!sym(name) == 0 ~ REF,
                                 !!sym(name) == 1 ~ alt_1,
                                 !!sym(name) == 2 ~ alt_2,
                                 !!sym(name) == 3 ~ alt_3,
                                 !!sym(name) == 4 ~ alt_4,
                                 !!sym(name) == 5 ~ alt_5,
                                 !!sym(name) == "." ~ ".",
                                 TRUE ~ "dropped"))}


letter_alleles <- lapply(allele_cols_names, function (col_name) {
  alt_allele_cols %>% 
  select(locus, rowname, REF, alt_1, alt_2, alt_3, alt_4, alt_5, col_name) %>% 
  alpha_alleles(col_name)}) %>% 
  reduce(full_join, by = c("locus", "rowname", "REF", "alt_1","alt_2", "alt_3", "alt_4", "alt_5")) %>% 
  select(-REF, -alt_1, -alt_2, -alt_3, -alt_4, -alt_5)

var_gq_all <- full_join(letter_alleles, gqs_to_rejoin, by = c("locus", "rowname")) %>% 
  select(-rowname)


# Check if any alternate alleles have been dropped
if(any(var_gq_all=="dropped")){message("Alleles have been dropped. Fix this")}

# End of section, remove excess stuff from environment
rm(multisib_family_variants_gq, flat_vcf, separated_alleles, sep_alleles_info, locus_max, added_counting_vcf, numeric_allele_rep,  alt_allele_cols, gqs_to_rejoin)


#########################
# Section 5: Merge info #
#########################


## Get rid of the problem loci for now
## WHY ARE THESE NOT THE SAME??????????

# Loci which were duplicated in the vcf and can't be resolved easily
# Get rid of them for now
problem_loci_names <- var_gq_all %>% 
  group_by(locus) %>% 
  count() %>% 
  filter(n >1) %>% 
  select(locus) %>% 
  pull()

me_to_merge <- me_attribution_named %>%
  filter(locus %in% problem_loci_names == FALSE) 

var_gq_to_merge <- var_gq_all %>% 
  filter(locus %in% problem_loci_names == FALSE)

# Make a giant tbl
giant <- full_join(me_to_merge, var_gq_to_merge, by = "locus")


write_delim(x = giant,
            path = paste0(me_vcf_path, "me_var_gq_family_", family_num, ".txt"),
            delim = "\t")



#############################
# Section 6: Locate mistake #
#############################
giant <- read_delim( paste0(me_vcf_path, "me_var_gq_family_", family_num, ".txt"),
                    delim = "\t") %>% 
  mutate_at(vars(ends_with("_gq")),
            funs(as.numeric(.)))

# # Take a subset of the giant to work with
# smaller <- giant[1:1000,] %>% 
#   mutate_at(vars(ends_with("_gq")),
#             funs(as.numeric(.)))

# Rename the me_attribution code columns
for_renaming <- multisib_single_family %>%
  select(ID, random_group)

lookup_rename <- function(df, column_lookup) {
  df2 <- df
  names(df2)[names(df2) %in% paste0("m_err_code_file_", column_lookup$random_group)] = paste0("m_err_code_file_", column_lookup$ID[match(names(df2)
                                                                                                                                         [names(df2) %in% paste0("m_err_code_file_", column_lookup$random_group)], paste0("m_err_code_file_", column_lookup$random_group))])
  df2
}

renamed_cols <-  lookup_rename(giant, for_renaming)

add_parents <- renamed_cols %>% 
  rename_at(vars(matches(father)),
            funs(str_c(., "_father", sep = ""))) %>% 
  rename_at(vars(matches(mother)),
            funs(str_c(., "_mother", sep = "")))




################# All chicks show error #####################
father_mistake_all <- add_parents %>% 
  filter(father_chick > 0 & 
           mother_chick == 0 & 
           unattributable == 0 &
           no_error == 0) %>% 
  mutate(individual_w_error = father) %>% 
  mutate(across(ends_with("_gq"), as.numeric)) %>% 
  mutate(chick_gq_sum = select(., ends_with("_gq")) 
          %>% 
           rowSums())

mother_mistake_all <- add_parents %>% 
  filter(father_chick == 0 & 
           mother_chick > 0 & 
           unattributable == 0 &
           no_error == 0) %>% 
  mutate(individual_w_error = mother) %>% 
  mutate(across(ends_with("_gq"), as.numeric)) %>% 
  mutate(chick_gq_sum = select(., ends_with("_gq")) 
         %>% rowSums())



add_conf_col <- function(df, which_parent, parent_name) {
  df %>% 
    mutate(trust = case_when(as.numeric(!!sym(paste0(parent_name, "_gq_", which_parent))) < chick_gq_sum ~ "high",
                             as.numeric(!!sym(paste0(parent_name, "_gq_", which_parent))) >= chick_gq_sum ~ "low"))
  
  
}

father_mistake_all <- add_conf_col(father_mistake_all, "father", father)
mother_mistake_all <- add_conf_col(mother_mistake_all, "mother", mother)

###################### Some chicks show error, some show no error ##########################

father_mistake_and_no_error <- add_parents %>% 
  filter(father_chick > 0 & 
           mother_chick == 0 & 
           unattributable == 0 &
           no_error > 0)


mother_mistake_and_no_error <- add_parents %>% 
  filter(father_chick == 0 & 
           mother_chick > 0 & 
           unattributable == 0 &
           no_error > 0)

make_subset_interesting_indivs <- function (row, df, parent) {
  chicks_w_no_error <-  df %>% 
    rownames_to_column() %>% 
    filter(rowname == row) %>% 
    select_if(. %in% 'no_error') %>% 
    rename_all(., funs(str_remove(., "m_err_code_file_"))) %>% 
    names()
  mother_col_suffices <- c("_allele1_mother", "_allele2_mother", "_gq_mother")
  father_col_suffices <- c("_allele1_father", "_allele2_father", "_gq_father")
  all_related_cols_excl_par <- c("locus",
                                 paste0(rep_len(chicks_w_no_error,length.out = length(chicks_w_no_error)),"_allele1"),
                                 paste0(rep_len(chicks_w_no_error,length.out = length(chicks_w_no_error)),"_allele2"),
                                 paste0(rep_len(chicks_w_no_error,length.out = length(chicks_w_no_error)),"_gq"))
  
  if (parent == "mother") {
    all_related_cols <-  c(all_related_cols_excl_par, paste0(mother, rep_len(mother_col_suffices, length.out = length(mother_col_suffices))))}
  if (parent == "father"){
    all_related_cols <- c(all_related_cols_excl_par, paste0(father, rep_len(father_col_suffices, length.out = length(father_col_suffices))))}
  
  rel_cols_df <- df %>% 
    rownames_to_column() %>%
    filter(rowname == row) %>% 
    select(all_related_cols) 
  return(rel_cols_df)
}

mmne <- lapply(1:nrow(mother_mistake_and_no_error), make_subset_interesting_indivs, df = mother_mistake_and_no_error, parent = "mother")
fmne <- lapply(1:nrow(father_mistake_and_no_error), make_subset_interesting_indivs, df = father_mistake_and_no_error, parent = "father")




assign_allele <- function(individual, parent_name, parent, list_df, row_or_sublist, all_chicks) {
  
  single_row <- list_df[[row_or_sublist]] %>% 
    mutate(!!sym(paste0(individual[1],"_thoughts")) := case_when(!!sym(paste0(individual[1],"_allele1")) == "." & !!sym(paste0(individual[1],"_allele2")) == "." ~ "no_info",
                                                                 str_detect(!!sym(paste0(individual[1],"_allele1")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
                                                                 str_detect(!!sym(paste0(individual[1],"_allele2")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
                                                                 str_detect(!!sym(paste0(individual[1],"_allele1")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
                                                                 str_detect(!!sym(paste0(individual[1],"_allele2")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
                                                                 TRUE ~ "do_not_match"))
  if (length(individual) >= 2){
    single_row <- single_row %>%
      mutate(!!sym(paste0(individual[2],"_thoughts")) := case_when(!!sym(paste0(individual[2],"_allele1")) == "." & !!sym(paste0(individual[2],"_allele2")) == "." ~ "no_info",
                                                                   str_detect(!!sym(paste0(individual[2],"_allele1")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[2],"_allele2")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[2],"_allele1")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[2],"_allele2")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
                                                                   TRUE ~ "do_not_match"))}
  
  if (length(individual) >= 3){
    single_row <- single_row %>% 
      mutate(!!sym(paste0(individual[3],"_thoughts")) := case_when(!!sym(paste0(individual[3],"_allele1")) == "." & !!sym(paste0(individual[3],"_allele2")) == "." ~ "no_info",
                                                                   str_detect(!!sym(paste0(individual[3],"_allele1")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[3],"_allele2")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[3],"_allele1")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[3],"_allele2")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
                                                                   TRUE ~ "do_not_match"))}
  if (length(individual) >= 4){
    single_row <- single_row %>% 
      mutate(!!sym(paste0(individual[3],"_thoughts")) := case_when(!!sym(paste0(individual[4],"_allele1")) == "." & !!sym(paste0(individual[4],"_allele2")) == "." ~ "no_info",
                                                                   str_detect(!!sym(paste0(individual[4],"_allele1")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[4],"_allele2")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[4],"_allele1")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[4],"_allele2")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
                                                                   TRUE ~ "do_not_match"))}
  if (length(individual) >= 5){
    single_row <- single_row %>% 
      mutate(!!sym(paste0(individual[3],"_thoughts")) := case_when(!!sym(paste0(individual[5],"_allele1")) == "." & !!sym(paste0(individual[3],"_allele2")) == "." ~ "no_info",
                                                                   str_detect(!!sym(paste0(individual[5],"_allele1")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[5],"_allele2")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[5],"_allele1")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[5],"_allele2")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
                                                                   TRUE ~ "do_not_match"))}
  if (length(individual) >= 6){
    single_row <- single_row %>% 
      mutate(!!sym(paste0(individual[3],"_thoughts")) := case_when(!!sym(paste0(individual[6],"_allele1")) == "." & !!sym(paste0(individual[3],"_allele2")) == "." ~ "no_info",
                                                                   str_detect(!!sym(paste0(individual[6],"_allele1")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[6],"_allele2")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[6],"_allele1")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
                                                                   str_detect(!!sym(paste0(individual[6],"_allele2")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
                                                                   TRUE ~ "do_not_match"))}
  
  if (length(individual) > 6){ message("Need to add possibility of there being more chicks to the function")}
  missing_chicks <-  all_chicks[!all_chicks %in% individual]
  missing_colname <- c(paste0(rep(missing_chicks), "_thoughts"))
  
  single_row <- single_row %>% 
    add_column(!!!set_names(as.list(rep("involved_in_error", length(missing_colname))), nm = missing_colname))
  
  sum_cbc <- single_row[1,] %>%  str_count("could_be_correct") %>% sum() %>% as.double()
  
  sum_ni <- single_row[1,] %>% str_count("no_info") %>% sum() %>% as.double()
  
  sum_iie <- single_row[1,] %>% str_count("involved_in_error") %>% sum() %>% as.double()
  
  
  
  single_row <- single_row %>% 
    mutate(gt_possible_counts = sum_cbc,
           gt_no_info = sum_ni,
           gt_error = sum_iie)  
  return(single_row)
}



all_chicks <- all_of_the_family[!all_of_the_family %in% c(mother, father)]




assess_no_error_chicks <- function (row_or_sublist, list_df, parent_name, parent, all_chicks) {  
  got_names <- list_df[[row_or_sublist]] %>% 
    select_at(vars(ends_with("_allele1"))) %>% 
    rename_at(vars(ends_with("_allele1")),
              funs(str_remove(.,"_allele1"))) %>% 
    names()
  all_chicks <- all_chicks
  trial <- assign_allele(got_names, parent_name, parent,  list_df, row_or_sublist, all_chicks)
  return(trial)
}

mmne_thoughts <- lapply(1:length(mmne), assess_no_error_chicks, list_df = mmne, parent_name = mother, parent = "mother", all_chicks = all_chicks) %>% 
  bind_rows() %>% 
  mutate_at(vars(ends_with("_thoughts")), funs(str_replace_na(., replacement = "involved_in_error"))) %>% 
  select_at(vars(matches("locus|_thoughts|gt_"))) %>% 
  mutate(individual_w_error = case_when(gt_error > gt_possible_counts ~ mother,
                                        gt_error == 1 & gt_possible_counts > 1 ~ "other",
                                        TRUE ~ "hard_to_attribute")) %>% 
  mutate(trust = case_when(individual_w_error == mother & gt_no_info < gt_error ~ "high",
                           individual_w_error == mother & gt_no_info > gt_error ~ "low",
                           individual_w_error == "other" & gt_no_info < gt_possible_counts ~ "high",
                           individual_w_error == "other" & gt_no_info > gt_possible_counts ~ "low",
                           TRUE ~ "low"))

fmne_thoughts <- lapply(1:length(fmne), assess_no_error_chicks, list_df = fmne, parent_name = father, parent = "father", all_chicks = all_chicks) %>% 
  bind_rows() %>% 
  mutate_at(vars(ends_with("_thoughts")), funs(str_replace_na(., replacement = "involved_in_error"))) %>% 
  select_at(vars(matches("locus|_thoughts|gt_"))) %>% 
  mutate(individual_w_error = case_when(gt_error > gt_possible_counts ~ father,
                                        gt_error == 1 & gt_possible_counts > 1 ~ "other",
                                        TRUE ~ "hard_to_attribute")) %>% 
  mutate(trust = case_when(individual_w_error == father & gt_no_info < gt_error ~ "high",
                           individual_w_error == father & gt_no_info > gt_error ~ "low",
                           individual_w_error == "other" & gt_no_info < gt_possible_counts ~ "high",
                           individual_w_error == "other" & gt_no_info > gt_possible_counts ~ "low",
                           TRUE ~ "low"))


if(length(all_chicks) == 3){
  mmne_attr <- mmne_thoughts %>% 
    mutate(individual_w_error = case_when(individual_w_error == "other" & !!sym(paste0(all_chicks[1], "_thoughts")) == "involved_in_error" ~ all_chicks[1],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[2], "_thoughts")) == "involved_in_error" ~ all_chicks[2],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[3], "_thoughts")) == "involved_in_error" ~ all_chicks[3],
                                          TRUE ~ individual_w_error))
  fmne_attr <- fmne_thoughts %>% 
    mutate(individual_w_error = case_when(individual_w_error == "other" & !!sym(paste0(all_chicks[1], "_thoughts")) == "involved_in_error" ~ all_chicks[1],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[2], "_thoughts")) == "involved_in_error" ~ all_chicks[2],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[3], "_thoughts")) == "involved_in_error" ~ all_chicks[3],
                                          TRUE ~ individual_w_error))}

if(length(all_chicks) == 4){
  mmne_attr <- mmne_thoughts %>% 
    mutate(individual_w_error = case_when(individual_w_error == "other" & !!sym(paste0(all_chicks[1], "_thoughts")) == "involved_in_error" ~ all_chicks[1],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[2], "_thoughts")) == "involved_in_error" ~ all_chicks[2],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[3], "_thoughts")) == "involved_in_error" ~ all_chicks[3],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[4], "_thoughts")) == "involved_in_error" ~ all_chicks[4],
                                          TRUE ~ individual_w_error))
  fmne_attr <- fmne_thoughts %>% 
    mutate(individual_w_error = case_when(individual_w_error == "other" & !!sym(paste0(all_chicks[1], "_thoughts")) == "involved_in_error" ~ all_chicks[1],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[2], "_thoughts")) == "involved_in_error" ~ all_chicks[2],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[3], "_thoughts")) == "involved_in_error" ~ all_chicks[3],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[4], "_thoughts")) == "involved_in_error" ~ all_chicks[4],
                                          TRUE ~ individual_w_error))}


if(length(all_chicks) == 5){
  mmne_attr <- mmne_thoughts %>% 
    mutate(individual_w_error = case_when(individual_w_error == "other" & !!sym(paste0(all_chicks[1], "_thoughts")) == "involved_in_error" ~ all_chicks[1],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[2], "_thoughts")) == "involved_in_error" ~ all_chicks[2],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[3], "_thoughts")) == "involved_in_error" ~ all_chicks[3],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[4], "_thoughts")) == "involved_in_error" ~ all_chicks[4],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[5], "_thoughts")) == "involved_in_error" ~ all_chicks[5],
                                          TRUE ~ individual_w_error))
  fmne_attr <- fmne_thoughts %>% 
    mutate(individual_w_error = case_when(individual_w_error == "other" & !!sym(paste0(all_chicks[1], "_thoughts")) == "involved_in_error" ~ all_chicks[1],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[2], "_thoughts")) == "involved_in_error" ~ all_chicks[2],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[3], "_thoughts")) == "involved_in_error" ~ all_chicks[3],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[4], "_thoughts")) == "involved_in_error" ~ all_chicks[4],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[5], "_thoughts")) == "involved_in_error" ~ all_chicks[5],
                                          TRUE ~ individual_w_error))}

if(length(all_chicks) == 6){
  mmne_attr <- mmne_thoughts %>% 
    mutate(individual_w_error = case_when(individual_w_error == "other" & !!sym(paste0(all_chicks[1], "_thoughts")) == "involved_in_error" ~ all_chicks[1],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[2], "_thoughts")) == "involved_in_error" ~ all_chicks[2],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[3], "_thoughts")) == "involved_in_error" ~ all_chicks[3],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[4], "_thoughts")) == "involved_in_error" ~ all_chicks[4],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[5], "_thoughts")) == "involved_in_error" ~ all_chicks[5],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[6], "_thoughts")) == "involved_in_error" ~ all_chicks[6],
                                          TRUE ~ individual_w_error))
  fmne_attr <- fmne_thoughts %>% 
    mutate(individual_w_error = case_when(individual_w_error == "other" & !!sym(paste0(all_chicks[1], "_thoughts")) == "involved_in_error" ~ all_chicks[1],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[2], "_thoughts")) == "involved_in_error" ~ all_chicks[2],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[3], "_thoughts")) == "involved_in_error" ~ all_chicks[3],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[4], "_thoughts")) == "involved_in_error" ~ all_chicks[4],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[5], "_thoughts")) == "involved_in_error" ~ all_chicks[5],
                                          individual_w_error == "other" & !!sym(paste0(all_chicks[6], "_thoughts")) == "involved_in_error" ~ all_chicks[6],
                                          TRUE ~ individual_w_error))}
if(length(all_chicks) %in% c(3,4,5,6) == FALSE){
  message("Wrong number of chicks")
}



##################################### One chick error ########################

one_chick_error <- add_parents %>% 
  filter(father_chick == 0 & 
           mother_chick == 0 & 
           unattributable == 1 &
           no_error == length(all_chicks) - 1)

assign_missingness <- function(individual){
  added_one_col <- one_chick_error %>% 
    mutate(!!sym(paste0(individual,"_thoughts")) := case_when(!!sym(paste0(individual,"_allele1")) == "." & !!sym(paste0(individual,"_allele2")) == "." ~ "no_info",
                                                              TRUE ~ "not_missing")) %>% 
    select_at(vars(matches("locus|_thoughts")))}

missingness <- lapply(all_chicks, assign_missingness) %>% 
  reduce(full_join, by = "locus") %>% 
  left_join(one_chick_error, by = "locus")



assign_mistake_to_chick <- function(row, all_chicks){
  missingness <- missingness %>% 
    rownames_to_column() %>% 
    filter(rowname == !!row)
  
  if (length(all_chicks) == 3){
    missingness <- missingness %>% 
      mutate(individual_w_error = case_when(!!sym(paste0("m_err_code_file_", all_chicks[1])) == "unattributable" ~ all_chicks[1],
                                            !!sym(paste0("m_err_code_file_", all_chicks[2])) == "unattributable" ~ all_chicks[2],
                                            !!sym(paste0("m_err_code_file_", all_chicks[3])) == "unattributable" ~ all_chicks[3]))}
  if (length(all_chicks) == 4){
    missingness <- missingness %>% 
      mutate(individual_w_error = case_when(!!sym(paste0("m_err_code_file_", all_chicks[1])) == "unattributable" ~ all_chicks[1],
                                            !!sym(paste0("m_err_code_file_", all_chicks[2])) == "unattributable" ~ all_chicks[2],
                                            !!sym(paste0("m_err_code_file_", all_chicks[3])) == "unattributable" ~ all_chicks[3],
                                            !!sym(paste0("m_err_code_file_", all_chicks[4])) == "unattributable" ~ all_chicks[4]))}
  if (length(all_chicks) == 5){
    missingness <- missingness %>% 
      mutate(individual_w_error = case_when(!!sym(paste0("m_err_code_file_", all_chicks[1])) == "unattributable" ~ all_chicks[1],
                                            !!sym(paste0("m_err_code_file_", all_chicks[2])) == "unattributable" ~ all_chicks[2],
                                            !!sym(paste0("m_err_code_file_", all_chicks[3])) == "unattributable" ~ all_chicks[3],
                                            !!sym(paste0("m_err_code_file_", all_chicks[4])) == "unattributable" ~ all_chicks[4],
                                            !!sym(paste0("m_err_code_file_", all_chicks[5])) == "unattributable" ~ all_chicks[5]))}
  if (length(all_chicks) == 6){
    missingness <- missingness %>% 
      mutate(individual_w_error = case_when(!!sym(paste0("m_err_code_file_", all_chicks[1])) == "unattributable" ~ all_chicks[1],
                                            !!sym(paste0("m_err_code_file_", all_chicks[2])) == "unattributable" ~ all_chicks[2],
                                            !!sym(paste0("m_err_code_file_", all_chicks[3])) == "unattributable" ~ all_chicks[3],
                                            !!sym(paste0("m_err_code_file_", all_chicks[4])) == "unattributable" ~ all_chicks[4],
                                            !!sym(paste0("m_err_code_file_", all_chicks[5])) == "unattributable" ~ all_chicks[5],
                                            !!sym(paste0("m_err_code_file_", all_chicks[6])) == "unattributable" ~ all_chicks[6]))}
  if(length(all_chicks) %in% c(3,4,5,6) == FALSE) {message("Wrong number of chicks")}
  
  sum_ni <- missingness[1,] %>% str_count("no_info") %>% sum() %>% as.double()
  final_missingness <- missingness %>% 
    mutate(sum_missing = sum_ni) %>% 
    select_at(vars(matches("locus|individual_w_error|sum_missing")))
  return(final_missingness)} 


one_chick_errors_assigned <- lapply(1:nrow(missingness), assign_mistake_to_chick, all_chicks = all_chicks) %>% 
  bind_rows() %>% 
  mutate(trust = case_when(sum_missing == 0 ~ "high",
                           sum_missing == length(all_chicks) - 1 ~ "low",
                           TRUE ~ "medium"))


###########################  Next ################

unknown_error <-  add_parents %>% 
  filter(locus %in% father_mistake_and_no_error$locus == FALSE) %>%  # 749
  filter(locus %in% father_mistake_all$locus == FALSE) %>% # 699
  filter(locus %in% mother_mistake_and_no_error$locus == FALSE) %>%  #523
  filter(locus %in% mother_mistake_all$locus == FALSE) %>% #490
  filter(locus %in% one_chick_errors_assigned$locus == FALSE) %>% # 280
  mutate(individual_w_error = "hard_to_attribute",
         trust = "low")




################### Put together ###########################

father_mistake_all <- father_mistake_all %>% 
  select(locus, individual_w_error, trust)

mother_mistake_all <- mother_mistake_all %>% 
  select(locus, individual_w_error, trust)

father_mistake_some <- fmne_attr %>% 
  select(locus, individual_w_error, trust)

mother_mistake_some <- mmne_attr %>% 
  select(locus, individual_w_error, trust)

one_chick_error <- one_chick_errors_assigned %>% 
  select(locus, individual_w_error, trust)

unknown_error <- unknown_error %>% 
  select(locus, individual_w_error, trust)

all_of_the_blame <- bind_rows(father_mistake_all, mother_mistake_all, father_mistake_some, mother_mistake_some, one_chick_error, unknown_error) 


table(all_of_the_blame$individual_w_error)
table(all_of_the_blame$trust)

high_quality_blame <- all_of_the_blame %>% 
  filter(individual_w_error != "hard_to_attribute" & trust == "high") %>% 
  select(-trust)

table(high_quality_blame$individual_w_error)

write_delim(high_quality_blame, path = paste0(me_vcf_path, "error_attributions_family_", family_num, ".txt"))



####################################################
# Re-add a single genotype col for each individual #
####################################################

onecol_alpha_alleles <- function( name, df) {
   mutate(df, !!paste0(name, "_combined") := str_c(!!sym(paste0(name, "_allele1")), !!sym(paste0(name, "_allele2")), sep = "/")) %>% 
    select_at(vars(matches("locus|rowname|_combined")))}


combined_letter_alleles <- lapply(all_of_the_family, onecol_alpha_alleles, df = letter_alleles) %>% 
  reduce(full_join, by = c("locus", "rowname")) %>% 
  rename_at(vars(ends_with("_combined")),
            funs(str_remove(., "_combined"))) %>% 
  gather(key = indiv, value = genotype, -rowname, -locus)


locus_indiv_gt <- inner_join(high_quality_blame, combined_letter_alleles, by = c(locus = "locus", individual_w_error = "indiv")) %>% 
  select(-rowname) %>% 
  filter(individual_w_error != mother) %>% 
  filter(individual_w_error != father)


write_delim(locus_indiv_gt, path = paste0(me_vcf_path, "chick_error_attributions_genotypes_family_", family_num, ".txt"))


#############################################
# Add a col for what the genotype should be #
#############################################


loci_of_interest <- locus_indiv_gt %>% 
  select(locus) %>% 
  pull()

parents_letters <- letter_alleles %>% 
  rename_at(vars(matches(father)),
            funs(str_c(., "_father", sep = ""))) %>% 
  rename_at(vars(matches(mother)),
            funs(str_c(., "_mother", sep = "")))

make_indiv_possibilities <- function(df, mother, father){
  df %>% 
    mutate(chick_possibility1 = str_c(!!sym(paste0(mother,"_allele1_mother")), !!sym(paste0(father,"_allele1_father")), sep = "/")) %>%
    mutate(chick_possibility2 = str_c(!!sym(paste0(mother,"_allele1_mother")), !!sym(paste0(father,"_allele2_father")), sep = "/")) %>%
    mutate(chick_possibility3 =  str_c(!!sym(paste0(mother,"_allele2_mother")), !!sym(paste0(father,"_allele1_father")), sep = "/")) %>%
    mutate(chick_possibility4 =  str_c(!!sym(paste0(mother,"_allele2_mother")), !!sym(paste0(father,"_allele2_father")), sep = "/")) %>% 
    select_at(vars(matches("locus|rowname|possibility")))
}

  
chick_possibilities <- make_indiv_possibilities(df = parents_letters, mother = mother, father = father) %>% 
   filter(locus %in% loci_of_interest) %>% 
   mutate(multiple_possibilities = case_when(chick_possibility1 == chick_possibility2 & 
                                               chick_possibility1 == chick_possibility3 & 
                                               chick_possibility1 == chick_possibility4 ~ "all_same",
                                             TRUE ~ "differences"))
  
correct_gt <- full_join(locus_indiv_gt, chick_possibilities, by = "locus") %>% 
  filter(multiple_possibilities == "all_same") %>% 
  mutate(missing = case_when(str_detect(chick_possibility1, "\\.") == TRUE ~ "missing",
                             TRUE ~ "non_missing")) %>% 
  filter(missing == "non_missing") %>% 
  select(locus, individual_w_error, chick_possibility1)


write_delim(correct_gt, paste0(me_vcf_path, "chick_error_attributions_correct_gt_family_", family_num, ".txt"))



# single_col_letter_alleles <- letter_alleles %>% 
#   mutate(!!sym(paste0(all_of_the_family, "_alleles")) := str_c(Wendy_allele1,Wendy_allele2, sep = "/"))
# left_join(high_quality_blame , flat_var_gq, by = c("locus", "individual")
# all_of_the_family

#################################
# Section 3: Genotype qualities #
#################################

# 
# #  SLOW old version
# # multisib_family_GQs <- lapply(multisib_single_family$random_group, function(group_id) {
# #   read_delim(file = paste0(me_gq_path, "part_", group_id, "_keep_indivs.GQ.FORMAT"), 
# #              delim = '\t') %>% 
# #     mutate(locus = str_c(CHROM, POS, sep = ":")) %>%
# #     filter(locus %in% me_sites) %>%
# #     select(-CHROM, -POS)})
# 
# # Read in GQ info and keep only for sites which had mendelian errors.
# multisib_family_GQs <- lapply(multisib_single_family$random_group, function(group_id) {
#   gqs_to_import <- paste0(me_gq_path, "part_", group_id, "_keep_indivs.GQ.FORMAT")
# 
#   preprocessing_command <- sprintf("grep --file %s %s",
#                                    paste0("me_locus_file_", family_num, ".txt"),
#                                    paste(gqs_to_import, collapse = " "))
# 
#   gqs <- data.table::fread(cmd = preprocessing_command) %>%
#     mutate(locus = str_c(CHROM, POS, sep = ":")) %>%
#     filter(locus %in% me_sites) %>%
#     select(-CHROM, -POS) %>%
#     rownames_to_column()
#   return(gqs)})
# 
# 
# # Flatten, remove duplicated columns (the parents have been read in multiple times)
# flat_gq <- multisib_family_GQs %>%
#   reduce(full_join, by = c("locus", "rowname")) %>% 
#   select(-ends_with(".x.x"), -ends_with(".y")) %>% 
#   rename_at(vars(ends_with(".x")), funs(str_remove(., "\\.x"))) %>% 
#   filter(locus %in% me_sites) %>% 
#   unique() 
# 
# # There are duplicated rows in the GQ file. I only want to keep the one that has the most non-missing GQ information
# # Currently do this by counting the number of dots ("."; missing data) and picking the entry with the lowest number of dots
# # This would miss info for an individual that is in one row but not the row that has the least missing data. 
# # Afaik there are no cases of this but I should possibly do a more thorough check.
# added_counting <- flat_gq %>%
#   unite(for_counting, matches("[^locus|^rowname]" ), sep = "") %>% 
#   mutate(num_dots = str_count(for_counting, "\\.")) %>% 
#   filter(locus %in% extra_difficulty == TRUE)
# 
# ranked_flat_gq_data <- full_join(flat_gq, added_counting, by = c(locus = "locus", rowname = "rowname")) %>% 
#   group_by(locus) %>% 
#   top_n(n = -1, wt = num_dots) %>% 
#   select(-for_counting, -rowname, -num_dots)
# 
# 
# # CHECK. Does the GQ data and the ME data have the same number of rows? NO, because of some problem loci which I am dealing with later.
# if (nrow(ranked_flat_gq_data) == nrow(joined_me_codes_summary)) {
#   message("The dataset has the correct number of rows")
# } else {
#     message("Something has gone wrong! The datasets have different numbers of rows")}
# 
# tidied_gq_data <- ranked_flat_gq_data %>%
#   rename_at(vars(contains(father)), funs(str_c(., '_(father)'))) %>% 
#   rename_at(vars(contains(mother)), funs(str_c(., '_(mother)'))) 
# 
# 
# # End of section, remove excess stuff from environment
# rm(multisib_family_GQs, flat_gq, added_counting, ranked_flat_gq_data)

#######################
# Section 4: Variants #
#######################

# # Import vcf lines for the mendelian error variants
# multisib_family_variants <- lapply(multisib_single_family$random_group, function(group_id) {
#   vcf_to_import <- paste0(me_vcf_path,"part_", group_id, "_keep_indivs.vcf")
#  
#   preprocessing_command <- sprintf("grep --file %s %s", 
#                                  paste0("me_locus_file_", family_num, ".txt"), 
#                                  paste(vcf_to_import, collapse = " "))
# 
#   vcf <- data.table::fread(cmd = preprocessing_command) %>%
#     select(-ID, -QUAL, -FILTER, -INFO, -FORMAT) %>% 
#     rename(CHROM = `#CHROM`) %>% 
#     rename_at(vars(-matches("CHROM|POS|REF|ALT")),
#             funs(str_c(., "_indiv", sep = ""))) %>% 
#     mutate_at(vars(ends_with("_indiv")),
#             funs(str_extract(., "^[:graph:]{3}(?=:)"))) %>%
#     mutate(locus = str_c(CHROM, POS, sep = ":")) %>%
#     select(-CHROM, -POS) %>%
#     rownames_to_column()
#   return(vcf)})
# 
# # Flatten, remove duplicated columns (the parents have been read in multiple times), remove duplicated rows
# # Why does filter locus %in% me sites remove any loci????????????????????????
# flat_var <- multisib_family_variants %>%
#   reduce(full_join, by = c("locus", "rowname")) %>% 
#   select(-ends_with(".x.x"), -ends_with(".y")) %>% 
#   rename_at(vars(ends_with(".x")),
#             funs(str_remove(.,"\\.x"))) 
# 
# separate_alleles <- function(tibble, indiv_name) {  
#     separate(data = tibble, col = paste0(indiv_name, "_indiv"), into = c(paste0(indiv_name, "_allele1"), paste0(indiv_name, "_allele2")), sep = "[/\\|]") %>% 
#       select_at(vars(matches("_allele[12]|locus|rowname")))
#   }
#   
# separated_alleles <- lapply(multisib_single_family$ID, function(individual){
#     separate_alleles(flat_var, individual)
#   }) %>% 
#   reduce(full_join, by = c("locus", "rowname")) 
#     
# sep_alleles_info <- full_join(flat_var, separated_alleles, by = c("locus", "rowname")) %>% 
#   select_at(vars(-ends_with("_indiv")))
# 
# # Find the highest alt allele. Could use this to determine how many columns to make.
# # Currently just using this to 
# ### FIND A BETTER WAY ####
# ### Apparently this doesn't even work ####
# allele_cols_names <-  sep_alleles_info %>% 
#     select_at(vars(matches("_allele[12]"))) %>% 
#     names()
#     
# locus_max <-   sep_alleles_info %>%  mutate(most = pmax(!!!rlang::syms(allele_cols_names)))
# most_alts <- as.numeric(max(locus_max$most))
# 
# 
# # Select one line for each locus based on having the least missingness
# added_counting_vcf <- sep_alleles_info %>%
#   unite(for_counting, matches("[^locus|^rowname|^REF|^ALT]" ), sep = "") %>% 
#   mutate(num_dots = str_count(for_counting, "\\.")) 
# 
# #################### NEED TO THINK MORE ABOUT WHY THERE ARE EXTRA LOCI ^^ HERE THAT SHOULDN'T BE ######  
# # Is grep treating the two values I give it as an or?
# 
# numeric_allele_rep <- full_join(sep_alleles_info, added_counting_vcf, by = c("locus", "rowname", "REF", "ALT")) %>% 
#   filter(locus %in% me_sites == TRUE) %>% 
#   group_by(locus) %>% 
#   top_n(n = -1, wt = num_dots) %>% 
#   select(-for_counting, -rowname, -num_dots) 
# 
# 
# # This is also still slightly longer than it should be because there are some duplicated loci, see:
# problem_loci_names <- table(numeric_allele_rep$locus) %>% 
#   as.data.frame() %>% 
#   as_tibble() %>% 
#   filter(Freq > 1) %>% 
#   select(Var1) %>% 
#   pull()
# 
# ## Fill the column with the actual alleles that individual has.
# # Some of the unused alts will be dropped. Need to double check this is working properly
# alt_allele_cols <- numeric_allele_rep %>% 
#   separate(col = ALT, 
#            into = paste0("alt_", 1:most_alts),
#            sep = ",",
#            fill = "right",
#            extra = "drop") %>% 
#   group_by(locus) 
# 
# 
# 
#  
# # Function to convert numeric allele representation into alphabetical representation
# alpha_alleles <- function(df, name) {
#   mutate(df, !!name := case_when(!!sym(name) == 0 ~ REF,
#                                  !!sym(name) == 1 ~ alt_1,
#                                  !!sym(name) == 2 ~ alt_2,
#                                  !!sym(name) == 3 ~ alt_3,
#                                  !!sym(name) == 4 ~ alt_4,
#                                  !!sym(name) == 5 ~ alt_5,
#                                  !!sym(name) == "." ~ ".",
#                                  TRUE ~ "dropped"))}
#   
#   
# letter_alleles <- lapply(allele_cols_names, function (col_name) {
#   alt_allele_cols %>% 
#     select(locus, REF, alt_1, alt_2, alt_3, alt_4, alt_5, col_name) %>% 
#     alpha_alleles(col_name)}) %>% 
#   reduce(full_join, by = c("locus", "REF", "alt_1","alt_2", "alt_3", "alt_4", "alt_5")) %>% 
#   select(-REF, -alt_1, -alt_2, -alt_3, -alt_4, -alt_5)
# 
# # Check if any alternate alleles have been dropped
# if(any(letter_alleles=="dropped")){message("Alleles have been dropped. Fix this")}
# 
# # End of section, remove excess stuff from environment
# rm(multisib_family_variants, flat_var, separated_alleles, sep_alleles_info, locus_max, added_counting_vcf, numeric_allele_rep,  alt_allele_cols)











# 
# check_repe group_by(locus) %>% 
#   count() %>% 
#   filter(n > 1) %>% 
#   select(locus) %>% 
#   pull()
# 
# 
# father_mistake_all %>% filter(locus %in% reconstituted_smaller == TRUE) # 50
# father_mistake_some %>% filter(locus %in% reconstituted_smaller == TRUE) #251
# one_chick_error %>% filter(locus %in% reconstituted_smaller == TRUE) #0
# unknown_error %>% filter(locus %in% reconstituted_smaller == TRUE) #510



































# assign_allele(individual = "Ruapuke",
#               parent_name = mother,
#               parent = "mother",
#               list_df = mmne,
#               row_or_sublist = 14,
#               all_chicks = all_chicks)
# assign_allele <- function(individual, parent_name, parent, list_df, row_or_sublist) {
#   
# single_row <- list_df[[row_or_sublist]] %>% 
#     mutate(!!sym(paste0(individual[1],"_thoughts")) := case_when(!!sym(paste0(individual[1],"_allele1")) == "." & !!sym(paste0(individual[1],"_allele2")) == "." ~ "no_info",
#                                                               str_detect(!!sym(paste0(individual[1],"_allele1")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
#                                                               str_detect(!!sym(paste0(individual[1],"_allele2")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
#                                                               str_detect(!!sym(paste0(individual[1],"_allele1")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
#                                                               str_detect(!!sym(paste0(individual[1],"_allele2")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
#                                                               TRUE ~ "do_not_match"))
# if (length(individual) == 2){
#   single_row <- single_row %>% 
#     mutate(!!sym(paste0(individual[2],"_thoughts")) := case_when(!!sym(paste0(individual[2],"_allele1")) == "." & !!sym(paste0(individual[2],"_allele2")) == "." ~ "no_info",
#                                                                  str_detect(!!sym(paste0(individual[2],"_allele1")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
#                                                                  str_detect(!!sym(paste0(individual[2],"_allele2")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
#                                                                  str_detect(!!sym(paste0(individual[2],"_allele1")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
#                                                                  str_detect(!!sym(paste0(individual[2],"_allele2")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
#                                                                  TRUE ~ "do_not_match"))}
# if (length(individual) == 3){
#   single_row <- single_row %>% 
#     mutate(!!sym(paste0(individual[3],"_thoughts")) := case_when(!!sym(paste0(individual[3],"_allele1")) == "." & !!sym(paste0(individual[3],"_allele2")) == "." ~ "no_info",
#                                                                  str_detect(!!sym(paste0(individual[3],"_allele1")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
#                                                                  str_detect(!!sym(paste0(individual[3],"_allele2")), !!sym(paste0(parent_name,"_allele1_", parent))) == TRUE ~ "could_be_correct",
#                                                                  str_detect(!!sym(paste0(individual[3],"_allele1")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
#                                                                  str_detect(!!sym(paste0(individual[3],"_allele2")), !!sym(paste0(parent_name,"_allele2_", parent))) == TRUE ~ "could_be_correct",
#                                                                  TRUE ~ "do_not_match"))}
# if (length(individual) > 3){ message("Need to add possibility of more having more chicks to the function")}
# 
# sum_cbc <- single_row[1,] %>%  str_count("could_be_correct") %>% sum() %>% as.double()
# 
# sum_ni <- single_row[1,] %>% str_count("no_info") %>% sum() %>% as.double()
# 
# single_row <- single_row %>% 
#   mutate(gt_possible_counts = sum_cbc,
#          gt_no_info = sum_ni) 
# return(single_row)
# }
# 
# assign_allele(individual = "Ruapuke",
#               parent_name = mother,
#               parent = "mother",
#               list_df = mmne,
#               row_or_sublist = 14)
# 
# 
# assess_no_error_chicks <- function (row_or_sublist, list_df, parent_name, parent) {  
#   got_names <- list_df[[row_or_sublist]] %>% 
#     select_at(vars(ends_with("_allele1"))) %>% 
#     rename_at(vars(ends_with("_allele1")),
#             funs(str_remove(.,"_allele1"))) %>% 
#   names()
#   
#   trial <- assign_allele(got_names, parent_name, parent,  list_df, row_or_sublist)
#   return(trial)
# }
# 
# mmne_thoughts <- lapply(1:length(mmne), assess_no_error_chicks, list_df = mmne, parent_name = mother, parent = "mother") %>% 
#   bind_rows() %>% 
#   mutate_at(vars(ends_with("_thoughts")), funs(str_replace_na(., replacement = "involved_in_error")))
# 
# fmne_thoughts <- lapply(1:length(fmne), assess_no_error_chicks, list_df = fmne, parent_name = father, parent = "father") %>% 
#   bind_rows() %>% 
#   mutate_at(vars(ends_with("_thoughts")), funs(str_replace_na(., replacement = "involved_in_error"))) %>% 
#   select_at(vars(matches("locus|_thoughts|gt_")))

# 
# ###################### Figure out what still needs doing ##########################
# 
# # so far have:
# father_mistake_and_no_error <- add_parents %>% 
#   filter(father_chick > 0 & 
#            mother_chick == 0 & 
#            unattributable == 0 &
#            no_error > 0) %>% 
#   select(locus) %>% 
#   pull()
# 
# father_mistake_all <- add_parents %>% 
#   filter(father_chick > 0 & 
#            mother_chick == 0 & 
#            unattributable == 0 &
#            no_error == 0) %>% 
#   mutate(individual_w_error = father) %>% 
#   mutate(chick_gq_sum = select(., ends_with("_gq")) 
#          %>% rowSums()) %>% 
#   select(locus) %>% 
#   pull() 
# 
# mother_mistake_and_no_error <- add_parents %>% 
#   filter(father_chick == 0 & 
#            mother_chick > 0 & 
#            unattributable == 0 &
#            no_error > 0) %>% 
#   select(locus) %>% 
#   pull()
# 
# mother_mistake_all <- add_parents %>% 
#   filter(father_chick == 0 & 
#            mother_chick > 0 & 
#            unattributable == 0 &
#            no_error == 0) %>% 
#   mutate(individual_w_error = mother) %>% 
#   mutate(chick_gq_sum = select(., ends_with("_gq")) 
#          %>% rowSums()) %>% 
#   select(locus) %>% 
#   pull()
# 
# 
# 
# # Next
# chick_error <-  add_parents %>% 
#   filter(locus %in% father_mistake_and_no_error == FALSE) %>%  # 749
#   filter(locus %in% father_mistake_all == FALSE) %>% # 699
#   filter(locus %in% mother_mistake_and_no_error == FALSE) %>%  #523
#   filter(locus %in% mother_mistake_all == FALSE) #490



# if(length(all_chicks) == 3){
#   blame <- missingness %>% 
#       mutate(individual_w_error = case_when(!!sym(paste0("m_err_code_file_", all_chicks[1])) == "unattributable" ~ all_chicks[1],
#                                             !!sym(paste0("m_err_code_file_", all_chicks[2])) == "unattributable" ~ all_chicks[2],
#                                             !!sym(paste0("m_err_code_file_", all_chicks[3])) == "unattributable" ~ all_chicks[3])) %>% 
#       mutate(trust = case_when(individual_w_error == all_chicks[1] & 
#                                  !!sym(paste0(all_chicks[2], "_thoughts")) == "not_missing" &
#                                  !!sym(paste0(all_chicks[3], "_thoughts")) == "not_missing" ~ "high",
#                                individual_w_error == all_chicks[2] & 
#                                  !!sym(paste0(all_chicks[1], "_thoughts")) == "not_missing" &
#                                  !!sym(paste0(all_chicks[3], "_thoughts")) == "not_missing" ~ "high",
#                                individual_w_error == all_chicks[3] & 
#                                  !!sym(paste0(all_chicks[2], "_thoughts")) == "not_missing" &
#                                  !!sym(paste0(all_chicks[1], "_thoughts")) == "not_missing" ~ "high",
#                                TRUE ~ "low"))}





# #noquote(paste0("individual_w_error == 'other' & ", rep_len(all_chicks, length.out = length(all_chicks)), "_thoughts == 'involved_in_error' ~ ", rep_len(all_chicks, length.out = length(all_chicks))))
# 
# dsf <- paste0("individual_w_error == 'other' & ", 
#        rep_len(all_chicks, length.out = length(all_chicks)),
#        "_thoughts == 'involved_in_error' ~ ", 
#        rep_len(all_chicks, length.out = length(all_chicks)), collapse = ", ")
# 
# 
# cols <- c("top_speed" = NA_real_ , "nhj" = NA_real_ ,"mpg" = NA_real_)
# 
# add_column(mtcars, !!!cols[setdiff(names(cols), names(mtcars))])
# 
# 
# # mutate(individual_w_error = case_when(individual_w_error == "other" & !!sym(paste0(all_chicks[1], "_thoughts")) == "involved_in_error" ~ all_chicks[1],
# #                                       individual_w_error == "other" & !!sym(paste0(all_chicks[2], "_thoughts")) == "involved_in_error" ~ all_chicks[2],
# #                                       individual_w_error == "other" & !!sym(paste0(all_chicks[3], "_thoughts")) == "involved_in_error" ~ all_chicks[3],
# #                                       TRUE ~ individual_w_error)) %>% 
# # mutate(individual_w_error = case_when(!!syms(paste0("individual_w_error == 'other' & ", 
# #                                                    rep_len(all_chicks, length.out = length(all_chicks)),
# #                                                    "_thoughts == 'involved_in_error' ~ ", 
# #                                                    rep_len(all_chicks, length.out = length(all_chicks)),
# #                                                    ",\n")),
# #                                       TRUE ~ "no")) %>% 
# mutate(individual_w_error = case_when(!!sym(dsf),
#                                       TRUE ~ "no"))
# 
# 
# fmne_father <- fmne_thoughts %>% 
#   filter(individual_w_error == father)
# 
# 
# 
# library(tidyverse)
# mtcars %>%
#   tbl_df() %>%
#   rownames_to_column("car") %>%
#   rowwise() %>%
#   mutate(top_speed = ifelse("top_speed" %in% names(.), top_speed, NA),
#          mpg = ifelse("mpg" %in% names(.), mpg, NA)) 
# 
# 
# 
# 
# 
# 
# # 
# dfnhs <- function(df, all_the_chicks, father_char) {
#   fmne_thoughts <- lapply(1:length(df), assess_no_error_chicks, list_df = df, parent_name = father_char, parent = "father", all_chicks = all_the_chicks) %>% 
#     bind_rows() %>% 
#     mutate_at(vars(ends_with("_thoughts")), funs(str_replace_na(., replacement = "involved_in_error"))) %>% 
#     select_at(vars(matches("locus|_thoughts|gt_"))) %>% 
#     mutate(individual_w_error = case_when(gt_error > gt_possible_counts ~ father_char,
#                                           gt_error == 1 & gt_possible_counts > 1 ~ "other",
#                                           TRUE ~ "mystery")) %>% 
#     mutate(trust = case_when(individual_w_error == father & gt_no_info < gt_error ~ "highish",
#                              individual_w_error == father & gt_no_info > gt_error ~ "lowish",
#                              individual_w_error == "other" & gt_no_info < gt_possible_counts ~ "highish",
#                              individual_w_error == "other" & gt_no_info > gt_possible_counts ~ "lowish",
#                              TRUE ~ "low"))  %>% 
#     mutate(individual_w_error = case_when(!!sym(paste0("individual_w_error == 'other' & ", 
#                                                        rep_len(all_the_chicks, length.out = length(all_the_chicks)),
#                                                        "_thoughts == 'involved_in_error' ~ ", 
#                                                        rep_len(all_the_chicks, length.out = length(all_the_chicks)), collapse = ", ")),
#                                           TRUE ~ "no"))
# }
# 
# 
# dfnhs(fmne, all_chicks, father)










# 
# 
# 
# 
# cols <- c(top_speed = NA_real_, nhj = NA_real_, mpg = NA_real_)
# 
# add_column(mtcars, !!!cols[setdiff(names(cols), names(mtcars))])
# 
# 
# assign_allele(got_names, mother, "mother",  mmne, 16)
# 
# do_a_thing(16, mmne)
# dsafkj <- onsrfe %>% slice(1) %>% unlist() %>% unname() %>%  str_count("could_be_correct") %>% sum()
# dsflh <- onsrfe %>% slice(1) %>% unlist() %>% unname() %>%  str_count("no_info") %>% sum()
# 
# 
# onsrfe %>% mutate(dkasdlu = dsafkj,
#                   sdf = dsflh)
# 
# 
# 
# got_names <- mmne[[16]] %>% 
#   select_at(vars(ends_with("_allele1"))) %>% 
#   rename_at(vars(ends_with("_allele1")),
#             funs(str_remove(.,"_allele1"))) %>% 
#   names()
# check_error_allele_same <- function(){}
# 
# potentially_incorrect_allele <- #Think maybe I need this??
# 
# 
# individuall <- c("Hananui", "Ruapuke")
# parent_name <- "Lisa"
# parent <- "mother"
#   peo <- "Hananui"
# 
#  
#  
# 
#  nm1 <- grep("no_error", colnames(mother_mistake_and_no_error), value = TRUE)
#  
#  mother_mistake_and_no_error %>% 
#    rename_at(vars(nm1), funs(gsub("mean", "", nm1)))
#  
#  
# no_error_indiv <- 
# 
#   paste(rep(vars, each = length(vis)), vis, sep = ".")
# 
# 
# 
# # if (mother_mistake_and_no_error[1]$paste0("m_err_code_file_", individual) == "no_error") {
# #   blahh <- !!paste0("m_err_code_file_", individual)}{
# #     blahh <- paste0("not_", individual)}
# # 
# 
# 
# 
# does_it_support <- function (individual, row) {
#   blah <- mother_mistake_and_no_error %>% 
#     rownames_to_column() %>% 
#     filter(rowname == row)
#   blahh <- vector()
#   
# }
# dsfkjgs <- does_it_support("Ruapuke", 1)
# 
# outer(seq_along(list), seq_along(list), 
#       FUN= Vectorize(function(i,j) add_function(list[[i]], list[[j]])))
# 
# 
# chicks_w_no_error <- mother_mistake_and_no_error %>% 
#   filter_all(any_vars(. %in% 'no_error'))
#   
#   
#   
#   
#   
#   %>% 
#   select_if(m_err_code_file_Ruapuke == "no_error" )
# 
# by(mother_mistake_and_no_error, 1:nrow(mother_mistake_and_no_error), function(row) {
#   dostuff)
# 
# 
# df %>% filter_all(any_vars(. %in% c('M017', 'M018')))
# 
# 
# 
# 
# mmne <- does_it_support("Ruapuke")
# #### REMEMBER Filter off anything that is on the x and y chromosomes






# asdfa <- smaller %>%
#   rename(`13` = m_err_code_file_13) %>% 
#   rename_at(vars(matches("m_err_code_file_")),
#             funs(lookup_rename(., for_renaming)))
# 
# 
# select(m_err_code_file_13) %>%
#   rename(!!paste0("m_err_code_file_", names(for_renaming)) := paste0("m_err_code_file_", for_renaming))



# asdkjg <- lookup_rename(asdfa, for_renaming)
# 
# names(df2)[names(df2) %in% for_renaming$old] = vm$new[match(names(df2)[names(df2) %in% vm$old], vm$old)]
# 
# renamer <- function(df, old_name, new_name) {
#   rename(df, !!old_name := !!new_name)}
# 
# asdfasfftesw <- lapply(for_renaming, function (named_vector) {
#   older <- paste0("m_err_code_file_", named_vector)
#   newer <- paste0("m_err_code_file_", names(named_vector))
#   print(names(named_vector))
#   smaller %>% 
#     renamer(older, newer)
# })
# 
#   


# r_alleles <- lapply(allele_cols_names, function (col_name) {
#   alt_allele_cols %>% 
#     select(locus, rowname, REF, alt_1, alt_2, alt_3, alt_4, alt_5, col_name) %>% 
#     alpha_alleles(col_name)}) %>% 
# reduce(full_join, by = c("locus", "rowname", "REF", "alt_1","alt_2", "alt_3", "alt_4", "alt_5")) %>% 
# select(-REF, -alt_1, -alt_2, -alt_3, -alt_4, -alt_5)
# 
# 
# rename_at(vars(starts_with("m_err_code_file_")),
#           funs(case_when(. == paste)))



# # gq_to_merge$locus %>% unique() %>% length()
# # [1] 18643
# # Current length 18648
# gq_to_merge <- tidied_gq_data %>% 
#   filter(locus %in% problem_loci_names == FALSE) 
# 
# extra_difficulty <- variants_to_merge %>% 
#   group_by(locus) %>% 
#   count()
# 
# # variants_to_merge$locus %>% unique() %>% length()
# # [1] 18643
# # Current length 18643. Also correct
# variants_to_merge <- letter_alleles %>% 
#   filter(locus %in% me_sites == TRUE) %>% 
#   filter(locus %in% problem_loci_names == FALSE) %>% # Longer
#   select(locus) %>% 
#   unique() #%>% 
#   filter(locus %in% me_to_merge == FALSE)


##### Things to do:
# Figure out why the resulting tbls are not the same length
# Change the _file_[:digit:] part of the 
# Figure out how to assign mistake to an individual - maybe the ones that are not errors, can say which indiv they support







####### Can start from here if I haven't deleted the file ######
# write_delim(numeric_allele_rep, "num_allele_temp.txt", delim = "\t")
# 
# most_alts <-5
# alpha_allele_rep <-  read_delim("num_allele_temp.txt", delim = "\t") %>% 
#   separate(col = ALT, 
#            into = paste0("alt_", 1:most_alts),
#            sep = ",",
#            fill = "right",
#            extra = "drop") %>% 
#   group_by(locus) 
#################################################################


#####################
# BROKEN, FIX THESE #
#####################
# 
# Should come after the M.E file write out
# # CHECK. There should be no loci where the each parent has different genotypes recorded
# # Genotypes are recorded as */* if there will be an error regardless of what they are. 
# # This can result in parent records appearing different even if they are the same, so I'm filtering these out
# 
# # This should be 0 (no differences in the same individual) 
# # Currently just doing a pairwise comparison between results from two random files
#  parent1s_to_compare <- paste0("parent1_file_", sample(relevant_file_nums, 2, replace = FALSE))
#   
# check_same_parent1 <- all_me_flat %>% 
#   select(parent1s_to_compare[1], parent1s_to_compare[2]) %>% 
#   filter(is.na(parent1s_to_compare[1] == FALSE)) %>% 
#   filter(is.na(parent1s_to_compare[2] == FALSE)) %>%
#   filter(parent1s_to_compare[1] != parent1s_to_compare[2]) %>% 
#   filter(parent1s_to_compare[1] != "*/*" & parent1s_to_compare[2] != "*/*")
# message(paste("There are", nrow(check_same_parent1), "differences between the columns corresponding to the same parent. There should be 0 differences."))
# 
# # This should be 0 (no differences in the same individual)
# parent2s_to_compare <- paste0("parent2_file_", sample(relevant_file_nums, 2, replace = FALSE))
# 
# check_same_parent2 <- all_me_flat %>% 
#   select(parent2s_to_compare[1], parent2s_to_compare[2]) %>% 
#   filter(is.na(parent2s_to_compare[1] == FALSE)) %>% 
#   filter(is.na(parent2s_to_compare[2] == FALSE)) %>%
#   filter(parent2s_to_compare[1] != parent2s_to_compare[2]) %>% 
#   filter(parent2s_to_compare[1] != "*/*" & parent2s_to_compare[2] != "*/*")
# message(paste("There are", nrow(check_same_parent2), "differences between the columns corresponding to the same parent. There should be 0 differences."))
# 
# # This should NOT be 0 (means this test is bad)
# chicks_to_compare <- paste0("chick_file_", sample(relevant_file_nums, 2, replace = FALSE))
# 
# check_different_chicks <- all_me_flat %>% 
#   select(chicks_to_compare[1], chicks_to_compare[2]) #%>% 
#   filter(is.na(chicks_to_compare[1] == FALSE)) #%>% 
#   filter(is.na(chicks_to_compare[2] == FALSE))# %>%
#   filter(chicks_to_compare[1] != chicks_to_compare[2]) %>% 
#   filter(chicks_to_compare[1] != "*/*" & chicks_to_compare[2] != "*/*")
# message(paste("There are", nrow(check_different_chicks), "differences between the offspring columns. There should NOT be 0 differences."))
# 


###############
# IN PROGRESS #
###############


#######
# OLD #
#######

# # THIS NEEDS TO BE LENGTH O. OTHERWISE THERE WILL BE TROUBLE!
# bound_correctly <- flat %>%
#   filter(locus != locus1)

# Probably best to pass the family ID from outside the script, then run the script once for each one 

# run_trios <-  read_delim("/cifs/biocldap/student_users/marissalelec/projects/kakapo-test/output/07_offspring_m_err/all_true_trio_parents.txt", delim = "\t", trim_ws = TRUE, col_names = FALSE)
# ajd <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/07_offspring_m_err/multitrio_files.txt", delim = ' ', col_names = F)


# mendel_file <- "/media/drive_6tb/projects/kakapo-genomics/output/07_offspring_m_err/all_details_1.mendel"
# GQ_file <- "/media/drive_6tb/projects/kakapo-genomics/output/07_offspring_m_err/part_1_keep_indivs.GQ.FORMAT"


# lapply(names(nested_input_files), function(data_subset){
#   subset_files <- nested_input_files[[data_subset]]
#   write_delim(as_tibble(subset_files[[1]]), paste0("part_",data_subset,"_keep_indivs.txt"), col_names = FALSE)
#   write_delim(as_tibble(subset_files[[2]]), paste0("part_",data_subset,"_names.txt"), col_names = FALSE)
#   write_delim(as_tibble(subset_files[[3]]), paste0("part_",data_subset,"_parents.txt"), col_names = FALSE)
#   write_delim(as_tibble(subset_files[[4]]), paste0("part_",data_subset,"_update_sex.txt"), col_names = FALSE)})

# 
# all_me_in_multisib_family <- lapply(dskjfg$random_group, function(group_id) {
#   read_delim(file = paste0("/media/drive_2tb/projects/kakapo-genomics/all_details_", group_id, ".lmendel"), 
#              delim = ' ', trim_ws = TRUE) %>% 
#     filter(N != 0)})

# lapply(dskjfg$random_group, read_delim(file = paste0("/media/drive_2tb/projects/kakapo-genomics/part_", dskjfg$random_group[], "_keep_indivs.GQ.FORMAT"), delim = '\t'))



# blah <- read_delim(mendel_file, delim = ' ', col_names = FALSE, skip = 1, trim_ws = TRUE)
# blahblah <- read_delim(GQ_file, delim = '\t') 
# 
# mendel_1 <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/07_offspring_m_err/all_details_1.mendel", delim = ' ', col_names = FALSE, skip = 1, trim_ws = TRUE) %>% 
#   select(X4)
# # mendel_30 <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/07_offspring_m_err/all_details_30.mendel", delim = ' ', col_names = FALSE, skip = 1, trim_ws = TRUE) %>% 
# #   select(X4)
# mendel_70 <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/07_offspring_m_err/all_details_70.mendel", delim = ' ', col_names = FALSE, skip = 1, trim_ws = TRUE) %>% 
#   select(X4)
# 
# 
# all_me_sites <- bind_rows(mendel_1, mendel_70) %>% 
#   unique() %>% 
#   pull()
# 
# 


