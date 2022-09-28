#!/usr/bin/env Rscript

library(tidyverse)

#############
# FUNCTIONS #
#############

recode_genos <- function(indiv, vcf){
  vcf %>% 
    select(all_of(indiv)) %>% 
    rownames_to_column() %>% 
    mutate(!!sym(indiv) := str_extract(!!sym(indiv),"[:graph:]{3}(?=:)")) %>% 
    mutate(!!paste0(sym(indiv),"_num_rep") := case_when(str_detect(!!sym(indiv), "\\.") == TRUE ~ "5",
                                                        !!sym(indiv) == "0/0" ~ "0",
                                                        !!sym(indiv) == "1/1" ~ "2",
                                                        !!sym(indiv) == "0/1" ~ "1",
                                                        !!sym(indiv) == "1/0" ~ "1",
                                                        TRUE ~ "need_more_options")) %>% # I shouldn't have cases of this and I don't when running from within R, but I do when running from the nf script (last genotype) is wrecked
    

    select(-!!sym(indiv))
  
}

# recode_section <- function(vcf_file, indiv_names){
#   
# }


###########
# GLOBALS #
###########

files <- commandArgs(trailingOnly = TRUE)

chromosome <- files[1]
bird_data_file <- files[2]

vcf_name <- paste0("filtered_chr_", chromosome, ".vcf")

# # Dev
# chromosome <- "CM013774.1"
# bird_data_file <- "data/bird_info_rships.csv"
# chromosome <- "CM013761.1"
# vcf_name <- paste0("/media/drive_6tb/projects/kakapo-genomics/output/10_fimpute_phasing/filtered_chr_", chromosome, ".vcf")


header_rows <- 113

########
# MAIN #
########


length_vcf <- length(readLines(vcf_name))
message("Length of vcf")
message(length_vcf)

n_snps_chr <- length_vcf-header_rows

message(n_snps_chr)


vcf_file <- read_delim(vcf_name, delim = '\t', skip = header_rows, col_names = FALSE) %>% 
      select(-num_range("X", 3:9))


# Make snp_info files
  chrom_pos <-  vcf_file %>% 
    select(X1, X2) %>% 
    rownames_to_column(var = "order") %>% 
    rename(CHROM = X1, POS = X2) %>% 
    mutate(ID = str_c("snp", order),
           chrom_number = "1") %>% 
    select(ID, chrom_number, POS, order)
  
  write_delim(chrom_pos, paste0("snps_info_", chromosome, ".txt" ))


indiv_names <- read_delim(vcf_name, delim = '\t', skip = header_rows -1, n_max = 1, col_names = FALSE ) %>% 
  select(-num_range("X", 1:9)) %>% 
  unlist(., use.names=FALSE)


# Make genotype info files
vcf <- vcf_file %>% 
  select(-X1, -X2) %>% 
  rename_at(vars(num_range("X", 10:178)), ~ indiv_names)

recoded <- lapply(indiv_names, recode_genos, vcf = vcf) %>%
  reduce(full_join, by = "rowname") %>%
  rename_at(vars(ends_with("_num_rep")),
            funs(str_remove(.,"_num_rep")))

recoded_long <- recoded %>%
  gather(variable, value, -rowname)

recoded_wide <- recoded_long %>%
  mutate(rowname = as.numeric(rowname)) %>%
  spread(rowname, value) %>%
  unite(col = genotypes, -variable, sep = "") %>%
  mutate(indiv = str_remove(variable, "_num_rep")) %>%
  mutate(chip_num = "1") %>%
  select(indiv, chip_num, genotypes)

write_delim(recoded_wide, paste0("chr_", chromosome,".txt"), delim = ":", col_names = TRUE)


# Make pedigree file
pedigree_file <- read_delim(bird_data_file, delim = ",") %>% 
  # add_row(ID = "Konini_3-4-16", mother = "Konini", father = "Blake", hatch_year = "2016", sex = "M") %>% 
  select(ID, father, mother, sex) %>% 
  mutate(mother = replace_na(mother, "0"),
         father = replace_na(father, "0"))


write_delim(pedigree_file, "fimpute_pedigree.txt", delim = "\t", col_names = TRUE)


# Make config file 
config_file <-  paste(paste0('title="Phasing_', chromosome, '";'),
                        paste0('genotype_file="fimpute_chr_', chromosome,  '.txt";'),
                        paste0('snp_info_file="snps_info_', chromosome,  '.txt";'),
                        'ped_file="fimpute_pedigree.txt";',
                        paste0('output_folder="output_', chromosome, '";'),
                        'parentage_test /ert_mm = 0.05;',
                       # 'turnoff_fam;', # Toggle family imputation
                        'turnoff_pop;',
                        'turnoff_random_fill;',
                        'save_hap_lib /id;',
                        'njob=5;',
                        'ped_depth=4;', sep = "\n")
  
  
  write_lines(config_file, paste0("fimpute_config_", chromosome,  ".txt"))
  
  

  
    
  
  