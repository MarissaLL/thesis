library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(readr)
library(purrr)


#############
# FUNCTIONS #
#############
add_name <- function(name, dataset){
  dataset[name] %>% 
    reduce(unlist(.)) %>% 
    mutate(name = str_remove(name, ".txt"))
}

get_per_chr_info <- function(idxstats, autosomes, w_chr, z_chr){
  dataset <- idxstats %>% 
    dplyr::rename(seq_name = X1, seq_len = X2, num_mapped_read_segs = X3, num_unmapped_read_segs = X4) %>% 
    mutate(perc_tot = (num_mapped_read_segs/sum(num_mapped_read_segs))*100) %>% 
    mutate(length_norm = (num_mapped_read_segs/seq_len)*1000) %>% 
    filter(str_detect(seq_name, "S") == TRUE) %>% 
    mutate(category = case_when(seq_name %in% autosomes ~ "autosome",
                                seq_name %in% w_chr ~ "W",
                                seq_name %in% z_chr ~ "Z",
                                TRUE ~ "discard")) %>% 
    filter(category != "discard") 
 return(dataset)
}



test_chr_differences <- function(idxstats, autosomes, w_chr, z_chr){
  dataset <- get_per_chr_info(idxstats, autosomes, w_chr, z_chr)
  write_delim(dataset, paste0(idxstats, "summary.txt"), delim = "\t")
  aov_res <- aov(dataset$length_norm ~ dataset$category)
  
  tukey_res <- TukeyHSD(aov_res)
  mean_cov <- dataset %>% 
    filter(category == "autosome") %>% 
    summarise(mean = mean(length_norm)) %>% 
    pull()

  
  cc <- as_tibble(tukey_res$`dataset$category`) %>% 
    mutate(name =  unique(dataset$name)) %>% 
    mutate(comparison = rownames(tukey_res$`dataset$category`)) %>% 
    mutate(difference = case_when(comparison == "W-autosome" ~ -diff,
                                  comparison == "Z-autosome" ~ -diff,
                                  comparison == "Z-W" ~ -diff,
                                  TRUE ~ diff),
           comparison =  case_when(comparison == "W-autosome" ~ "autosome-W",
                                   comparison == "Z-autosome" ~ "autosome-Z",
                                   comparison == "Z-W" ~ "W-Z",
                                   TRUE ~ comparison)) %>% 
    select(-lwr, -upr, -diff) %>% 
    pivot_longer(cols = c(difference, `p adj`), names_to = "type", values_to = "value") %>% 
    pivot_wider(id_cols = c(name, type), names_from = comparison, values_from = value)
  
  differences <- cc %>% 
    filter(type == "difference") %>% 
    rename(WZ_diff = `W-Z`, autoZ_diff = `autosome-Z`, autoW_diff = `autosome-W`)
  
  p_vals <- cc %>% 
    filter(type == "p adj") %>% 
    rename(WZ_pval = `W-Z`, autoZ_pval = `autosome-Z`, autoW_pval = `autosome-W`) 
  
  res <- full_join(differences, p_vals, by = "name") %>% 
    select(-type.x, -type.y) %>% 
    mutate(mean_cov = mean_cov)
  
  return(res)
}

# Difference is <0 when the value for the second is larger
assign_sex <- function(dataset){
  dataset %>% 
    mutate(WZ_evidence = case_when(WZ_diff < 0 & WZ_pval < 0.05 ~ "M",
                                   WZ_diff > 0 ~ "F",
                                   TRUE ~ "uninformative"),
           autoZ_evidence = case_when(autoZ_diff < 0 & autoZ_pval < 0.05 ~ "M",
                                      autoZ_diff > 0 & autoZ_pval < 0.05 ~ "F",
                                      TRUE ~ "uninformative"),
           WZ_diff_only = case_when(WZ_diff < 0 ~ "M", 
                                      TRUE ~ "uninformative")) %>% 
    mutate(result = case_when(WZ_evidence == "M" & autoZ_evidence != "F" ~ "male",
                              WZ_evidence == "F" & autoZ_evidence != "M" ~ "female",
                              WZ_evidence != "F" & autoZ_evidence == "M" ~ "male",
                              WZ_evidence != "M" & autoZ_evidence == "F" ~ "female",
                              TRUE ~ "low_conf"),
           confidence = case_when(result == "male" | result == "female" ~ "high",
                                  result == "low_conf" ~ "low"),
           result_all = case_when(result != "low_conf" ~ result,
                                  result == "low_conf" & WZ_diff == "M" ~ "male",
                                  result == "low_conf" & WZ_diff != "M" ~ "male"))
}



###########
# GLOBALS #
###########


# idxstats_filepath <- "~/projects/kakapo-test/"
# idxstats_filepath <- "/media/drive_6tb/projects/kakapo-genomics/genome_files/old_assembly"
#idxstats_filepath <- "/media/drive_6tb/projects/kakapo-genomics/genome_files/new_assembly/"
# idxstats_filepath <- "/media/drive_6tb/projects/kakapo-genomics/output/17_processing_np_chicks/sex_test/"
idxstats_filepath <- "./"

# Kakapo
autosomes <- paste0(rep("S"),c(1,2,4:12,14:20))
w_chr <- "S13"
z_chr <- "S3"

# # Takahe
# idxstats_filepath <- "./"
 
# autosomes <- c("SUPER_7", "SUPER_8", "SUPER_9", "SUPER_12", "SUPER_13", "SUPER_14", 
#                "SUPER_15", "SUPER_16", "SUPER_17", "SUPER_18", "SUPER_2", "SUPER_19", 
#                "SUPER_20", "SUPER_21", "SUPER_22", "SUPER_23", "SUPER_24", "SUPER_25", 
#                "SUPER_26", "SUPER_27", "SUPER_1", "SUPER_28", "S29", "S30", 
#                "S31", "SUPER_3", "SUPER_30", "SUPER_29", "SUPER_4", "SUPER_31", 
#                "SUPER_10", "SUPER_11", "SUPER_5", "SUPER_32", "SUPER_6")
# w_chr <- c("SUPER_W_unloc_13", "SUPER_W_unloc_14", "SUPER_W", "SUPER_W_unloc_1", 
#            "SUPER_W_unloc_3", "SUPER_W_unloc_4", "SUPER_W_unloc_5", "SUPER_W_unloc_2", 
#            "SUPER_W_unloc_6", "SUPER_W_unloc_7", "SUPER_W_unloc_8", "SUPER_W_unloc_9", 
#            "SUPER_W_unloc_10", "SUPER_W_unloc_11", "SUPER_W_unloc_12")
# z_chr <- c("SUPER_Z_unloc_1", "SUPER_Z_unloc_2", "SUPER_Z_unloc_3", "SUPER_Z_unloc_4", 
#            "SUPER_Z_unloc_5", "SUPER_Z")



########
# MAIN #
########


all_idxstats <- list.files(path = paste0(idxstats_filepath), pattern = ".*\\.idxstats", full.names = TRUE) %>%
  setNames(., sub("\\.idxstats", "", basename(.))) %>% 
  map(., read_delim, delim = "\t", col_names = FALSE) 


tidy_idxstats <- lapply(names(all_idxstats), add_name, dataset = all_idxstats)

results_test <- lapply(tidy_idxstats, get_per_chr_info, autosomes = autosomes, w_chr = w_chr, z_chr = z_chr) %>%
  reduce(bind_rows) %>% write_delim("coverage_stats.txt", delim = "\t")
 

results_output <- lapply(tidy_idxstats, test_chr_differences, autosomes = autosomes, w_chr = w_chr, z_chr = z_chr) %>% 
  reduce(bind_rows)

results_output2 <- assign_sex(results_output) %>% 
  rename(sex = result_all) %>% 
  select(name, sex, confidence)
write_delim(results_output2, "assigned_sex.txt", delim = "\t")
