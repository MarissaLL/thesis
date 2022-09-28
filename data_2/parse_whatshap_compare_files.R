library(tidyverse)

# indivs <- read_delim(file = "names.txt", delim = ", ", col_names = FALSE) %>% 
#   str_replace("\\[", "") %>% 
#   str_replace("\\]", "")
#   
inputs <- commandArgs(trailingOnly = TRUE)

n_vcfs <- inputs[2]


indivs <- read_delim(inputs[1], delim = " ", col_names = FALSE) %>% 
  unlist(., use.names=FALSE) %>% 
  str_remove("whatshap_compare_") %>% 
  str_remove(".txt")

parse_wc_file_2_vcf <- function(indiv){
  wc_file <- read_lines(paste0("whatshap_compare_", indiv, ".txt")) %>% 
    trimws()
  file0 <- str_remove(wc_file[3], "file0 =\\s+")
  file1 <- str_remove(wc_file[4], "file1 =\\s+")
  variant_count_file0 <-  str_remove(wc_file[7], "file0:\\s+[:digit:]+\\s+/\\s+")
  het_count_file0 <- str_extract(wc_file[7], "(?<=file0:\\s)\\s+[:digit:]+")
  variant_count_file1 <-  str_remove(wc_file[8], "file1:\\s+[:digit:]+\\s+/\\s+")
  het_count_file1 <- str_extract(wc_file[8], "(?<=file1:\\s)\\s+[:digit:]+")
  union_all_vars <- str_remove(wc_file[9], "UNION:\\s+[:digit:]+\\s+/\\s+")
  union_het_vars <- str_extract(wc_file[9], "(?<=UNION:\\s)\\s+[:digit:]+")
  intersection_all_vars <- str_remove(wc_file[10], "INTERSECTION:\\s+[:digit:]+\\s+/\\s+")
  intersection_het_vars <- str_extract(wc_file[10], "(?<=INTERSECTION:\\s)\\s+[:digit:]+")
  hets_in_common <- str_remove(wc_file[12], "common heterozygous variants:\\s+")
  non_singleton_file0 <- str_remove(wc_file[14], "non-singleton blocks in file0:\\s+")
  blocks_cover_file0 <- str_remove(wc_file[15], "--> covered variants:\\s+")
  non_singleton_file1 <- str_remove(wc_file[16], "non-singleton blocks in file1:\\s+")
  blocks_cover_file1 <- str_remove(wc_file[17], "--> covered variants:\\s+")
  n_intersection_blocks <- str_remove(wc_file[18], "non-singleton intersection blocks:\\s+")
  intersection_blocks_cover <- str_remove(wc_file[19],  "--> covered variants:\\s+")
  phased_var_pairs <- str_remove(wc_file[21], "phased pairs of variants assessed:\\s+")
  switch_errors <- str_remove(wc_file[22], "switch errors:\\s+")
  switch_error_rate <- str_extract(wc_file[23], "(?<=switch error rate:)\\s+[:digit:]+\\.[:digit:]+(?=%)")
  switch_flip_decomposition <- str_remove(wc_file[24], "switch/flip decomposition:\\s+")
  switch_flip_rate <- str_extract(wc_file[25],  "(?<=switch/flip rate:)\\s+[:graph:]+(?=%)")
  blockwise_hamming <- str_remove(wc_file[26], "Block-wise Hamming distance:\\s+")
  blockwise_hamming_perc <- str_extract(wc_file[27], "(?<=Block-wise Hamming distance \\[%\\]:)\\s+[:graph:]+(?=%)")
  largest_intersection_block <- str_remove(wc_file[29], "phased pairs of variants assessed:\\s+")
  largest_switch_errors <- str_remove(wc_file[30], "switch errors:\\s+")
  largest_switch_rate <- str_extract(wc_file[31], "(?<=switch error rate:)\\s+[:graph:]+(?=%)")
  largest_switch_flip_decomposition <- str_remove(wc_file[32], "switch/flip decomposition:\\s+")
  largest_switch_flip_rate <- str_extract(wc_file[33], "(?<=switch/flip rate:)\\s+[:graph:]+(?=%)")
  largest_hamming_distance <- str_remove(wc_file[34], "Hamming distance:\\s+")
  largest_hamming_distance_perc <- str_extract(wc_file[35], "(?<=Hamming distance \\[%\\]:)\\s+[:graph:]+(?=%)")
  
  stats_tibble <- tibble(indiv, file0,file1, variant_count_file0, het_count_file0, variant_count_file1, het_count_file1, union_all_vars, union_het_vars, 
         intersection_all_vars, intersection_het_vars, hets_in_common, non_singleton_file0, blocks_cover_file0, non_singleton_file1, 
         blocks_cover_file1, n_intersection_blocks, intersection_blocks_cover, phased_var_pairs, switch_errors, switch_error_rate, switch_flip_decomposition, 
         switch_flip_rate, blockwise_hamming, blockwise_hamming_perc, largest_intersection_block, largest_switch_errors, largest_switch_rate, 
         largest_switch_flip_decomposition, largest_switch_flip_rate, largest_hamming_distance, largest_hamming_distance_perc)
  

  return(stats_tibble)
}

### This doesn't include all the fields the 2 vcf one does because I haven't used them. Will add if needed.
parse_wc_file_3_vcf <- function(indiv){
  wc_file <- read_lines(paste0("whatshap_compare_", indiv, ".txt")) %>% 
    trimws()
  # General
  file0 <- str_remove(wc_file[3], "file0 =\\s+")
  file1 <- str_remove(wc_file[4], "file1 =\\s+")
  file2 <- str_remove(wc_file[5], "file2 =\\s+")
  variant_count_file0 <-  str_remove(wc_file[8], "file0:\\s+[:digit:]+\\s+/\\s+")
  het_count_file0 <- str_extract(wc_file[8], "(?<=file0:\\s)\\s+[:digit:]+")
  variant_count_file1 <-  str_remove(wc_file[9], "file1:\\s+[:digit:]+\\s+/\\s+")
  het_count_file1 <- str_extract(wc_file[9], "(?<=file1:\\s)\\s+[:digit:]+")
  variant_count_file2 <-  str_remove(wc_file[10], "file2:\\s+[:digit:]+\\s+/\\s+")
  het_count_file2 <- str_extract(wc_file[10], "(?<=file2:\\s)\\s+[:digit:]+")
  union_all_vars <- str_remove(wc_file[11], "UNION:\\s+[:digit:]+\\s+/\\s+")
  union_het_vars <- str_extract(wc_file[11], "(?<=UNION:\\s)\\s+[:digit:]+")
  intersection_all_vars <- str_remove(wc_file[12], "INTERSECTION:\\s+[:digit:]+\\s+/\\s+")
  intersection_het_vars <- str_extract(wc_file[12], "(?<=INTERSECTION:\\s)\\s+[:digit:]+")
  # File0 and file1 comparison
  hets_in_common_01 <- str_remove(wc_file[14], "common heterozygous variants:\\s+")
  blocks_cover_file0_01 <- str_remove(wc_file[17], "--> covered variants:\\s+")
  blocks_cover_file1_01 <- str_remove(wc_file[19], "--> covered variants:\\s+")
  intersection_blocks_cover_01 <- str_remove(wc_file[21],  "--> covered variants:\\s+")
  phased_var_pairs_01 <- str_remove(wc_file[23], "phased pairs of variants assessed:\\s+")
  switch_errors_01 <- str_remove(wc_file[24], "switch errors:\\s+")
  switch_error_rate_01 <- str_extract(wc_file[25], "(?<=switch error rate:)\\s+[:digit:]+\\.[:digit:]+(?=%)")
  blockwise_hamming_01 <- str_remove(wc_file[28], "Block-wise Hamming distance:\\s+")
  # File0 and file2 comparison
  hets_in_common_02 <- str_remove(wc_file[39], "common heterozygous variants:\\s+")
  blocks_cover_file0_02 <- str_remove(wc_file[42], "--> covered variants:\\s+")
  blocks_cover_file2_02 <- str_remove(wc_file[44], "--> covered variants:\\s+")
  intersection_blocks_cover_02 <- str_remove(wc_file[46],  "--> covered variants:\\s+")
  phased_var_pairs_02 <- str_remove(wc_file[48], "phased pairs of variants assessed:\\s+")
  switch_errors_02 <- str_remove(wc_file[49], "switch errors:\\s+")
  switch_error_rate_02 <- str_extract(wc_file[50], "(?<=switch error rate:)\\s+[:digit:]+\\.[:digit:]+(?=%)")
  blockwise_hamming_02 <- str_remove(wc_file[53], "Block-wise Hamming distance:\\s+")
  # File1 and file2 comparison
  hets_in_common_12 <- str_remove(wc_file[64], "common heterozygous variants:\\s+")
  blocks_cover_file1_12 <- str_remove(wc_file[67], "--> covered variants:\\s+")
  blocks_cover_file2_12 <- str_remove(wc_file[69], "--> covered variants:\\s+")
  intersection_blocks_cover_12 <- str_remove(wc_file[71],  "--> covered variants:\\s+")
  phased_var_pairs_12 <- str_remove(wc_file[73], "phased pairs of variants assessed:\\s+")
  switch_errors_12 <- str_remove(wc_file[74], "switch errors:\\s+")
  switch_error_rate_12 <- str_extract(wc_file[75], "(?<=switch error rate:)\\s+[:digit:]+\\.[:digit:]+(?=%)")
  blockwise_hamming_12 <- str_remove(wc_file[78], "Block-wise Hamming distance:\\s+")
  #Overall agree/disagree
  covered_3_intersection <- str_remove(wc_file[98], "--> covered variants:\\s+")
  compared_pairs_3 <- str_remove(wc_file[99], "Compared pairs of variants:\\s+")
  n_agree_all <- str_extract(wc_file[101], "(?<=\\{file0,file1,file2\\} vs. \\{\\}:)\\s+[:digit:]+")
  perc_agree_all <- str_extract(wc_file[101], "[:digit:]+\\.[:digit:]+(?=%$)")
  file01v2 <-  str_extract(wc_file[103], "(?<=\\{file0,file1\\} vs. \\{file2\\}:)\\s+[:digit:]+")
  perc_file01v2  <- str_extract(wc_file[103], "[:digit:]+\\.[:digit:]+(?=%$)")
  file0v12 <-  str_extract(wc_file[104], "(?<=\\{file0\\} vs. \\{file1,file2\\}:)\\s+[:digit:]+")
  perc_file0v12  <- str_extract(wc_file[104], "[:digit:]+\\.[:digit:]+(?=%$)")
  
  stats_tibble <- tibble(indiv, file0,file1, file2, variant_count_file0, het_count_file0, variant_count_file1, het_count_file1, 
                         variant_count_file2, het_count_file2, union_all_vars, union_het_vars, intersection_all_vars, intersection_het_vars, 
                         # File0 and file1 comparison
                         hets_in_common_01, blocks_cover_file0_01, blocks_cover_file1_01, intersection_blocks_cover_01, phased_var_pairs_01, switch_errors_01,
                         switch_error_rate_01 , blockwise_hamming_01,
                         # File0 and file2 comparison
                         hets_in_common_02, blocks_cover_file0_02, blocks_cover_file2_02, intersection_blocks_cover_02, phased_var_pairs_02, switch_errors_02,
                         switch_error_rate_02 , blockwise_hamming_02,
                         # File1 and file2 comparison
                         hets_in_common_12, blocks_cover_file1_12, blocks_cover_file2_12, intersection_blocks_cover_12, phased_var_pairs_12, switch_errors_12,
                         switch_error_rate_12 , blockwise_hamming_12,
                         #Overall agree/disagree
                         covered_3_intersection, compared_pairs_3, n_agree_all, perc_agree_all,
                         file01v2, perc_file01v2, file0v12, perc_file0v12
                         )
  
  
  return(stats_tibble)
}



if(n_vcfs == 2){
  all_so_far <- lapply(indivs, parse_wc_file_2_vcf) %>% 
    reduce(bind_rows) %>% 
    mutate(across(.cols = c(variant_count_file0, het_count_file0, variant_count_file1, het_count_file1, union_all_vars, union_het_vars, 
                            intersection_all_vars, intersection_het_vars, hets_in_common, non_singleton_file0, blocks_cover_file0, non_singleton_file1, 
                            blocks_cover_file1, n_intersection_blocks, intersection_blocks_cover, phased_var_pairs, switch_errors, switch_error_rate,  
                            switch_flip_rate, blockwise_hamming, blockwise_hamming_perc, largest_intersection_block, largest_switch_errors, largest_switch_rate, 
                            largest_switch_flip_rate, largest_hamming_distance, largest_hamming_distance_perc), 
                  .fns = as.numeric)) %>% 
    mutate(non_intersection_blocks_file0 = blocks_cover_file0 - intersection_blocks_cover,
           non_intersection_blocks_file1 = blocks_cover_file1 - intersection_blocks_cover)
}

if(n_vcfs == 3){
  all_so_far <- lapply(indivs, parse_wc_file_3_vcf) %>% 
    reduce(bind_rows) %>% 
    mutate(across(.cols = c(-indiv, -file0, -file1, -file2), 
                  .fns = as.numeric))
}





write_delim(all_so_far, "all_whatshap_compare.txt", delim = "\t")

