library(tidyverse)
library(VennDiagram)

#############
# FUNCTIONS #
#############

get_data <- function(file){
  read_delim(file, delim = "\t") %>% 
    mutate(file0 = str_extract(file_name0, "^[:alpha:]+(?=_)"),
           file1 = str_extract(file_name1, "^[:alpha:]+(?=_)"),
           both = paste0(file0, "_", file1)) %>% 
    select(-c(dataset_name0, dataset_name1, file_name0, file_name1, file0, file1))
}

parse_wc_lines <- function(wc_file){
  chromosome <- str_remove(wc_file[1], "-+$")
  variant_count_file0 <-  str_remove(wc_file[3], "file0:\\s+[:digit:]+\\s+/\\s+")
  het_count_file0 <- str_extract(wc_file[3], "(?<=file0:\\s)\\s+[:digit:]+")
  variant_count_file1 <-  str_remove(wc_file[4], "file1:\\s+[:digit:]+\\s+/\\s+")
  het_count_file1 <- str_extract(wc_file[4], "(?<=file1:\\s)\\s+[:digit:]+")
  union_all_vars <- str_remove(wc_file[5], "UNION:\\s+[:digit:]+\\s+/\\s+")
  union_het_vars <- str_extract(wc_file[5], "(?<=UNION:\\s)\\s+[:digit:]+")
  intersection_all_vars <- str_remove(wc_file[6], "INTERSECTION:\\s+[:digit:]+\\s+/\\s+")
  intersection_het_vars <- str_extract(wc_file[6], "(?<=INTERSECTION:\\s)\\s+[:digit:]+")
  hets_in_common <- str_remove(wc_file[8], "common heterozygous variants:\\s+")
  non_singleton_file0 <- str_remove(wc_file[10], "non-singleton blocks in file0:\\s+")
  blocks_cover_file0 <- str_remove(wc_file[11], "--> covered variants:\\s+")
  non_singleton_file1 <- str_remove(wc_file[12], "non-singleton blocks in file1:\\s+")
  blocks_cover_file1 <- str_remove(wc_file[13], "--> covered variants:\\s+")
  n_intersection_blocks <- str_remove(wc_file[14], "non-singleton intersection blocks:\\s+")
  intersection_blocks_cover <- str_remove(wc_file[15],  "--> covered variants:\\s+")
  phased_var_pairs <- str_remove(wc_file[17], "phased pairs of variants assessed:\\s+")
  switch_errors <- str_remove(wc_file[18], "switch errors:\\s+")
  switch_error_rate <- str_extract(wc_file[19], "(?<=switch error rate:)\\s+[:digit:]+\\.[:digit:]+(?=%)")
  switch_flip_decomposition <- str_remove(wc_file[20], "switch/flip decomposition:\\s+")
  switch_flip_rate <- str_extract(wc_file[21],  "(?<=switch/flip rate:)\\s+[:graph:]+(?=%)")
  blockwise_hamming <- str_remove(wc_file[22], "Block-wise Hamming distance:\\s+")
  blockwise_hamming_perc <- str_extract(wc_file[23], "(?<=Block-wise Hamming distance \\[%\\]:)\\s+[:graph:]+(?=%)")
  largest_intersection_block <- str_remove(wc_file[25], "phased pairs of variants assessed:\\s+")
  largest_switch_errors <- str_remove(wc_file[26], "switch errors:\\s+")
  largest_switch_rate <- str_extract(wc_file[27], "(?<=switch error rate:)\\s+[:graph:]+(?=%)")
  largest_switch_flip_decomposition <- str_remove(wc_file[28], "switch/flip decomposition:\\s+")
  largest_switch_flip_rate <- str_extract(wc_file[29], "(?<=switch/flip rate:)\\s+[:graph:]+(?=%)")
  largest_hamming_distance <- str_remove(wc_file[30], "Hamming distance:\\s+")
  largest_hamming_distance_perc <- str_extract(wc_file[31], "(?<=Hamming distance \\[%\\]:)\\s+[:graph:]+(?=%)")
  
  stats_tibble <- tibble(chromosome, variant_count_file0, het_count_file0, variant_count_file1, het_count_file1, union_all_vars, union_het_vars, 
                         intersection_all_vars, intersection_het_vars, hets_in_common, non_singleton_file0, blocks_cover_file0, non_singleton_file1, 
                         blocks_cover_file1, n_intersection_blocks, intersection_blocks_cover, phased_var_pairs, switch_errors, switch_error_rate, switch_flip_decomposition, 
                         switch_flip_rate, blockwise_hamming, blockwise_hamming_perc, largest_intersection_block, largest_switch_errors, largest_switch_rate, 
                         largest_switch_flip_decomposition, largest_switch_flip_rate, largest_hamming_distance, largest_hamming_distance_perc)
  
  
  return(stats_tibble)
}

parse_wc_lines2 <- function(wc_file, n_compared){
  n_pair_comparisons <- factorial(n_compared)/(2*factorial(n_compared-2))
  chromosome <- str_remove(wc_file[1], "-+$")
  variant_counts <-  lapply(c(3:(2+n_compared)), function(x){str_remove(wc_file[x], "file[:digit:]:\\s+[:digit:]+\\s+/\\s+")}) %>% unlist()
  het_counts <- lapply(c(3:(2+n_compared)), function(x){str_extract(wc_file[x], "(?<=file[:digit:]:\\s)\\s+[:digit:]+")}) %>% unlist()
  union_all_vars <- str_remove(wc_file[3+n_compared], "UNION:\\s+[:digit:]+\\s+/\\s+")
  union_het_vars <- str_extract(wc_file[3+n_compared], "(?<=UNION:\\s)\\s+[:digit:]+")
  intersection_all_vars <- str_remove(wc_file[4+n_compared], "INTERSECTION:\\s+[:digit:]+\\s+/\\s+")
  intersection_het_vars <- str_extract(wc_file[4+n_compared], "(?<=INTERSECTION:\\s)\\s+[:digit:]+")
  if (n_compared > 2){
    hets_in_common <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +6 + n_compared), 
                             function(x){str_remove(wc_file[x], "common heterozygous variants:\\s+")}) %>% unlist()
    non_singleton_file0 <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +8 + n_compared), 
                                  function(x){str_remove(wc_file[x], "non-singleton blocks in file[:digit:]:\\s+")}) %>% unlist()
    blocks_cover_file0 <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +9 + n_compared), 
                                 function(x){str_remove(wc_file[x], "--> covered variants:\\s+")}) %>% unlist()
    non_singleton_file1 <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +10 + n_compared), 
                                  function(x){str_remove(wc_file[x], "non-singleton blocks in file[:digit:]:\\s+")}) %>% unlist()
    blocks_cover_file1 <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +11 + n_compared), 
                                 function(x){str_remove(wc_file[x], "--> covered variants:\\s+")}) %>% unlist()
    n_intersection_blocks <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +12 + n_compared), 
                                    function(x){str_remove(wc_file[x], "non-singleton intersection blocks:\\s+")}) %>% unlist()
    intersection_blocks_cover <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +13 + n_compared),
                                        function(x){str_remove(wc_file[x],  "--> covered variants:\\s+")}) %>% unlist()
    phased_var_pairs <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +15 + n_compared),
                               function(x){str_remove(wc_file[x], "phased pairs of variants assessed:\\s+")}) %>% unlist()
    switch_errors <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +16 + n_compared),
                            function(x){str_remove(wc_file[x], "switch errors:\\s+")}) %>% unlist()
    switch_error_rate <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +17 + n_compared),
                                function(x){str_extract(wc_file[x], "(?<=switch error rate:)\\s+[:digit:]+\\.[:digit:]+(?=%)")}) %>% unlist()
    switch_flip_decomposition <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +18 + n_compared),
                                        function(x){str_remove(wc_file[x], "switch/flip decomposition:\\s+")}) %>% unlist()
    switch_flip_rate <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +19 + n_compared),
                               function(x){str_extract(wc_file[x],  "(?<=switch/flip rate:)\\s+[:graph:]+(?=%)")}) %>% unlist()
    blockwise_hamming <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +20 + n_compared),
                                function(x){str_remove(wc_file[x], "Block-wise Hamming distance:\\s+")}) %>% unlist()
    blockwise_hamming_perc <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +21 + n_compared),
                                     function(x){str_extract(wc_file[x], "(?<=Block-wise Hamming distance \\[%\\]:)\\s+[:graph:]+(?=%)")}) %>% unlist()
    largest_intersection_block <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +23 + n_compared),
                                         function(x){str_remove(wc_file[x], "phased pairs of variants assessed:\\s+")}) %>% unlist()
    largest_switch_errors <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +24 + n_compared),
                                    function(x){str_remove(wc_file[x], "switch errors:\\s+")}) %>% unlist()
    largest_switch_rate <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +25 + n_compared),
                                  function(x){str_extract(wc_file[x], "(?<=switch error rate:)\\s+[:graph:]+(?=%)")}) %>% unlist()
    largest_switch_flip_decomposition <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +26 + n_compared),
                                                function(x){str_remove(wc_file[x], "switch/flip decomposition:\\s+")}) %>% unlist()
    largest_switch_flip_rate <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +27 + n_compared),
                                       function(x){str_extract(wc_file[x], "(?<=switch/flip rate:)\\s+[:graph:]+(?=%)")}) %>% unlist()
    largest_hamming_distance <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +28 + n_compared),
                                       function(x){str_remove(wc_file[x], "Hamming distance:\\s+")}) %>% unlist()
    largest_hamming_distance_perc <- lapply(c((seq(from = 0, to = n_pair_comparisons-1)*25) +29 + n_compared),
                                            function(x){str_extract(wc_file[x], "(?<=Hamming distance \\[%\\]:)\\s+[:graph:]+(?=%)")}) %>% unlist()
    
    #Overall agree/disagree
    multiway_all <- str_remove(wc_file[(n_pair_comparisons-1)*25 + 31 + n_compared], "common heterozygous variants:\\s+")
    multiway_intersection <- str_remove(wc_file[(n_pair_comparisons-1)*25 + 34 + (n_compared*3)], "--> covered variants:\\s+")
    compared_pairs_multiway <- str_remove(wc_file[(n_pair_comparisons-1)*25 + 35 + (n_compared*3)], "Compared pairs of variants:\\s+")
    n_agree_all <- str_extract(wc_file[(n_pair_comparisons-1)*25 + 37 + (n_compared*3)], "(?<=:)\\s+[:digit:]+")
    perc_agree_all <- str_extract(wc_file[(n_pair_comparisons-1)*25 + 37 + (n_compared*3)], "[:digit:]+\\.[:digit:]+(?=%$)")
    disagreement_names <- lapply(c(seq(from = (((n_pair_comparisons-1)*25) +39 + (n_compared*3)), to =length(wc_file)-1)),
                                 function(x){str_extract(wc_file[x], "[:graph:]+\\s+vs.\\s+[:graph:]*(?=:)")}) %>% unlist()
    disagreement_n <- lapply(c(seq(from = (((n_pair_comparisons-1)*25) +39 + (n_compared*3)), to =length(wc_file)-1)),
                             function(x){str_extract(wc_file[x], "(?<=:)\\s+[:digit:]+(?=\\s+)")}) %>% unlist()
    disagreement_perc <- lapply(c(seq(from = (((n_pair_comparisons-1)*25) +39 + (n_compared*3)), to =length(wc_file)-1)),
                                function(x){str_extract(wc_file[x], "(?<=\\s)[:graph:]+$")}) %>% unlist()
    
  }
  
  tests <- tibble(names = c("chromosome", paste0("variant_counts_",seq(1:n_compared)), paste0("het_counts_", seq(1:n_compared)), 
                   "union_all_vars", "union_het_vars", "intersection_all_vars", "intersection_het_vars", 
                   paste0("hets_in_common_",seq(1:n_pair_comparisons)), paste0("non_singleton_file0_",seq(1:n_pair_comparisons)),  
                   paste0("blocks_cover_file0_",seq(1:n_pair_comparisons)), paste0("non_singleton_file1_",seq(1:n_pair_comparisons)), 
                   paste0("blocks_cover_file1_",seq(1:n_pair_comparisons)), paste0("n_intersection_blocks_",seq(1:n_pair_comparisons)), 
                   paste0("intersection_blocks_cover_",seq(1:n_pair_comparisons)), paste0("phased_var_pairs_",seq(1:n_pair_comparisons)),  
                   paste0("switch_errors_", seq(1:n_pair_comparisons)), paste0("switch_error_rate_",seq(1:n_pair_comparisons)), 
                   paste0("switch_flip_decomposition_",seq(1:n_pair_comparisons)), paste0("switch_flip_rate_",seq(1:n_pair_comparisons)), 
                   paste0("blockwise_hamming_",seq(1:n_pair_comparisons)), paste0("blockwise_hamming_perc_",seq(1:n_pair_comparisons)), 
                   paste0("largest_intersection_block_",seq(1:n_pair_comparisons)), paste0("largest_switch_errors_",seq(1:n_pair_comparisons)), 
                   paste0("largest_switch_rate_",seq(1:n_pair_comparisons)),  paste0("largest_switch_flip_decomposition_",seq(1:n_pair_comparisons)), 
                   paste0("largest_switch_flip_rate_",seq(1:n_pair_comparisons)), paste0("largest_hamming_distance_",seq(1:n_pair_comparisons)), 
                   paste0("largest_hamming_distance_perc_", seq(1:n_pair_comparisons)) , "multiway_all", "multiway_intersection", "compared_pairs_multiway",
                   "n_agree_all", "perc_agree_all", paste0("disagreement_names_", seq(1:length(disagreement_names))), 
                   paste0("disagreement_n_", seq(1:length(disagreement_n))), paste0("disagreement_perc_", seq(1:length(disagreement_perc)))),
    values = c(chromosome, variant_counts, het_counts, union_all_vars, union_het_vars, intersection_all_vars,intersection_het_vars, hets_in_common, non_singleton_file0, blocks_cover_file0, non_singleton_file1, blocks_cover_file1, n_intersection_blocks, intersection_blocks_cover, phased_var_pairs, switch_errors, switch_error_rate, switch_flip_decomposition,switch_flip_rate, blockwise_hamming, blockwise_hamming_perc, largest_intersection_block, largest_switch_errors, largest_switch_rate,  largest_switch_flip_decomposition, largest_switch_flip_rate, largest_hamming_distance, largest_hamming_distance_perc , multiway_all, multiway_intersection, compared_pairs_multiway, n_agree_all, perc_agree_all, disagreement_names, disagreement_n, disagreement_perc)) %>% 
    pivot_wider(names_from = names, values_from = values)
    return(tests)
}

read_pairwise_log <- function(indiv, pair){
  wc_file <- read_file(paste0(indiv, "_", pair, ".log"))
  
  test <- str_split(wc_file, "Chromosome") %>% 
    unlist()
  
  collated <- lapply(2:length(test), function(chr_chunk){
    str_split(test[chr_chunk], "\n", simplify = TRUE) %>% 
      trimws() %>% 
      parse_wc_lines()}) %>% 
    reduce(bind_rows) %>% 
    mutate(pair = pair,
           indiv = indiv,
           switch_only = str_extract(switch_flip_decomposition, "[:digit:]+(?=/)"),
           flip_only = str_extract(switch_flip_decomposition, "(?<=/)[:digit:]+$"),
           across(c("variant_count_file0", "het_count_file0", "variant_count_file1", 
             "het_count_file1", "union_all_vars", "union_het_vars", "intersection_all_vars", 
             "intersection_het_vars", "hets_in_common", "non_singleton_file0", 
             "blocks_cover_file0", "non_singleton_file1", "blocks_cover_file1", 
             "n_intersection_blocks", "intersection_blocks_cover", "phased_var_pairs", 
             "switch_errors", "switch_error_rate","switch_flip_rate", "switch_only", "flip_only", "blockwise_hamming", "blockwise_hamming_perc", 
             "largest_intersection_block", "largest_switch_errors", "largest_switch_rate", "largest_switch_flip_rate", 
             "largest_hamming_distance", "largest_hamming_distance_perc"), as.numeric)) %>% 
    select(-switch_flip_decomposition)
  return(collated)
}


###########
# GLOBALS #
###########
res_dir <- "/media/drive_6tb/projects/kakapo-genomics/output/29_whatshap_compare/Aparima/"
# res_dir <- "~/Downloads/whatshap_compares/"
########
# MAIN #
########

#############################################
# Pairwise comparisons between phasing sets #
#############################################

pairs <- c("beagle_whatshap", "alphapeel_fimpute", "alphapeel_whatshap", "beagle_fimpute", "alphapeel_beagle", "fimpute_whatshap")

indiv <- read_lines("~/projects/kakapo-genomics/data/all_seq_birdnames.txt")



setwd("~/Downloads/whatshap_compares/")
all <- lapply(pairs, function(pair){lapply(indiv, read_pairwise_log, pair = pair) %>% 
  reduce(bind_rows)}) %>% 
  reduce(bind_rows) %>% 
  mutate(chromosome = str_remove_all(chromosome, " ")) %>% 
  mutate(chromosome = factor(chromosome, levels = paste0("S",c(1,2,4:12,14:26))),
         switch_only_rate = (switch_only/phased_var_pairs)*100,
         flip_only_rate = (flip_only/phased_var_pairs)*100) %>% 
  rename(combined_switch_flip_rate = switch_error_rate)


ggplot(all, aes(x = pair, y = phased_var_pairs, colour = chromosome))+ 
  geom_point(position = position_jitter(width = 0.2))+
  labs(y = "Shared phased variant pairs", x = "Pairwise comparison", colour = "Chromosome") +
  theme_classic()+ 
  theme(legend.position = "right",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA,
                                        colour = "black",
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.6)))+
  scale_x_discrete(labels = function(labels) {
    fixedLabels <- c()
    for (l in 1:length(labels)) {
      fixedLabels[l] <- paste0(ifelse(l %% 2 == 0, '', '\n'), labels[l])
    }
    return(fixedLabels)
  })


ggplot(subset(all, phased_var_pairs > 100), aes(x = pair, y = switch_only_rate))+
  geom_point(position = position_jitter(width = 0.2))+
  labs(x = "Pairwise comparison", y = "Switch errors (% of phased pairs)" )+
  ylim(-2,50)+
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA,
                                        colour = "black",
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.6)))+
  scale_x_discrete(labels = function(labels) {
    fixedLabels <- c()
    for (l in 1:length(labels)) {
      fixedLabels[l] <- paste0(ifelse(l %% 2 == 0, '', '\n'), labels[l])
    }
    return(fixedLabels)
  })


ggplot(subset(all, phased_var_pairs > 100), aes(x = pair, y = flip_only_rate))+
  geom_point(position = position_jitter(width = 0.2))+
  labs(x = "Pairwise comparison", y = "Flip error rate (% of phased variants)")+
  ylim(-2,50)+
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA,
                                        colour = "black",
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.6)))+
  scale_x_discrete(labels = function(labels) {
    fixedLabels <- c()
    for (l in 1:length(labels)) {
      fixedLabels[l] <- paste0(ifelse(l %% 2 == 0, '', '\n'), labels[l])
    }
    return(fixedLabels)
  })


ggplot(subset(all, phased_var_pairs >100), aes(x = pair, y = blockwise_hamming_perc))+ 
  geom_point(position = position_jitter(width = 0.2))+
  labs(x = "Pairwise comparison", y = "Blockwise percentage dissimilarity (Hamming distance)")+
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA,
                                        colour = "black",
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.6)))+
  scale_x_discrete(labels = function(labels) {
    fixedLabels <- c()
    for (l in 1:length(labels)) {
      fixedLabels[l] <- paste0(ifelse(l %% 2 == 0, '', '\n'), labels[l])
    }
    return(fixedLabels)
  })

all %>% 
  group_by(pair) %>% 
  summarise(phased_var_pairs = mean(phased_var_pairs), switch_only_rate = mean(switch_only_rate, na.rm = TRUE), 
            flip_only_rate = mean(flip_only_rate, na.rm = TRUE), blockwise_hamming_perc = mean(blockwise_hamming_perc, na.rm = TRUE))


# Look into this some more - any specific founders? others? that have bad phasing
ggplot(subset(all, phased_var_pairs >100), aes(x = indiv, y = switch_only_rate, colour = pair))+ 
  geom_point(position = position_jitter(width = 0.2))

###########################
# Compared all approaches #
###########################

res_dir <- "/cifs/biocldap/archive/deardenlab/marissalelec/projects/kakapo-genomics/output/whatshap_compare/all/"

listed_files <- list.files(path = res_dir, 
             pattern = "whatshap_compare_all_", full.names = TRUE) 

all_compared <- listed_files %>%
             map(function(x){
               indiv <- basename(x) %>% 
                 str_remove("whatshap_compare_all_") %>% 
                 str_remove(".log")
               message(indiv)
               test <- read_file(x) %>% 
                 str_split("Chromosome") %>% 
                 unlist()
               collated <- lapply(2:length(test), function(chr_chunk){
                 str_split(test[chr_chunk], "\n", simplify = TRUE) %>% 
                   trimws() %>% 
                   parse_wc_lines2(n_compared = 4) %>% 
                   mutate(ID = indiv)
               })}) %>%
             reduce(bind_rows) %>% 
  mutate(chromosome = factor(trimws(chromosome), levels = c(paste0("S",c(1,2,4:12,14:26)))))



ggplot(subset(all_compared, is.na(perc_agree_all)==FALSE), aes(x = ID, 
                                                               y = as.numeric(perc_agree_all), 
                                                               colour = chromosome))+
  geom_point()+
  scale_colour_viridis_d() +
  labs(y = "Phased variant pairs that agree in all datasets (%)", x = "Individual", colour = "Chromosome") +
  theme_classic()+ 
  geom_hline(yintercept = 98) +
  theme(legend.position = "none",
        legend.text = element_text(size = rel(1.6)),
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA,
                                        colour = "black",
                                        size = 0.7 ),
        axis.text.y = element_text(size = rel(1.6)),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = rel(1)))



# Those that didn't have high agreement in all four phased datasets
low_agreement <- all_compared %>% 
  filter(as.numeric(perc_agree_all) < 98 ) %>% 
  pivot_longer(cols = matches("disagreement_names_")) %>% 
  mutate(actual_perc = case_when(name == "disagreement_names_1" ~ disagreement_perc_1,
                                 name == "disagreement_names_2" ~ disagreement_perc_2,
                                 name == "disagreement_names_3" ~ disagreement_perc_3,
                                 name == "disagreement_names_4" ~ disagreement_perc_4,
                                 name == "disagreement_names_5" ~ disagreement_perc_5,
                                 name == "disagreement_names_6" ~ disagreement_perc_6,
                                 name == "disagreement_names_7" ~ disagreement_perc_7),
         actual_n = case_when(name == "disagreement_names_1" ~ disagreement_n_1,
                              name == "disagreement_names_2" ~ disagreement_n_2,
                              name == "disagreement_names_3" ~ disagreement_n_3,
                              name == "disagreement_names_4" ~ disagreement_n_4,
                              name == "disagreement_names_5" ~ disagreement_n_5,
                              name == "disagreement_names_6" ~ disagreement_n_6,
                              name == "disagreement_names_7" ~ disagreement_n_7)) %>% 
  select(chromosome, ID, perc_agree_all, value, actual_perc, actual_n) %>% 
  filter(is.na(actual_perc) == FALSE) %>% 
  mutate(translation  = case_when(value == "{file0,file1,file2} vs. {file3}" ~ "whatshap_disagrees",
                                  value == "{file0,file1,file3} vs. {file2}" ~ "fimpute_disagrees",
                                  value == "{file0,file2,file3} vs. {file1}" ~ "alphapeel_disagrees",
                                  value == "{file0} vs. {file1,file2,file3}" ~ "beagle_disagrees",
                                  TRUE ~ "even_split")) %>% 
  select(-value, -actual_n) %>% 
  mutate(perc_disagreement = 100 - as.numeric(perc_agree_all) ) %>% 
  mutate(actual_perc = as.numeric(str_remove(actual_perc, "%"))) %>% 
  mutate(contribution = (as.numeric(actual_perc)/perc_disagreement)*100) %>% 
  mutate(group = str_c(ID,chromosome, sep = "_")) %>% 
  mutate(translation = factor(translation, levels = c("alphapeel_disagrees", "beagle_disagrees", "fimpute_disagrees", "whatshap_disagrees", "even_split"), 
                              labels = c("Alphapeel", "Beagle", "FImpute3", "WhatsHap", "Even split")))



for_pplt <- low_agreement %>%
  select(group, translation, contribution) %>% 
  group_by(group, translation) %>% 
  summarise(contribution = sum(contribution)) %>% 
  pivot_wider(names_from = translation, values_from = contribution) %>% 
  arrange(desc(Alphapeel)) %>% 
  ungroup() %>% 
  mutate(order = seq(1,max(row_number()))) %>% 
  select(group, order)

actual_pplt <- low_agreement %>% 
  left_join(for_pplt) 

ggplot(actual_pplt, aes(x = reorder(group, 
                                      order), y = contribution, fill = translation))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_viridis_d() +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c(0,25,50,75,100))+
  labs(y = "Percentage of disagreeing phased pairs", x = "Chromosome x Individual", fill = "Dataset which disagrees") +
  theme_classic()+ 
  theme(legend.position = "right",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = rel(1.6)),
        axis.text.x = element_blank())

## All agreement - not filtered for low
agreement <- all_compared %>% 
  pivot_longer(cols = matches("disagreement_names_")) %>% 
  mutate(actual_perc = case_when(name == "disagreement_names_1" ~ disagreement_perc_1,
                                 name == "disagreement_names_2" ~ disagreement_perc_2,
                                 name == "disagreement_names_3" ~ disagreement_perc_3,
                                 name == "disagreement_names_4" ~ disagreement_perc_4,
                                 name == "disagreement_names_5" ~ disagreement_perc_5,
                                 name == "disagreement_names_6" ~ disagreement_perc_6,
                                 name == "disagreement_names_7" ~ disagreement_perc_7),
         actual_n = case_when(name == "disagreement_names_1" ~ disagreement_n_1,
                              name == "disagreement_names_2" ~ disagreement_n_2,
                              name == "disagreement_names_3" ~ disagreement_n_3,
                              name == "disagreement_names_4" ~ disagreement_n_4,
                              name == "disagreement_names_5" ~ disagreement_n_5,
                              name == "disagreement_names_6" ~ disagreement_n_6,
                              name == "disagreement_names_7" ~ disagreement_n_7)) %>% 
  select(chromosome, ID, perc_agree_all, value, actual_perc, actual_n) %>% 
  filter(is.na(actual_perc) == FALSE) %>% 
  mutate(translation  = case_when(value == "{file0,file1,file2} vs. {file3}" ~ "whatshap_disagrees",
                                  value == "{file0,file1,file3} vs. {file2}" ~ "fimpute_disagrees",
                                  value == "{file0,file2,file3} vs. {file1}" ~ "alphapeel_disagrees",
                                  value == "{file0} vs. {file1,file2,file3}" ~ "beagle_disagrees",
                                  TRUE ~ "even_split")) %>% 
  select(-value, -actual_n) %>% 
  mutate(perc_disagreement = 100 - as.numeric(perc_agree_all) ) %>% 
  mutate(actual_perc = as.numeric(str_remove(actual_perc, "%"))) %>% 
  mutate(contribution = (as.numeric(actual_perc)/perc_disagreement)*100) %>% 
  mutate(group = str_c(ID,chromosome, sep = "_")) %>% 
  mutate(translation = factor(translation, levels = c("alphapeel_disagrees", "beagle_disagrees", "fimpute_disagrees", "whatshap_disagrees", "even_split"), 
                              labels = c("Alphapeel", "Beagle", "FImpute3", "WhatsHap", "Multiple")))

# agreement <- read_delim("~/Downloads/phasing_agreement.txt", delim = "\t")

ggplot(agreement, aes(x = translation, y = actual_perc, colour = translation))+
  geom_boxplot(outlier.size = -1, colour = "black")+
  geom_point(position = position_jitter(width = 0.2))+
  labs(x = "Method",  y = "Disagreement caused") +
  ylim(-1,50)+
  scale_colour_viridis_d()+
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2))) 


#####################################
# Per individual disagreement rates #
#####################################

per_indiv_rates <- all_compared %>% 
  select(chromosome, ID, perc_agree_all) %>% 
  mutate(threshold = case_when(perc_agree_all > 98 ~ "pass",
                               is.na(perc_agree_all) == TRUE ~ "unknown",
                               TRUE ~ "fail")) %>% 
  group_by(ID, threshold) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = threshold, values_from = n) %>% 
  mutate(across(c(fail, pass, unknown), ~replace_na(., 0))) %>% 
  mutate(disagree_rate = fail/(fail+pass)*100)

gen1 <- tibble(ID = read_lines("/media/drive_6tb/projects/kakapo-genomics/cohort_gen1.txt"),
               gen = 1) 
gen2 <- tibble(ID = read_lines("/media/drive_6tb/projects/kakapo-genomics/cohort_gen2.txt"),
               gen = 2)
gen3 <- tibble(ID = read_lines("/media/drive_6tb/projects/kakapo-genomics/cohort_gen3.txt"),
               gen = 3)

indiv_gen <- bind_rows(gen1, gen2, gen3) %>% 
  mutate(gen = factor(gen, levels =c("1", "2", "3")))

disagree_by_gen <- left_join(per_indiv_rates, indiv_gen)

ggplot(disagree_by_gen, aes(x = gen, y = disagree_rate, colour = gen))+ 
  geom_boxplot(colour = "black", outlier.size = -1)+
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  ylim(-1,100)+
  labs(x = "Generation", y = "Disagreement (%)") +
  theme_classic()+
  scale_color_viridis_d()+
  theme(legend.position = "none",
        legend.text = element_text(size = rel(1.6)),
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA,
                                        colour = "black",
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.6)))


blah <- low_agreement %>% 
  mutate(remove_focal = as.numeric(perc_agree_all) + actual_perc) %>% 
  filter(remove_focal > 98) %>% 
  group_by(translation) %>% 
  summarise(n())


# File order
# file0 = beagle
# file1 = alphapeel
# file2 = fimpute
# file3 = whatshap



# wc_file <- read_file("Aparima_all_compared.log")
# 
# wc_file <- read_file("/cifs/biocldap/scratch/deardenlab/marissalelec/whatshap_compare/whatshap_compare_all_Gulliver.log")
# 
# 
# test <- str_split(wc_file, "Chromosome") %>% 
#   unlist()
# 
# 
# collated <- lapply(2:length(test), function(chr_chunk){
#   str_split(test[chr_chunk], "\n", simplify = TRUE) %>% 
#     trimws() %>% 
#     parse_wc_lines2(n_compared = 4)
# }) %>% 
#   reduce(bind_rows)
  



# # Testing stuff on a single sample.  Results files
# results <- list.files(path = res_dir, 
#                       pattern = "Aparima_pairwise_", full.names = TRUE) %>%
#   setNames(., sub("\\.txt$", "", basename(.))) %>% 
#   map(get_data) %>% 
#   reduce(bind_rows)
# 
# intersection <- read_delim(paste0(res_dir, "whatshap_compare_Aparima_alphapeel_first.txt"), delim = "\t")
# 
# 
# 
# ggplot(results, aes(x = both, y = all_switch_rate, colour = both))+ geom_point()
# 
# 
# ggplot(results, aes(x = both, y = covered_variants, colour = both)) + geom_point()
# 
# 
# ggplot(results, aes(x = chromosome, y = covered_variants, colour = both)) + geom_point()
# 
# 
# ggplot(results, aes(x = both, y = all_switchflip_rate, colour = both)) + geom_point()
# 
# ggplot(results, aes(x = both, y = blockwise_hamming_rate, colour = both)) + geom_point()

