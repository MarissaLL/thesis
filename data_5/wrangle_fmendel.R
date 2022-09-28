library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(readr)
library(purrr)

wrangle_results <- function(file){
  res <- read_delim(file, delim = "\t", col_names = FALSE) %>%
  filter(X1 != "--") %>% 
  mutate(num = rep_len(c("rate","errors"), n_combos*2)) %>%
  mutate(num2 = rep(1:n_combos, each =2)) %>% 
  pivot_wider(names_from =num, values_from = X1) %>% 
  mutate(parents = str_extract(rate, "[:alpha:]+__[:alpha:]+")) %>% 
  mutate(geno_rate = str_extract(rate, "(?<=is\\s)[:graph:]+(?=.$)")) %>% 
  mutate(errors = as.numeric(str_extract(errors, "(?<=:\\s)[:digit:]+(?=\\sMendel)"))) %>% 
  select(parents, geno_rate, errors) %>% 
  mutate(genotyped_variants = 1747971 * as.numeric(geno_rate)) %>% 
  mutate(error_rate = as.numeric(errors)/genotyped_variants) %>%
  arrange(errors)
}

n_combos <- 312 # WH parents
#n_combos <- 8393 # All parents
#n_combos <- 9810


args <- commandArgs(trailingOnly = TRUE)
fmendel_file <- args[1]

file_basename <- str_remove(fmendel_file, ".fmendel")

wrangle_results(fmendel_file) %>%
	write_delim(paste0(file_basename, ".txt"), delim = "\t")
 
