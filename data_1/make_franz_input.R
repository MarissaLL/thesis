library(tidyverse)

cmdargs <- commandArgs(trailingOnly = TRUE)


file <- cmdargs[1]

outfile <- paste0(str_extract(file, "^[:graph:]*(?=raw)"), "franz")

pedigree <- read_delim("bird_info_rships_nov20.csv", delim = ",") %>% 
  select(ID, sex)

data <- read_delim(file, delim = " ", skip = 1, col_names = FALSE) %>%
  #select(-c(FID, SEX, PHENOTYPE, MAT, PAT)) %>% 
  select(-c(X2,X3, X4, X5,X6)) %>%
  mutate(across(-c(X1), ~ case_when(.x == 0 ~ "0/0",
                                           .x == 1 ~"0/1",
                                           .x == 2 ~ "1/1",
                                           TRUE ~ "?/?"))) %>% 
  left_join(pedigree, by =c(X1 = "ID")) %>% 
  mutate(ones_col = 1, birth = "?", death = "?") %>% 
  select(X1, ones_col, birth, death, sex, everything()) %>% 
  mutate(X1 = str_c(X1, "              ")) %>% 
  mutate(X1 = str_trunc(X1, width = 10, ellipsis = ""))  %>% 
  write_delim(outfile, col_names = FALSE)

