#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(changepoint))
suppressPackageStartupMessages(library(tidyverse,
                                       warn.conflicts = FALSE))





args <- commandArgs(trailingOnly = TRUE)
chromosome <- args[1]
ncycles <- args[2]
inherited_chr <- args[3]
seq_birds_file <- args[4]
parent <- str_remove(string = inherited_chr, pattern = "_only")
message(paste("Processing chromosome", chromosome, "ncycles", ncycles, parent))


# chromosome <- "S4"
# ncycles <- 100

## This comes from alphapeel_seg_haplos_added.R
added_haplos <- read_delim(paste0("added_haplos_", chromosome,"_ncycles", ncycles,".txt"), delim = "\t")


all_seq <- read_delim(seq_birds_file, delim = "\t", col_names = FALSE) %>% pull() %>% str_sort()

# Birds that didn't get phased through to grandparents based on alphapeel_seg_haplos_added.R
unphased <- c("Alice", "Arab", "Barnard", "Basil", "Bella", "Ben", "Bill", 
              "Blades", "Bonus", "Boss", "Cyndy", "Felix", "Flossie", "Fuchsia", 
              "Gumboots", "Gunner", "Heather", "Jean", "Jimmy", "Joe", "John-girl", 
              "Lee", "Lionel", "Lisa", "Luke", "Maggie", "Margaret-Maree", 
              "Merv", "Nog", "Nora", "Ox", "Pegasus", "Piripi", "Ralph", "Rangi", 
              "Richard_Henry", "Ruth", "Sandra", "Sass", "Smoko", "Solstice", 
              "Stumpy", "Sue", "Suzanne", "Waynebo", "Whiskas")


first_test <- all_seq[all_seq %in% unphased == FALSE]


phased_only <- added_haplos %>% 
  filter(chick %in% first_test == TRUE)

only_095 <- phased_only %>% 
  filter(sum_val >= 0.95)

### Remove low quality phasing individuals ###

# Test filter - fewer than 10% of the SNPs that the best has on that chr - all exclusions seem reasonable
# 20% seems a little harsh, though there definitely are a fair few bad ones in there. Maybe they can be excluded through other means though
num_snps_095 <- only_095 %>% 
  filter(chick %in% first_test == TRUE) %>% 
  filter(focal_parent == inherited_chr) %>% 
  group_by(chick) %>% 
  summarise(nsnp = n())


most <- max(num_snps_095$nsnp)

filtered_out <- num_snps_095 %>% 
  filter(nsnp < most/5)

miss_filter <- num_snps_095 %>% 
  mutate(missing_filter = case_when(nsnp < most/5 ~ "fail",
                                    nsnp > most/5 ~ "pass"))


## % uncertain in the sliding window
# seems like some of the lowest uncertainty individuals are just the ones with hardly any snps, even though the % somewhat accounts for thats
find_perc_uncertain <-  function(indiv, focal_parent){
  aaa <- only_095 %>% 
    filter(focal_parent == !!(focal_parent)) %>% 
    filter(chick == !!(indiv)) %>% 
    mutate(inherited2 = case_when(inherited == "mat" ~ 1,
                                  inherited == "pat" ~ 0)) %>% 
    group_by(chick) %>% 
    arrange(pos) %>%
    mutate(lagged_mean = rollmean(x = inherited2, 50, align = "right", fill = NA)) %>% 
    mutate(good_rgn_lagged = case_when(lagged_mean <= 0.2 ~ "paternal",
                                       lagged_mean >= 0.8 ~ "maternal",
                                       TRUE ~ "unclear"))
  
  bbb <- tibble(indiv = indiv,
                unclear = length(aaa$good_rgn_lagged[aaa$good_rgn_lagged == "unclear"]),
                clear = length(aaa$good_rgn_lagged[aaa$good_rgn_lagged != "unclear"]),
                chromosome = chromosome,
                focal_parent = focal_parent)
  return(bbb)
}

perc_uncertain <- lapply(first_test, find_perc_uncertain, inherited_chr) %>% 
  reduce(bind_rows) %>% 
  mutate(perc_unclear = unclear/(clear + unclear))


### Number of phase switches in the sw data, or mean length of phase block
### Why is this completely different from the S9 plots from the same indivs 

find_total_ps <- function(indiv, focal_parent){
  aaa <- only_095 %>% 
    filter(focal_parent == !!(focal_parent)) %>% 
    filter(chick == indiv) %>% 
    mutate(inherited2 = case_when(inherited == "mat" ~ 1,
                                  inherited == "pat" ~ 0)) %>% 
    group_by(chick) %>% 
    arrange(pos) %>%
    mutate(lagged_mean = rollmean(x = inherited2, 50, align = "right", fill = NA)) %>% 
    mutate(good_rgn_lagged = case_when(lagged_mean <= 0.2 ~ "paternal",
                                       lagged_mean >= 0.8 ~ "maternal",
                                       TRUE ~ "unclear")) %>% 
    arrange(pos)
  
  # rle_lengths <- rle(aaa$good_rgn_lagged)[[1]][rle(aaa$good_rgn_lagged)[[1]] > 50]
  bbb <- tibble(indiv = indiv,
                tot_ps = length(rle(aaa$good_rgn_lagged)[[1]]),
                small_ps = length(rle(aaa$good_rgn_lagged)[[1]][rle(aaa$good_rgn_lagged)[[1]] < 50]),
                large_ps = length(rle(aaa$good_rgn_lagged)[[1]][rle(aaa$good_rgn_lagged)[[1]] > 100]),
                med_ps = length(rle(aaa$good_rgn_lagged)[[1]][rle(aaa$good_rgn_lagged)[[1]] > 50 & rle(aaa$good_rgn_lagged)[[1]] < 100]),
                mean_psl = mean(rle(aaa$good_rgn_lagged)[[1]]),
                chromosome = chromosome,
                focal_parent = focal_parent)
}


ps <- lapply(first_test, find_total_ps, inherited_chr) %>% 
  reduce(bind_rows)


all_info <- full_join(miss_filter, perc_uncertain, by = c(chick = "indiv")) %>% 
  full_join(ps, by = c(chromosome = "chromosome", focal_parent = "focal_parent", chick = "indiv"))


### With missing. Could possibly do with being a little harsher on the %unclear or mean_psl  
failed <- all_info %>% 
  mutate(fails = case_when(missing_filter == "fail" ~ "fail",
                           (perc_unclear > 0.75) == TRUE ~ "fail",
                           (mean_psl < 150) == TRUE ~ "fail",
                           is.na(nsnp) == TRUE & tot_ps == 0 ~ "fail",
                           TRUE ~ "pass"))

write_delim(failed, paste0("filter_res_", chromosome, "_ncycles_", ncycles, "_", parent, ".txt"))

message(paste0(table(failed$fails)[1], " failed. ", table(failed$fails)[2], " passed."))

indivs_to_keep <- failed %>% 
  filter(fails == "pass") %>% 
  select(chick)

write_delim(indivs_to_keep, paste0("indivs_to_keep_", chromosome, "_ncycles_", ncycles, "_", parent, ".txt"))
