library(tidyverse)

#############
# FUNCTIONS #
#############

decode_haplotypes <- function (indiv){
  coded_haplos %>% 
    mutate(!!paste0(sym(indiv),"_pat") := case_when(!!sym(indiv) == 0 ~ "A1",
                                                    !!sym(indiv) == 1 ~ "A3",
                                                    !!sym(indiv) == 2 ~ "A2",
                                                    !!sym(indiv) == 3 ~ "A1",
                                                    !!sym(indiv) == 4 ~ "A2",
                                                    !!sym(indiv) == 5 ~ "miss",
                                                    !!sym(indiv) == 6 ~ "A1",
                                                    !!sym(indiv) == 7 ~ "A2",
                                                    !!sym(indiv) == 8 ~ "miss",
                                                    !!sym(indiv) == 9 ~ "miss")) %>% 
    mutate(!!paste0(sym(indiv),"_mat") := case_when(!!sym(indiv) == 0 ~ "A1",
                                                    !!sym(indiv) == 1 ~ "A3",
                                                    !!sym(indiv) == 2 ~ "A2",
                                                    !!sym(indiv) == 3 ~ "A2",
                                                    !!sym(indiv) == 4 ~ "A1",
                                                    !!sym(indiv) == 5 ~ "miss",
                                                    !!sym(indiv) == 6 ~ "miss",
                                                    !!sym(indiv) == 7 ~ "miss",
                                                    !!sym(indiv) == 8 ~ "A1",
                                                    !!sym(indiv) == 9 ~ "A2"))
}

check_if_haplos_are_same <- function(combo){
  combo_part1 <- str_extract(combo, "[:graph:]*(?=_._)")
  combo_part2 <- str_extract(combo, "(?<=_._)[:graph:]*")
  comparison <- combined %>% 
    select(combo_part1, combo_part2) %>% 
    mutate(matching = case_when(!!sym(combo_part1) == !!sym(combo_part2) ~ "same",
                                !!sym(combo_part1) != !!sym(combo_part2) ~ "diff"))
  paste(sum(str_count(comparison$matching, "same")),
        sum(str_count(comparison$matching, "diff")), sep = "___")
}

return_haplo_differences <- function(combo){
  combo_part1 <- str_extract(combo, "[:graph:]*(?=_._)")
  combo_part2 <- str_extract(combo, "(?<=_._)[:graph:]*")
  comparison <- combined %>% 
    select(combo_part1, combo_part2) %>% 
    mutate(matching = case_when(!!sym(combo_part1) == !!sym(combo_part2) ~ "same",
                                !!sym(combo_part1) != !!sym(combo_part2) ~ "diff"))
  return(comparison)
}


########
# MAIN #
########

# Just doing this for CM013779.1 because it has relatively few snps, and a low rate of errors (89 out of 1407) compared to some of the other chromosomes.


coded_haplos <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/10_fimpute_phasing/output_CM013779.1/genotypes_imp.txt", delim = ' ', trim_ws = TRUE) %>% 
  filter(ID == "Hoki" | ID == "Zephyr" | ID == "Felix") %>% 
  separate(col = `Calls...`, into = paste0("pos_", c(1:1407)), sep = c(1:1406)) %>% 
  select(-Chip) %>% 
  t() %>% 
  as_tibble(.name_repair = "unique") %>% 
  rename(Felix = `...1`, Hoki = `...2`, Zephyr = `...3`) %>% 
  filter(Felix != "Felix")


hoki_haplos <- decode_haplotypes("Hoki") %>% 
  rownames_to_column()
zephyr_haplos <- decode_haplotypes("Zephyr") %>% 
  rownames_to_column()
felix_haplos <- decode_haplotypes("Felix") %>% 
  rownames_to_column()


combined <- left_join(hoki_haplos, zephyr_haplos, by = "rowname") %>% 
            left_join(felix_haplos, by = "rowname") %>% 
  select(-ends_with(".x"), -ends_with(".y"))



# Ideas 
# Use the combination thing to come up with all pairwise combos of indiv_mat, and indiv_pat for all indivs
# Case_when same ~ same, diff ~ diff, count number of diffs

indiv <-  c("Hoki", "Zephyr", "Felix")
suffix <- c("_mat", "_pat")


combinations_to_check <- crossing(haplo_1 = paste(rep(indiv, each = length(suffix)), suffix, sep = ""), 
                                  haplo_2 = paste(rep(indiv, each = length(suffix)), suffix, sep = "")) %>% 
  mutate(combo = str_c(haplo_1, haplo_2, sep = "_._")) %>% 
  select(combo) %>% 
  pull()


check <- check_if_haplos_are_same("Felix_pat_._Hoki_pat")

############# Count the number of SNPs that are different between the maternally inherited haplotype in the chick and the maternal haplotype it is derived from (same for paternal)
haplos_same_diff <-  lapply(combinations_to_check, check_if_haplos_are_same) %>% 
  unlist() 


summarise_diffs <- enframe(x = haplos_same_diff, value = "matches") %>% 
  mutate(combos = combinations_to_check) %>% 
  select(matches, combos) %>% 
  separate(matches, 
           into = c("same", "diff"),
           sep = "___") %>% 
  mutate(same = as.numeric(same),
         diff = as.numeric(diff))




############# Return the actual haplotype differences. Need to check if I am doing this again that there are no differences where one of them is just missing

# haplo_diff <- return_haplo_differences("Hoki_mat_._Zephyr_pat") %>% 
#   rownames_to_column() %>% 
#   filter(matching == "diff") #60
# 
# haplo_diff <- return_haplo_differences("Hoki_pat_._Felix_mat") %>%
#   rownames_to_column %>%
#   filter(matching == "diff") #299

haplo_diff_mat <- return_haplo_differences("Hoki_mat_._Zephyr_pat") %>% 
  rownames_to_column()

haplo_diff_pat <- return_haplo_differences("Hoki_pat_._Felix_mat") %>% 
  rownames_to_column

ggplot(haplo_diff_mat, aes(x = matching, y = rowname)) + 
  geom_point(shape = "_") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggplot(haplo_diff_pat, aes(x = matching, y = rowname)) + 
  geom_point(shape = "_") +
 
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

hoki_mat_c <- str_c(combined$Hoki_mat,collapse = "")
zephyr_mat_c <- str_c(combined$Zephyr_pat,collapse = "")






# See if haplo_diff_pat differences instead match the other paternal haplotype
is_it_the_other_one_pat <- combined %>% 
  select(Hoki_pat, Felix_mat, Felix_pat) %>% 
  mutate(matching_exp = case_when(Hoki_pat == Felix_mat ~ "same",
                              Hoki_pat != Felix_mat ~ "diff")) %>% 
  mutate(matching_other = case_when(Hoki_pat == Felix_pat ~ "same",
                                    Hoki_pat != Felix_pat ~ "diff")) %>% 
  mutate(all_same = case_when(matching_exp == "same" & matching_other == "same" ~ "all_identical",
                              matching_exp == "same" & matching_other == "diff" ~ "only_exp",
                              matching_exp == "diff" & matching_other == "same" ~ "only_other",
                              matching_exp == "diff" & matching_other == "diff" ~ "non_match")) %>% 
  rownames_to_column() %>% 
  mutate(rowname = as.numeric(rowname)) %>% 
  mutate(all_same = fct_relevel(all_same, c("all_identical","only_exp", "only_other", "non_match")))



ggplot(is_it_the_other_one_pat, aes(x = matching_exp, y = rowname, colour = all_same)) + 
  geom_point(shape = "_", size = 10) +
  labs(colour = "Matching with haplotype",
       x = "Paternal") +
  theme_classic() +
  theme(legend.position = "blank",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text()) +
  scale_colour_manual(values = c("#ca766aff", "#dbcbbaff", "#452736ff", "#6997b6ff"))




is_it_the_other_one_mat <- combined %>% 
  select(Hoki_mat, Zephyr_mat, Zephyr_pat) %>% 
  mutate(matching_exp = case_when(Hoki_mat == Zephyr_pat ~ "same",
                                  Hoki_mat != Zephyr_pat ~ "diff")) %>% 
  mutate(matching_other = case_when(Hoki_mat == Zephyr_mat ~ "same",
                                    Hoki_mat != Zephyr_mat ~ "diff")) %>% 
  mutate(all_same = case_when(matching_exp == "same" & matching_other == "same" ~ "all_identical",
                              matching_exp == "same" & matching_other == "diff" ~ "only_exp",
                              matching_exp == "diff" & matching_other == "same" ~ "only_other",
                              matching_exp == "diff" & matching_other == "diff" ~ "non_match")) %>% 
  rownames_to_column() %>% 
  mutate(rowname = as.numeric(rowname)) %>% 
  mutate(all_same = fct_relevel(all_same, c("all_identical","only_exp", "only_other", "non_match")))



ggplot(is_it_the_other_one_mat, aes(x = matching_exp, y = rowname, colour = all_same)) + 
  geom_point(shape = "_", size = 10) +
  labs(colour = "Matching with haplotype",
       x = "Maternal") +
  theme_classic() +
  theme(legend.position = "blank",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_colour_manual(values =  c("#ca766aff", "#dbcbbaff", "#452736ff", "#6997b6ff"))



# Not separated by same/diff

ggplot(is_it_the_other_one_pat, aes(x = "Paternal", y = rowname, colour = all_same)) + 
  geom_point(shape = "_", size = 10) +
  labs(colour = "Matching with haplotype",
       x = "Paternal") +
  theme_classic() +
  theme(legend.position = "blank",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text()) +
  scale_colour_manual(values = c("#ca766aff", "#dbcbbaff", "#452736ff", "#6997b6ff"))


ggplot(is_it_the_other_one_mat, aes(x = "Maternal", y = rowname, colour = all_same)) + 
  geom_point(shape = "_", size = 10) +
  labs(colour = "Matching with haplotype",
       x = "Maternal") +
  theme_classic() +
  theme(legend.position = "blank",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_colour_manual(values =  c("#ca766aff", "#dbcbbaff", "#452736ff", "#6997b6ff"))



# 0 & A1A1 \\
# 1 & Unphased heterozygous\\
# 2 & A2A2\\
# 3 & A1A2\\
# 4 & A2A1\\
# 5 & missing\\
# 6 & A1–\\
# 7 & A2–\\
# 8 & –A1\\
# 9 & –A2\\
# 
# 
# decoding_haplos <- coded_haplos %>% 
#   mutate(Felix_pat = case_when(Felix == 0 ~ "A1",
#                                Felix == 1 ~ "A3",
#                                Felix == 2 ~ "A2",
#                                Felix == 3 ~ "A1",
#                                Felix == 4 ~ "A2",
#                                Felix == 5 ~ "miss",
#                                Felix == 6 ~ "A1",
#                                Felix == 7 ~ "A2",
#                                Felix == 8 ~ "miss",
#                                Felix == 9 ~ "miss")) %>% 
#   mutate(Felix_mat = case_when(Felix == 0 ~ "A1",
#                                Felix == 1 ~ "A3",
#                                Felix == 2 ~ "A2",
#                                Felix == 3 ~ "A2",
#                                Felix == 4 ~ "A1",
#                                Felix == 5 ~ "miss",
#                                Felix == 6 ~ "miss",
#                                Felix == 7 ~ "miss",
#                                Felix == 8 ~ "A1",
#                                Felix == 9 ~ "A2"))
