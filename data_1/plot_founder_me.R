library(tidyverse)
library(patchwork)

setwd("/media/drive_6tb/projects/kakapo-genomics/")

fmendel <- read_delim("output/041_test_dubious_trios_2020/all_candidates_fmendel_allfounders.txt", delim = " ", col_names = FALSE, trim_ws = TRUE) %>% 
  mutate(X1 = as.numeric(X1))

run_trios <- read_delim("output/041_test_dubious_trios_2020/candidate_trios_allfounders.txt", delim = " ", trim_ws = TRUE)

fmendel_true_trios <- read_delim("output/041_test_dubious_trios_2020/all_candidates_fmendel_step.txt", delim = " ", col_names = FALSE, trim_ws = TRUE) %>% 
  mutate(X1 = as.numeric(X1))

run_true_trios <- read_delim("output/041_test_dubious_trios_2020/candidate_trios_step.txt", delim = " ", trim_ws = TRUE)

# me_summary_trios <- read_delim("output/041_test_dubious_trios_2020/logs/mendelian_errors_summary_step.txt", delim = "\t")


# Note that M.E were also calculated for RH but have been excluded because he is known to be from a different population.
founder_set_1 <- c("Alice", "Arab", "Barnard", "Basil", "Bella", "Ben", "Bill", 
                   "Blades", "Bonus", "Boss", "Cyndy", "Felix", "Flossie", "Fuchsia", 
                   "Gumboots", "Gunner", "Jean", "Jimmy", "Joe", "Lee", "Lionel", 
                   "Lisa", "Luke")
founder_set_2 <- c("Maggie", "Margaret-Maree", "Merty", "Merv", "Nog", "Nora", 
                   "Ox", "Piripi", "Ralph", "Rangi", "Ruth", "Sandra", 
                   "Sarah", "Sass", "Smoko", "Solstice", "Sue", "Suzanne", "Waynebo", 
                   "Wendy", "Whiskas")



# According to m_err.log there were 200000 genotyped loci. Calculate mendelian error rate from this.
fmendel_by_chick <- left_join(fmendel, run_trios, by = c(X1 = "family_id")) %>% 
  mutate(m_err_rate = (X5/200000)*100) %>% 
  mutate(plot_grp = case_when(chick %in% founder_set_1 == TRUE ~ 1,
                               chick %in% founder_set_2 == TRUE ~ 2,
                              TRUE ~ 3))



fmendel_true <- left_join(fmendel_true_trios, run_true_trios, by = c(X1 = "family_ID")) %>% 
  mutate(m_err_rate = (X5/200000)*100) %>% 
  filter(step == "real") %>%
  mutate(chick = "Known trios")



fmendel_all1  <- bind_rows(fmendel_by_chick,fmendel_true) %>% 
  filter(is.na(plot_grp) == TRUE | plot_grp == 1) %>% 
  mutate(plot_grp = 1)

fmendel_all1$chick = factor(fmendel_all1$chick, levels = c("Known trios", "Alice", "Arab", "Barnard", "Basil", "Bella", "Ben", "Bill", 
                                             "Blades", "Bonus", "Boss", "Cyndy", "Felix", "Flossie", "Fuchsia", 
                                             "Gumboots", "Gunner", "Jean", "Jimmy", "Joe", "Lee", "Lionel", 
                                             "Lisa", "Luke"))

fmendel_all2  <- bind_rows(fmendel_by_chick, fmendel_true) %>% 
  filter(is.na(plot_grp) == TRUE | plot_grp == 2) %>% 
  mutate(plot_grp = 2)

fmendel_all2$chick = factor(fmendel_all2$chick, levels = c("Known trios", "Maggie", "Margaret-Maree", "Merty", "Merv", "Nog", "Nora", 
                                                           "Ox", "Piripi", "Ralph", "Rangi", "Ruth", "Sandra", 
                                                           "Sarah", "Sass", "Smoko", "Solstice", "Sue", "Suzanne", "Waynebo", 
                                                           "Wendy", "Whiskas"))

  
a <- ggplot(fmendel_all1, aes(x = chick, y = m_err_rate)) +
  # geom_hline(yintercept = 0.6, colour = "red")+
  # geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2), colour = "#452736") +

  labs(x = NULL,
       y = "Mendelian errors (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

b <- ggplot(fmendel_all2, aes(x = chick, y = m_err_rate)) +
  # geom_hline(yintercept = 0.6, colour = "red")+
  # geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2), colour = "#452736") +
  
  labs(x = "Founder",
       y = "Mendelian errors (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

a / b


