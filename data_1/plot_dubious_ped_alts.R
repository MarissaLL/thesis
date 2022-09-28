#!/usr/bin/env Rscript

library(tidyverse)
library(plotly)
library(xtable)

###########
# GLOBALS #
###########

fmendel <- read_delim("output/04_test_dubious_trios/all_candidates_fmendel.txt", delim = " ", col_names = FALSE, trim_ws = TRUE) %>% 
  mutate(X1 = as.numeric(X1))

run_trios <- read_delim("output/04_test_dubious_trios/candidate_trios_test.txt", delim = " ", trim_ws = TRUE)

step_trio_res_tbl <-  read_delim("output/03_mendelian_error_distribution/two_backup/backup/combined_stepwise_trio_results.txt", delim = " ")


########
# MAIN #
########



# According to m_err.log there were 190020 genotyped loci. Calculate mendelian error rate from this.
fmendel_by_chick <- left_join(fmendel, run_trios, by = c(X1 = "family_id", X2 = "ID1", X3 = "ID")) %>% 
  mutate(m_err_rate = (X5/190020)*100)


# Plot mendelian error rates for all the tested trios
ggplot(fmendel_by_chick, aes(x = chick, y = m_err_rate)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2)) +
  labs(x = "Chick",
       y = "Mendelian errors (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))


# Merge results with the stepwise trio results and add groups to allow points to be coloured differently
stepwise_trio_results <- step_trio_res_tbl %>% 
  mutate(point_group = case_when(step == "True trio" & m_err_rate > 2 ~ "outliers",
                                  TRUE ~ "expected_dist"))

all_likely <- fmendel_by_chick %>%  
  rename(father.y = X2, mother.y = X3, ID = chick) %>% 
  mutate(point_group = "new", 
         step = paste0("alternative_", ID))


combined_results <- bind_rows(stepwise_trio_results, all_likely) %>% 
  mutate(final_colour = case_when(point_group == "expected_dist" ~ "expected_dist",
                                  point_group == "outliers" & ID == "Egilsay" ~ "egilsay",
                                  point_group == "outliers" & ID == "Heather" ~ "heather",
                                  point_group == "outliers" & ID == "Ian" ~ "ian",
                                  point_group == "outliers" & ID == "Taeatanga" ~ "taeatanga",
                                  point_group == "new" & ID == "Egilsay" ~ "egilsay",
                                  point_group == "new" & ID == "Heather" ~ "heather",
                                  point_group == "new" & ID == "Ian" ~ "ian",
                                  point_group == "new" & ID == "Taeatanga" ~ "taeatanga"))


# Rearrange the order of groups for the plot
combined_results$step <- factor(combined_results$step, levels = c("alternative_Egilsay", 
                                              "alternative_Heather", 
                                              "alternative_Ian", 
                                              "alternative_Taeatanga", 
                                              "True trio", 
                                              "Sire replaced", 
                                              "Dam replaced", 
                                              "Both replaced"), 
                      labels = c("Egilsay", 
                                 "Heather", 
                                 "Ian", 
                                 "Taeatanga", 
                                 "True trio", 
                                 "Sire replaced", 
                                 "Dam replaced", 
                                 "Both replaced"))

# Plot the candidates for the dubious trios against the stepwise trio results
ggplot(combined_results, aes(x = step , y = m_err_rate)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2),  aes(colour = final_colour)) +
  labs(x = "",
       y = "Mendelian errors (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14)) 
# +  geom_hline(yintercept=1.23)


# Run an ANOVA and post-hoc tukey test
anova_res <- aov(m_err_rate ~ step, data = combined_results)
plot(anova_res)
summary(anova_res)
TukeyHSD(anova_res)



# Plot the candidates for the dubious trios with plotly to interactively look at the points
plotly_plot <- ggplot(fmendel_by_chick, aes(x = chick, y = m_err_rate, colour = X2, size = X3)) +
  # geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2)) +
  labs(x = "Chick",
       y = "Mendelian errors (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))



ggplotly(plotly_plot, tooltip = c("colour", "size"))





########### 3.7K possibilities #################
# This one ran until it crashed  

fmendel <- read_delim("output/04_test_dubious_trios/one3724_fmendels.txt", delim = " ", trim_ws = TRUE, col_names = FALSE)
run_trios <- read_delim("output/04_test_dubious_trios/candidate_trios_test.txt", delim = " ", trim_ws = TRUE)
step_trio_res_tbl <-  read_delim("output/03_mendelian_error_distribution/deepvariant.r1_all_step_trios.txt", delim = " ")



fmendel_by_chick <- left_join(fmendel, run_trios, by = c(X1 = "family_id", X2 = "ID1", X3 = "ID")) %>% 
  mutate(m_err_rate = (X5/190020)*100)

ggplot(fmendel_by_chick, aes(x = chick, y = m_err_rate)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2)) +
  labs(x = "Chick",
       y = "Mendelian errors (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))+
  geom_hline(aes(yintercept = 1.2))


fmendel_by_chick %>% 
  filter(m_err_rate < 1.2) %>% select(chick, X3, X2) %>% 
  xtable()


################## Chick by chick



egilsay_fmendel <- read_delim("output/04_test_dubious_trios/egilsay_candidates_fmendel.txt", delim = " ", trim_ws = TRUE, col_names = FALSE) %>% 
  mutate(chick = "Egilsay", m_err_rate = (X5/189905)*100)

heather_fmendel <- read_delim("output/04_test_dubious_trios/heather_candidates_fmendel.txt", delim = " ", trim_ws = TRUE, col_names = FALSE) %>% 
  mutate(chick = "Heather", m_err_rate = (X5/189905)*100)

ian_fmendel <- read_delim("output/04_test_dubious_trios/ian_candidates_fmendel.txt", delim = " ", trim_ws = TRUE, col_names = FALSE) %>% 
  mutate(chick = "Ian", m_err_rate = (X5/189905)*100)

taeatanga_fmendel <- read_delim("output/04_test_dubious_trios/taeatanga_candidates_fmendel.txt", delim = " ", trim_ws = TRUE, col_names = FALSE) %>% 
  mutate(chick = "Taeatanga", m_err_rate = (X5/189905)*100)

elliott_fmendel <- read_delim("output/04_test_dubious_trios/elliott_candidates_fmendel.txt", delim = " ", trim_ws = TRUE, col_names = FALSE) %>% 
  mutate(chick = "Elliott", m_err_rate = (X5/189905)*100)

candidates <- bind_rows(egilsay_fmendel, heather_fmendel, ian_fmendel,  elliott_fmendel) # also taeatanga_fmendel,


ggplot(candidates, aes(x = chick, y = m_err_rate)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2)) +
  labs(x = "Chick",
       y = "Mendelian errors (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))+
  geom_hline(aes(yintercept = 1.2))


rimu_fmendel <- read_delim("output/04_test_dubious_trios/rimu_candidates_fmendel.txt", delim = " ", trim_ws = TRUE, col_names = FALSE) %>% 
  mutate(chick = "Rimu", m_err_rate = (X5/189905)*100)

jamieson_fmendel <- read_delim("output/04_test_dubious_trios/jamieson_candidates_fmendel.txt", delim = " ", trim_ws = TRUE, col_names = FALSE) %>% 
  mutate(chick = "Jamieson", m_err_rate = (X5/189905)*100)

konini_3416_fmendel <- read_delim("output/04_test_dubious_trios/konini_3-4-16_candidates_fmendel.txt", delim = " ", trim_ws = TRUE, col_names = FALSE) %>% 
  mutate(chick = "Konini_3-4-16", m_err_rate = (X5/189905)*100)

te_awa_fmendel <- read_delim("output/04_test_dubious_trios/te_awa_candidates_fmendel.txt", delim = " ", trim_ws = TRUE, col_names = FALSE) %>% 
  mutate(chick = "Te_Awa", m_err_rate = (X5/189905)*100)

paternity_candidates <-  bind_rows(rimu_fmendel, jamieson_fmendel, konini_3416_fmendel) %>% 
  rename(step = chick) %>% 
  mutate(point_group = "test")

paternity_and_step <- bind_rows(paternity_candidates, stepwise_trio_results) %>% 
  filter(point_group != "outliers") %>% 
  mutate(reordered = fct_relevel(step, c("Rimu", "Jamieson", "Konini_3-4-16", "True trio", "Sire replaced", "Dam replaced", "Both replaced")))

ggplot(paternity_candidates, aes(x = step, y = m_err_rate)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2)) +
  labs(x = "Chick",
       y = "Mendelian errors (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))+
  geom_hline(aes(yintercept = 1.2))


ggplot(paternity_and_step, aes(x = reordered, y = m_err_rate)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2)) +
  labs(x = "Trio",
       y = "Mendelian errors (%)") +
  ylim(0,20) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))





### All the potential mistakes 

dubious_final <- bind_rows(egilsay_fmendel, elliott_fmendel, heather_fmendel, ian_fmendel, taeatanga_fmendel, rimu_fmendel, jamieson_fmendel) 
 
ggplot(dubious_final, aes(x = chick, y = m_err_rate)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2)) +
  labs(x = "Chick",
       y = "Mendelian errors (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))+
  geom_hline(aes(yintercept = 0.12))


## Compare to step results


step_trio_res_tbl <-  read_delim("output/03_mendelian_error_distribution/two_backup/backup/combined_stepwise_trio_results.txt", delim = " ") %>% 
  rename(step_or_chick = step)

all_mistakes <- bind_rows(egilsay_fmendel, elliott_fmendel, heather_fmendel, ian_fmendel, rimu_fmendel, jamieson_fmendel) %>% 
  rename(step_or_chick = chick)

comparison <- bind_rows(step_trio_res_tbl, all_mistakes) %>% 
  mutate(colour_group = case_when(step_or_chick == "True trio" & m_err_rate > 2 ~ "outlier",
                                  TRUE ~ "normal")) 



plot_order <- c("Egilsay",
                "Elliott",
                "Heather",
                "Ian",
                "Taeatanga",
                "Rimu",
                "Jamieson",
                "True trio",
                "Father replaced",
                "Mother replaced",
                "Both replaced")

comparison$step_or_chick <-  factor(comparison$step_or_chick, levels = c("Rimu",
                                                                         "Jamieson",
                                                                         "Egilsay",
                                                                         "Elliott",
                                                                         "Heather",
                                                                         "Ian",
                                                                         "Taeatanga",
                                                                         "True trio",
                                                                         "Sire replaced",
                                                                         "Dam replaced",
                                                                         "Both replaced"),
                                    labels = c("Rimu", "Jamieson", "Egilsay", "Elliott", 
                                               "Heather", "Ian", "Taeatanga", "True trio",
                                               "Father replaced", "Mother replaced",
                                               "Both replaced"))



ggplot(comparison, aes(x = step_or_chick, y = m_err_rate)) +
  geom_boxplot(data = subset(comparison, step_or_chick != "Rimu" & step_or_chick != "Jamieson"), outlier.size = -1) +
  geom_point(data = subset(comparison, colour_group == "normal" & step_or_chick != "Rimu" & step_or_chick != "Jamieson"),position = position_jitter(width= 0.2), colour = "#452736") +
  geom_point(data = subset(comparison, colour_group == "outlier" & step_or_chick != "Rimu" & step_or_chick != "Jamieson"),position = position_jitter(width= 0.2), colour = "red") +
  labs(x = "Trio",
       y = "Mendelian errors (%)") +
  ylim(0,20)+
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14)) # +
  # geom_hline(aes(yintercept = 1.2))



#### For plotly to not crash
dubuious_final_small <- dubious_final %>%  
  filter(m_err_rate <3)
  

ah <- ggplot(dubious_final_small, aes(x = chick, y = m_err_rate, colour = X2, size = X3)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2)) +
  labs(x = "Chick",
       y = "Mendelian errors (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))+
  geom_hline(aes(yintercept = 1.2))


ggplotly(ah)


