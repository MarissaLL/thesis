#!/usr/bin/env Rscript

setwd("~/projects/kakapo-genomics/")

library(tidyverse)
library(xtable)

trios_run <- read_delim("deepvariant.r1_all_trios.txt", delim = " ")

mendelian_error_results <- read_delim("prev_me_dists/all_fmendel_5k200k.txt", delim = " ", col_names = FALSE) %>% 
  rename(FID = X1, father = X2, mother = X3, num_offspring = X4, num_m_err = X5) %>% 
  mutate(FID = as.numeric(FID)) %>% 
  distinct(FID,.keep_all = TRUE)

set_to_zero <- read_delim("prev_me_dists/set_to_zero_5k_200k.txt", delim = " ", col_names = FALSE) %>% 
  select(X1) %>% 
  mutate(FID = 1:5000) %>% 
  rename(num_set_zero = X1)

me_res <- inner_join(mendelian_error_results, set_to_zero, by = c("FID" = "FID")) %>% 
  mutate(tot_var = 189938- num_set_zero) %>% 
  mutate(m_err_rate = (num_m_err/tot_var)*100)


mendelian_error_results %>% 
  filter(as.numeric(num_offspring) != 1)


ggplot(me_res, aes(x = m_err_rate)) +
  geom_histogram(binwidth = 0.5)+
  labs(x = "Mendelian errors (%)",
       y = "Number of trios") 






##### Added later ############ 
mean <- mean(me_res$m_err_rate)
sd1 <- mean - sd(me_res$m_err_rate)
sd2 <- mean - 2* sd(me_res$m_err_rate)
sd3 <- mean + sd(me_res$m_err_rate)
sd4 <- mean + 2* sd(me_res$m_err_rate)

ggplot(me_res, aes(x = m_err_rate)) +
  geom_histogram(binwidth = 0.5)+
  labs(x = "Mendelian errors (%)",
       y = "Number of trios") +
  geom_vline(xintercept = mean, colour = "blue")+
  geom_vline(xintercept = sd1)+
  geom_vline(xintercept = sd2)+
  geom_vline(xintercept = sd3)+
  geom_vline(xintercept = sd4)


##########################################################################

# These are not from the same run, and therefore family ids that are the same do not correspond to the same families
# step_trios_run <- read_delim("output/03_mendelian_error_distribution/two_backup/backup/all_fmendel.txt", delim = " ", col_names = FALSE)
# 
# me_step_res <- read_delim("stepwise_all_fmendel.txt", delim = " ", col_names = FALSE) %>% 
#   rename(FID = X1, father = X2, mother = X3, num_offspring = X4, num_m_err = X5) %>% 
#   mutate(FID = as.numeric(FID)) %>% 
#   distinct(FID,.keep_all = TRUE)
#   
# 
# # step_set_to_zero <- read_delim("stepwise_set0.txt", delim = " ", col_names = FALSE) %>% 
# #   select(X1) %>% 
# #   rename(num_set_zero = X1) %>% 
# #   mutate(FID = 1:510)
# 
# step_trio_res <- left_join(me_step_res, step_trios_run, by = c(FID = "family_ID")) %>%
#     mutate(m_err_rate = (num_m_err/190020)*100) 
#     
# 
# 
# # Don't know if I should use this because I'm not confident of the order
# # me_step_res <- inner_join(me_step_results, step_set_to_zero, by = c("FID" = "FID")) %>% 
# #   mutate(tot_var = 189938 - num_set_zero) %>% 
# #   mutate(m_err_rate = (num_m_err/tot_var)*100)
# 
# tidy_step_export <- left_join(me_step_res, step_trios_run, by = c("FID" = "family_ID", "father" = "father", "mother" = "mother"))
# 
# 
# me_step_res %>% 
#   filter(as.numeric(num_offspring) != 1)
# 
# 
# ggplot(step_trio_res, aes(x = m_err_rate)) +
#   geom_histogram(binwidth = 0.05)+  
#   labs(x = "Mendelian errors (%)",
#        y = "Number of trios")
# 
# 
# step_trio_res$step <- factor(step_trio_res$step, levels = c("real",  "repl_father", "repl_mother", "repl_both"), labels = c("True trio", "Father replaced", "Mother replaced", "Both replaced"))
# 
# write_delim(step_trio_res, "output/03_mendelian_error_distribution/combined_stepwise_trio_results.txt")
step_trio_res <- read_delim("output/03_mendelian_error_distribution/two_backup/backup/combined_stepwise_trio_results.txt", delim = " ")

############# This is new #########

# I'm not sure if the mean and standard deviation really works. Joseph suggests rank first
# Maybe 
step_res <- step_trio_res %>% 
  filter(step == "True trio") %>% 
  select(m_err_rate)

themean <- mean(step_res$m_err_rate)
sd1 <- themean - sd(step_res$m_err_rate)
sd2 <- themean - 2* sd(step_res$m_err_rate)
sd3 <- themean + sd(step_res$m_err_rate)
sd4 <- themean + 2* sd(step_res$m_err_rate)

ggplot(step_res, aes(x = m_err_rate)) +
  geom_histogram(binwidth = 0.05) +
  labs(x = "Mendelian errors (%)",
       y = "Number of trios") +
  geom_vline(xintercept = themean, colour = "blue")+
  geom_vline(xintercept = c(sd1, sd2, sd3, sd4))

########### Add z-scores ##################

themean <- mean(step_res$m_err_rate)
sd1 <- sd(step_res$m_err_rate)
z_res <- step_res %>%
  mutate(mean = themean, sd = sd1) %>% 
  mutate(z_score = (m_err_rate - mean)/sd)


###########################################################




ggplot(step_trio_res, aes(x = step , y = m_err_rate)) +
  geom_boxplot() +
  geom_point(position = "jitter") +
  labs(x = "Type",
       y = "Mendelian errors (%)")


# Check significance of the differences between groups
anova_res <- aov(num_m_err ~ step, data=step_trio_res)
plot(anova_res)
summary(anova_res)
c <- TukeyHSD(anova_res)



check_outliers <- step_trio_res %>% 
  mutate(change = case_when(step == "True trio" & m_err_rate > 2 ~ "bright",
                            TRUE ~ "dim")) 

check_outliers$step <- factor(check_outliers$step, levels = c("True trio", "Sire replaced", "Dam replaced", "Both replaced"), labels = c("True trio", "Father replaced", "Mother replaced", "Both replaced"))






ggplot(check_outliers, aes(x = step , y = m_err_rate)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2),  aes(colour = change)) +
  labs(x = "Trio type",
       y = "Mendelian errors (%)") +
  theme_classic() +
  scale_colour_manual(values =  c("red","#452736")) +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))







dubious_real_trios <-  step_trio_res %>% 
  filter(step == "True trio" & m_err_rate > 2) %>% 
  select(FID, ID, mother.x, father.x,num_m_err, m_err_rate) %>% 
  rename(mother = mother.x, father = father.x)

# write_tsv(dubious_real_trios, "dubious_real_trios.txt")
  

xtable(dubious_real_trios)


names <- dubious_real_trios %>% 
  select(mother, father, ID) %>% 
  gather() %>% 
  select(value) %>% 
  unique() %>% 
  as.vector()

names2 <- c("Suzanne", "Flossie", "Sue", "Arab",  "Ben", "Boss",  "Basil", "Egilsay", "Heather", "Ian", "Taeatanga"   )

kin <- read_delim("~/Downloads/plink2.kin0", delim = "\t") %>% 
  rename(ID1 = `#ID1`) %>% 
  filter(ID1 %in% names2) 

kin2 <-  kin %>% 
  group_by(ID1) %>% 
  group_split()


list2env(kin2, envir=.GlobalEnv)


not_incl <- kin %>% 
  filter(ID2 %in% c("Egilsay", "Ian", "Heather", "Taetanga")) %>% 
  filter(KINSHIP > 0.2)

pedigree_data <- read_csv("old_files/formatted_bird_list_manual.csv") %>% 
  mutate(ID = str_replace_all(ID, " ", "_")) %>% 
  mutate(mother = str_replace_all(mother, " ", "_")) %>% 
  mutate(father = str_replace_all(father, " ", "_"))



# Looking at the four trios that seemed wrong: Egilsay, Heather, Ian, Taeatanga. Generating a table of indivs with kinship over 0.2 for each

make_table <-  function(indiv, section) {
  indiv_name <- bind_rows(section, not_incl) %>% 
    filter(ID1 == indiv | ID2 == indiv) %>% 
    filter(KINSHIP > 0.177) %>% 
    arrange(desc(KINSHIP)) %>% 
    mutate(ID = if_else(ID1 == indiv, ID2, ID1)) %>% 
    select(ID, KINSHIP)
  
  
  indiv_info <- inner_join(indiv_name, pedigree_data, by = c(ID = "ID")) %>%
    select(ID, hatch_year, sex, KINSHIP)


  print(xtable(indiv_info,
               align =  c("l","l", "l", "l", "l"),
               digits = c(0,0,0,0,4),
               caption = indiv),
        include.rownames = FALSE)

}

make_table("Egilsay", `5`)
make_table("Heather", `7`)
make_table("Ian", `8`)
make_table("Taeatanga", `11`)

# Could make a thing where I assign the number of errors that could be the result of the mother being wrong, the father being wrong or both?


to_find <-  pedigree_data %>% 
  filter(father == "Unk")

####################### Another go with step.nf output

fmendel_stepnf <- read_delim("output/03_mendelian_error_distribution/all_stepwise_fmendel.txt", delim = " ", trim_ws = TRUE, col_names = FALSE)
trios_stepnf <-  read_delim("output/03_mendelian_error_distribution/deepvariant.r1_all_step_trios.txt", delim = " ", trim_ws = TRUE)

fmendel_res <- fmendel_stepnf %>% 
  rename(family_ID = X1, father = X2, mother = X3, num_offspring = X4, num_m_err = X5)

all_step_trios <- left_join(fmendel_res, trios_stepnf) %>% 
  mutate(m_err_rate = (num_m_err/94949)*1000)

all_step_trios$step <- factor(all_step_trios$step, levels = c("real",  "repl_father", "repl_mother", "repl_both"), labels = c("True trio", "Father replaced", "Mother replaced", "Both replaced"))



ggplot(all_step_trios, aes(x = step, y = m_err_rate)) +
  geom_point(aes(subset(all_step_trios, ID == "Zephyr" & step == "True trio")))+
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2)) +
  labs(x = "Type",
       y = "Mendelian errors (%)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))



anova_res <- aov(num_m_err ~ step, data=all_step_trios)
plot(anova_res)
summary(anova_res)
TukeyHSD(anova_res)





ggplot(all_step_trios, aes(x = m_err_rate)) +
  geom_histogram(binwidth = 0.05)+
  labs(x = "Mendelian errors (%)",
       y = "Number of trios") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))


all_step_trios %>% 
  filter(step == "True trio" & m_err_rate > 2)


# Colour the five outliers differently 
check_outliers <- all_step_trios %>% 
  mutate(point_group = case_when(step == "True trio" & m_err_rate > 2  ~ "outliers",
                            TRUE ~ "expected_dist")) %>% 
  filter(ID != "Taeatanga")

ggplot(subset(check_outliers, point_group == 'expected_dist'), aes(x = step , y = m_err_rate)) +
  geom_boxplot(outlier.size = -1) +
  geom_point(position = position_jitter(width= 0.2)) +
  geom_point(data = subset(check_outliers, point_group == 'outliers'), 
             colour = 'red',
             position = position_jitter(width= 0.2)) +
  labs(x = "Trio",
       y = "Mendelian errors (% of total variants)") +
  ylim(0,20) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
      axis.text = element_text(size = 14)) #+
  # geom_hline(yintercept=0.15)





# To tidy data:
# step_trios_run <- read_delim("deepvariant.r1_all_step_trios.txt", delim = " ")
# 
# me_step_results <- read_delim("stepwise_all_fmendel.txt", delim = " ", col_names = FALSE) %>% 
#   rename(FID = X1) %>% 
#   mutate(FID = as.numeric(FID)) %>% 
#   distinct(FID,.keep_all = TRUE)
# 
# step_trio_res <- left_join(me_step_results, step_trios_run, by = c(FID = "family_ID")) %>% 
#   mutate(m_err_rate = (X5/190020)*100)  %>% 
#   select(-X2, -X3) %>% 
#   rename(num_offspring = X4, num_mendel_err = X5)
# 
# write_delim(step_trio_res, "tidy_stepwise_trio_results.tsv", delim = "\t")
# 
# 
# step_trio_res
