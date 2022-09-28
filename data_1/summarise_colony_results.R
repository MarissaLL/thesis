library(patchwork)
library(tidyverse)

#############
# FUNCTIONS #
#############

assess_colony_results <- function(n_snps, replicate, dataset){
  colony_res <- read_delim(paste0(colony_dir, dataset, "/colony_out_sample", n_snps, "_rep", replicate, ".BestConfig"), delim = " ", trim_ws = TRUE) %>% 
    select(-ClusterIndex) %>% 
    pivot_longer(cols = -OffspringID, names_to = "Parent", values_to = "Parent_name")
  
  comparison <- full_join(colony_res, real_ped, by = c(OffspringID = "ID", "Parent")) %>% 
    mutate(result = case_when(Parent_name.x == Parent_name.y ~ "Correct",
                              (str_detect(Parent_name.x, "#")| str_detect(Parent_name.x, "\\*")| is.na(Parent_name.x)) == TRUE & is.na(Parent_name.y) == TRUE ~ "Both unassigned",
                              (str_detect(Parent_name.x, "#")| str_detect(Parent_name.x, "\\*") | is.na(Parent_name.x)) == TRUE & is.na(Parent_name.y) == FALSE ~ "Unassigned",
                              (str_detect(Parent_name.x, "#")| str_detect(Parent_name.x, "\\*") | is.na(Parent_name.x)) == FALSE & is.na(Parent_name.y) == TRUE ~ "False assignment",
                              (str_detect(Parent_name.x, "#")| str_detect(Parent_name.x, "\\*") | is.na(Parent_name.x)) == FALSE & is.na(Parent_name.y) == FALSE ~ "Incorrect",
                              TRUE ~ "unexpected_case")) %>% 
    mutate(n_snps = n_snps, 
           replicate = replicate,
           dataset = dataset)
  return(comparison)
}

###########
# GLOBALS #
###########
 seq_birds <- read_lines("/media/drive_6tb/projects/kakapo-genomics/data/all_seq_birdnames.txt")
 
 colony_dir <- "/media/drive_6tb/projects/kakapo-genomics/output/39_colony_retry/"

########
# MAIN #
########
# Because Egilsay and Elliot are misnamed in the dataset, swap their parents around in the ped so this doesn't appear as an error in the sequoia stats.
real_ped <- read_delim("/media/drive_6tb/projects/kakapo-genomics/data/bird_info_rships_nov20_wo_founderrships.csv", delim = ",") %>% 
  mutate(mother = case_when(ID == "Elliott" ~ "Suzanne",
                            ID == "Egilsay" ~ "Sue",
                            TRUE ~ mother),
         father = case_when(ID == "Elliott" ~ "Arab",
                            ID == "Egilsay" ~ "Blades",
                            TRUE ~ father)) %>% 
  select(ID, mother, father) %>% 
  dplyr::rename(MotherID = mother, FatherID = father) %>%
  pivot_longer(cols = -ID, names_to = "Parent", values_to = "Parent_name") 





colony_res <- mapply(assess_colony_results, rep(seq(100,1000,100),each = 10), rep(c(1:10), times = 10), rep(c("gwas_set", "deepvariant"), each = 100), SIMPLIFY = FALSE) %>% 
  reduce(bind_rows) #%>% 
  # group_by(n_snps, replicate, dataset, result) %>%
  # summarise(n = n()) %>%
  # ungroup() 

################# PLOT AS SCATTER (OLD) ##########

# colony_res_dots <- colony_res %>% 
#   pivot_wider(id_cols = c(n_snps, replicate, dataset), names_from = result, values_from = n , values_fill = 0) %>%
#   mutate(perc_assigned = (Correct+Incorrect)/(Correct+Incorrect+Unassigned)*100) %>%
#   mutate(correctness = (Correct)/(Correct+Incorrect)*100)
# 
# a <- ggplot(colony_res_dots, aes(x = n_snps, y = perc_assigned, colour = dataset))+
#   geom_point(position = position_jitter(10), size = 3)+
#   scale_colour_manual(values =  c( "#ca766aff","#452736", "red"),
#                       labels = c( "New", "Original"),
#                       breaks = c( "gwas_set", "deepvariant"),
#                       name = "Dataset") + # "#A9E4EF","#dbcbbaff","#b69a96ff",
#   labs(x = "SNPs sampled", y = "Percent of true parents assigned") +
#   theme_classic()+
#   theme(legend.position = "none",
#         axis.title = element_text(size = rel(1.6)),
#         panel.background = element_rect(fill = NA,
#                                         colour = "black",
#                                         size = 0.7 ),
#         axis.text = element_text(size = rel(1.2)))
# 
# 
# b <- ggplot(colony_res_dots, aes(x = n_snps, y = correctness, colour = dataset))+
#   geom_point(position = position_jitter(25), size = 3)+
#   scale_colour_manual(values =  c( "#ca766aff","#452736", "red"),
  #                     labels = c( "New", "Original"),
  #                     breaks = c( "gwas_set", "deepvariant"),
  #                     name = "Dataset") + # "#A9E4EF","#dbcbbaff","#b69a96ff",
  # labs(x = "SNPs sampled", y = "Percent of assignments that were correct") +
  # ylim(95,100)+
  # theme_classic()+
  # theme(legend.position = "right",
  #       axis.title = element_text(size = rel(1.6)),
  #       panel.background = element_rect(fill = NA,
  #                                       colour = "black",
  #                                       size = 0.7 ),
  #       axis.text = element_text(size = rel(1.2)))

# a + b

########## PLOT AS STACKED BARS ############

colony_res_stacked <- colony_res %>% 
  mutate(plot = case_when(result == "Correct" | result == "Unassigned" | result == "Incorrect" ~ "parent",
                          result == "Both unassigned" | result == "False assignment" ~ "no_parent")) %>% 
  #unite(col = combined_x, n_snps, replicate, sep = "_")
  group_by(dataset,result, n_snps,plot) %>% 
  summarise(n = n())

dv_parent <- colony_res_stacked %>% 
  filter(dataset == "deepvariant" & result %in% c("Correct", "Unassigned", "Incorrect"))

gwas_parent <- colony_res_stacked %>% 
  filter(dataset == "gwas_set" & result %in% c("Correct", "Unassigned", "Incorrect"))


dv_noparent <- colony_res_stacked %>% 
  filter(dataset == "deepvariant" & result %in% c("Both unassigned", "False assignment"))

gwas_noparent <- colony_res_stacked %>% 
  filter(dataset == "gwas_set" & result %in% c("Both unassigned", "False assignment"))

ggplot(dv_parent, aes(x = n_snps, y = n, fill = result))+
         geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values =  c("#452736")) +
  scale_x_continuous(breaks = seq(100, 1000, by = 100))+
  labs(x = "Number of SNPs",  y = "Proportion of assignments", fill = "Result") +
  theme_classic()+ 
  theme(legend.position = "right",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.6)))

ggplot(gwas_parent, aes(x = n_snps, y = n, fill = result))+
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values =  c("#dbcbbaff", "#ca766aff","#452736"))+
  scale_x_continuous(breaks = seq(100, 1000, by = 100))+
  labs(x = "Number of SNPs",  y = "Proportion of assignments", fill = "Result") +
  theme_classic()+ 
  theme(legend.position = "right",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.6)))

ggplot(dv_noparent, aes(x = n_snps, y = n, fill = result))+
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values =  c("#dbcbbaff", "#ca766aff","#452736"))+
  scale_x_continuous(breaks = seq(100, 1000, by = 100))+  
  labs(x = "Number of SNPs",  y = "Proportion of assignments", fill = "Result") +
  theme_classic()+ 
  theme(legend.position = "right",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.6)))

ggplot(gwas_noparent, aes(x = n_snps, y = n, fill = result))+
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values =  c("#dbcbbaff", "#ca766aff","#452736"))+
  scale_x_continuous(breaks = seq(100, 1000, by = 100))+
  labs(x = "Number of SNPs",  y = "Proportion of assignments", fill = "Result") +
  theme_classic()+ 
  theme(legend.position = "right",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.6)))


# Look for unexpected assigned relationships in the best run
colony_unexp <- colony_res %>% 
  filter(dataset == "gwas_set" & n_snps == 700 & result == "False assignment") %>% 
  group_by(OffspringID, Parent_name.x, Parent_name.y) %>% 
  summarise(n = n())
  