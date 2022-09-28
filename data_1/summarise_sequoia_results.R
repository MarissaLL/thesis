library(patchwork)
library(tidyverse)

#############
# FUNCTIONS #
#############

import_sequoia_results <- function(n_snps, replicate, dataset){
  sequoia_res <- read_delim(paste0(sequoia_dir, "sequoia_sample_", 
                                 n_snps, "_rep", replicate, "_", 
                                 dataset, ".txt"), delim = "\t", trim_ws = TRUE) %>% 
    select(id, dam, sire) %>% 
    rename(mother = dam, father = sire) %>% 
    pivot_longer(-id, names_to = "Parent", values_to = "Parent_name")
  
  comparison <- left_join(sequoia_res, real_ped, by = c(id = "ID", Parent = "Parent")) %>% 
    mutate(result = case_when(Parent_name.x == Parent_name.y ~ "Correct",
                              (str_detect(Parent_name.x, "M0")| str_detect(Parent_name.x, "F0")| is.na(Parent_name.x)) == TRUE & is.na(Parent_name.y) == TRUE ~ "Both unassigned",
                              (str_detect(Parent_name.x, "M0")| str_detect(Parent_name.x, "F0") | is.na(Parent_name.x)) == TRUE & is.na(Parent_name.y) == FALSE ~ "Unassigned",
                              (str_detect(Parent_name.x, "M0")| str_detect(Parent_name.x, "F0") | is.na(Parent_name.x)) == FALSE & is.na(Parent_name.y) == TRUE ~ "False assignment",
                              (str_detect(Parent_name.x, "M0")| str_detect(Parent_name.x, "F0") | is.na(Parent_name.x)) == FALSE & is.na(Parent_name.y) == FALSE ~ "Incorrect",
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
 
 sequoia_dir <- "/media/drive_6tb/projects/kakapo-genomics/output/37_sequoia_pt2/"

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
  pivot_longer(-ID, names_to = "Parent", values_to = "Parent_name") 



sequoia_res <- mapply(import_sequoia_results, rep(seq(100,1000,100),each = 10), rep(c(1:10), times = 10), rep(c("gwas_set", "deepvariant"), each = 100), SIMPLIFY = FALSE) %>% 
  reduce(bind_rows)

###### Plot with stacked bars #######


sequoia_res_stacked <- sequoia_res %>% 
  mutate(plot = case_when(result == "Correct" | result == "Unassigned" | result == "Incorrect" ~ "parent",
                          result == "Both unassigned" | result == "False assignment" ~ "no_parent"))%>% 
  group_by(dataset,result, n_snps,plot) %>% 
  summarise(n = n())

dv_parent <- sequoia_res_stacked %>% 
  filter(dataset == "deepvariant" & result %in% c("Correct", "Unassigned", "Incorrect"))

gwas_parent <- sequoia_res_stacked %>% 
  filter(dataset == "gwas_set" & result %in% c("Correct", "Unassigned", "Incorrect"))


dv_noparent <- sequoia_res_stacked %>% 
  filter(dataset == "deepvariant" & result %in% c("Both unassigned", "False assignment"))

gwas_noparent <- sequoia_res_stacked %>% 
  filter(dataset == "gwas_set" & result %in% c("Both unassigned", "False assignment"))


ggplot(dv_parent, aes(x = n_snps, y = n, fill = result))+
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values =  c("#ca766aff","#452736")) +
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

# Find potential missed rships
sequoia_unexp <- sequoia_res %>% 
  filter(n_snps ==200 & 
           dataset == "gwas_set" & 
           result == "False assignment") %>% 
  group_by(id,Parent_name.x) %>% 
  summarise(n = n()) %>% 
  filter(n > 4) # Has to be in at least 50% of replicates
