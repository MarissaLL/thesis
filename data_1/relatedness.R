library(tidyverse)
library(patchwork)

#############
# FUNCTIONS #
#############

organise_king_for_plot <- function(version, rships_file){
  snps <- read_delim(paste0("/media/drive_6tb/projects/kakapo-genomics/output/28_relatedness/", version, ".kin0"), delim = "\t") %>% 
    rename(IID1 = `#IID1`) %>% 
    left_join(rships_file, by = c(IID1 = "indiv1", IID2 = "indiv2")) %>% 
    mutate(rship2 = case_when(is.na(rship) == TRUE ~ "unrel",
                              rship == "offspring_father" ~ "parent_offspring",
                              rship == "offspring_mother" ~ "parent_offspring",
                              rship == "full_sibs" ~ "full_sibs",
                              rship == "half_sibs" ~ "half_sibs",
                              TRUE ~ "ERROR"))
  snps$rship2 <- factor(snps$rship2, levels = c("unrel", "half_sibs", "full_sibs", "parent_offspring"), labels = c("Unrelated", "Half-siblings", "Full-siblings", "Parent-offspring"))
  return(snps)
}

organise_rxy_for_plot <- function(dataset, rships_file){
  
  rel <- read_delim(paste0("output/28_relatedness/", dataset, ".rxy"), delim = "\t") %>% 
    select(a,b,rab)
  
  names <- read_lines(paste0("output/28_relatedness/",dataset, ".samplenames")) %>% 
    str_split(" ") %>% 
    unlist() %>% 
    unique()
  
  res <- tibble(name1 = as_factor(names),
                val = seq(0,length(names)-1)) 
  
  rxy_vals <- rel %>% 
    left_join(res, by = c(a = "val")) %>% 
    left_join(res, by = c(b = "val")) %>% 
    filter(is.na(rab) == FALSE) %>% 
    select(-a, -b)
  
  rxy_res <- rxy_vals %>% 
  left_join(rships_file, by = c(name1.x = "indiv1", name1.y = "indiv2")) %>% 
    mutate(rship2 = case_when(is.na(rship) == TRUE ~ "unrel",
                              rship == "offspring_father" ~ "parent_offspring",
                              rship == "offspring_mother" ~ "parent_offspring",
                              rship == "full_sibs" ~ "full_sibs",
                              rship == "half_sibs" ~ "half_sibs",
                              TRUE ~ "ERROR"))
  rxy_res$rship2 <- factor(rxy_res$rship2, levels = c("unrel", "half_sibs", "full_sibs", "parent_offspring"), labels = c("Unrelated", "Half-siblings", "Full-siblings", "Parent-offspring"))
  return(rxy_res)
}

########
# MAIN #
########

pairwise_rships <- read_delim("/media/drive_6tb/projects/kakapo-genomics/data/pairwise_relationships_wo_founder.txt", delim = '\t') %>% 
  # Change around Egilsay and Elliott parents because they aren't corrected in this SNP set yet
  mutate(indiv1 = case_when(indiv1 == "Egilsay" ~ "E1",
                            indiv1 == "Elliott" ~ "E2",
                            TRUE ~ indiv1),
         indiv2 = case_when(indiv2 == "Egilsay" ~ "E1",
                            indiv2 == "Elliott" ~ "E2",
                            TRUE ~ indiv2)) %>% 
  mutate(indiv1 = case_when(indiv1 == "E1" ~ "Elliott",
                            indiv1 == "E2" ~ "Egilsay",
                            TRUE ~ indiv1),
         indiv2 = case_when(indiv2 == "E1" ~ "Elliott",
                            indiv2 == "E2" ~ "Egilsay",
                            TRUE ~ indiv2))
  

# KING relatedness
new_snps <- organise_king_for_plot("new", pairwise_rships) 
original_snps <- organise_king_for_plot("original", pairwise_rships)


a <- ggplot(new_snps, aes(x = rship2, y = KINSHIP, colour = rship2))+ 
  geom_point(position = position_jitter(0.3)) +
  geom_boxplot(outlier.size = -1, colour = "black", alpha = 0)+
  scale_colour_manual(values = c("#dbcbbaff","#b69a96ff","#ca766aff","#452736"))+
  labs(x = "Relationship",  y = "Kinship") +
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)))


b <- ggplot(original_snps, aes(x = rship2, y = KINSHIP, colour = rship2)) +
  geom_point(position = position_jitter(0.3)) +
  geom_boxplot(outlier.size = -1, colour = "black", alpha = 0)+
  scale_colour_manual(values = c("#dbcbbaff","#b69a96ff","#ca766aff","#452736"))+
  labs(x = "Relationship",  y = "Kinship") +
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)))

a + b


# Rxy relatedness
new_rxy <- organise_rxy_for_plot("gwas_set_intermediate_ld", pairwise_rships)
dv_rxy <- organise_rxy_for_plot("deepvariant.r1_intermediate_ld", pairwise_rships)


c <- ggplot(new_rxy, aes(x = rship2, y = rab, colour = rship2)) +
  geom_point(position = position_jitter(0.3)) +
  geom_boxplot(outlier.size = -1, colour = "black", alpha = 0)+
  scale_colour_manual(values = c("#dbcbbaff","#b69a96ff","#ca766aff","#452736"))+
  labs(x = "Relationship",  y = "Kinship") +
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)))

d <- ggplot(dv_rxy, aes(x = rship2, y = rab, colour = rship2)) +
  geom_point(position = position_jitter(0.3)) +
  geom_boxplot(outlier.size = -1, colour = "black", alpha = 0)+
  scale_colour_manual(values = c("#dbcbbaff","#b69a96ff","#ca766aff","#452736"))+
  labs(x = "Relationship",  y = "Kinship") +
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)))

c + d

# Get unexpectedly high relatedness
king_mean_sd <- new_snps %>% 
  filter(rship2 == "Parent-offspring") %>% 
  summarise(mean = mean(KINSHIP), sd = sd(KINSHIP))

king_unexp <- new_snps %>% 
  filter(rship2 == "Unrelated" & KINSHIP > (king_mean_sd$mean - king_mean_sd$sd))

rxy_mean_sd <- new_rxy %>% 
  filter(rship2 == "Parent-offspring") %>% 
  summarise(mean = mean(rab), sd = sd(rab))

rxy_unexp <- new_rxy %>% 
  filter(rship2 == "Unrelated" & rab > (rxy_mean_sd$mean - rxy_mean_sd$sd))
