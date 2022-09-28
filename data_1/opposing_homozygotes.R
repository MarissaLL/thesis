library(tidyverse)


# opp_hom_new <- read_delim("~/Documents/phd_python/opposing_homozygotes_newdata.csv", delim = ",", col_names = TRUE) %>% 
#   dplyr::select(indiv1, indiv2, opp_hom)


#############
# FUNCTIONS #
#############

prep_oh_data <- function(data, rship_file){
  opp_hom <- read_delim(paste0("/media/drive_6tb/projects/kakapo-genomics/output/27_opposing_homozygotes/", data, ".txt"), delim = ",") %>% 
  select(-X1)

  switched <-  opp_hom %>%
    rename(indiv1 = indiv2, indiv2 = indiv1)

  all_oh <- bind_rows(opp_hom, switched)

  big <- full_join(all_oh, rship_file, by = c(indiv1 = "indiv1", indiv2 = "indiv2")) %>% 
    mutate(status = case_when(is.na(rship) == TRUE ~ "unrel",
                              rship == "full_sibs" ~"full_sibs",
                              rship == "half_sibs" ~ "half_sibs",
                              TRUE ~ "parent_offspring" )) %>% 
    filter(is.na(opp_hom) == FALSE)
  return(big)
}


########
# MAIN #
########

pairwise_rships <- read_delim("~/projects/kakapo-genomics/data/pairwise_relationships_wo_founder.txt", delim = "\t") %>% 
  # Change around Egilsay and Elliott parents because the sample swap isn't corrected in this SNP set yet
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

new_oh <- prep_oh_data("opposing_homozygotes_new_autosome", pairwise_rships)


opp_plot <- ggplot(new_oh, aes(x = indiv1, y = opp_hom)) +
  geom_point(data = subset(new_oh, status == "unrel"), colour = "#dbcbbaff")+ 
  geom_point(data = subset(new_oh, status == "half_sibs"), colour = "#b69a96ff")+ 
  geom_point(data = subset(new_oh, status == "full_sibs"), colour = "#ca766aff")+
  geom_point(data = subset(new_oh, status == "parent_offspring"), colour = "#452736")+
  labs(x = "Individual",  y = "Number of opposing homozygotes") +
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

opp_plot



new_oh$status <- factor(new_oh$status, levels = c("unrel", "half_sibs", "full_sibs", "parent_offspring"), labels = c("Unrelated", "Half-siblings", "Full-siblings", "Parent-offspring"))

ggplot(new_oh, aes(x = status, y = opp_hom, colour = status)) +
  # geom_signif(comparisons = split(t(combn(levels(for_boxplot$status), 2)), 
  #                                 seq(nrow(t(combn(levels(for_boxplot$status), 2))))),
  #             # list(c("Unrelated", "Half-siblings"), c("Half-siblings", "Full-siblings"), c("Full-siblings", "Parent-offspring")), 
  #             map_signif_level = TRUE, tip_length = 0, textsize=6, test = "wilcox.test",
  #             step_increase = 0.07)+
  geom_point(position = position_jitter(width= 0.36))+
  geom_boxplot(outlier.size = -1, colour = "black", alpha = 0)+
  scale_colour_manual(values = c("#dbcbbaff","#b69a96ff","#ca766aff","#452736"))+
  labs(x = "Relationship",  y = "Number of opposite homozygotes") +
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)),
        axis.ticks.x = element_blank())

# Original dataset
dv_oh <- prep_oh_data("opposing_homozygotes_original_autosome", pairwise_rships)

dv_oh$status <- factor(dv_oh$status, levels = c("unrel", "half_sibs", "full_sibs", "parent_offspring"), labels = c("Unrelated", "Half-siblings", "Full-siblings", "Parent-offspring"))

ggplot(dv_oh, aes(x = status, y = opp_hom, colour = status)) +
  geom_point(position = position_jitter(width= 0.36))+
  geom_boxplot(outlier.size = -1, colour = "black", alpha = 0)+
  scale_colour_manual(values = c("#dbcbbaff","#b69a96ff","#ca766aff","#452736"))+
  labs(x = "Relationship",  y = "Number of opposite homozygotes") +
  theme_classic()+ 
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)),
        axis.ticks.x = element_blank())


#### Check significance of the differences between relatedness groups
fit <- aov(opp_hom ~ status, data=new_oh)
plot(fit)
summary(fit)
TukeyHSD(fit)


### Look at any unrelated pairs within the range of oh for parent-offspring. Founders?

oh_stats <- new_oh  %>% 
  filter(status == "Parent-offspring") %>% 
  summarise(mean = mean(opp_hom), sd = sd(opp_hom)) 

oh_unexp <- new_oh %>% 
  filter(status == "Unrelated" & opp_hom < (oh_stats$mean + oh_stats$sd))

## Old unused
# test <- new_oh %>% 
#   filter(opp_hom < max_oh) %>% 
#   filter(status == "unrel") %>% 
#   select(-status, -rship) %>% 
#   arrange(opp_hom) 
# test_output <- test %>% 
#   mutate(sel = rep_along(along = opp_hom, x = c(1,2))) %>% 
#   filter(sel == 1) %>% 
#   select(-sel)
# 
# print(xtable(test_output, digits = 0), include.rownames = FALSE)
# 
# 
# test %>% 
#   filter(indiv1 %in% founders == F) %>% 
#   filter(indiv2 %in% founders == F)
# 

### Test if lowest in group
# to_test <- test %>% 
#   select(indiv1) %>% 
#   pull() %>% 
#   unique()
# 
# best <- function(name, test_subset, main){
#   main %>% 
#     filter(indiv1 == name) %>% 
#     arrange(opp_hom) %>% 
#     left_join(test_subset, by = c(indiv1 = "indiv1", indiv2 = "indiv2"))
# }

# lapply(to_test, best, test_subset = test, main = big)

### To make a tikz file for latex
# require( tikzDevice )
# tikz( '~/Documents/phd_docs/figs/opp_hom_all_rships.tex', width = 7, height = 4 )
# print(opp_plot)
# dev.off()







####################################
### With older relationship file ###
####################################

# rships <- read_delim("~/projects/kakapo-genomics/data/bird_info_rships.csv", delim = ",") %>% 
#   dplyr::select(ID, mother, father) %>% 
#   pivot_longer(cols = -ID, names_to = "type", values_to = "parent")
# 
# 
# rswitched <-  rships %>% 
#   rename(ID = parent, parent = ID)
# 
# r_all = bind_rows(rships, rswitched)


# ggplot(test, aes(x = X1, y = X6)) +
#   geom_point()
# 

# 
# blergh <- full_join(all_test, r_all, by = c(indiv1 = "ID", indiv2 = "parent")) %>% 
#   mutate(status = case_when(is.na(type) == TRUE ~ "unrel",
#                             TRUE ~ "rel")) %>% 
#   filter(is.na(opp_hom) == FALSE)
# 
# 
# ggplot(blergh, aes(x = indiv1, y = opp_hom, colour = status)) +
#   geom_point()+ 
#   labs(x = "Individual",  y = "Number of opposing homozygotes") +
#   theme_classic()+ 
#   theme(legend.position = "none",
#         axis.title = element_text(size = 18),
#         panel.background = element_rect(fill = NA, 
#                                         colour = "black", 
#                                         size = 0.7 ),
#         axis.text = element_text(size = 14), 
#         axis.text.x = element_blank(), 
#         axis.ticks.x = element_blank())


# ### Zoomed to the bottom
# ggplot(blergh, aes(x = indiv1, y = opp_hom, colour = status)) +
#   geom_point()+ 
#   labs(x = "Individual",  y = "Number of opposing homozygotes") +
#   ylim(0, 50000) +
#   theme_classic()+ 
#   theme(legend.position = "none",
#         axis.title = element_text(size = 18),
#         panel.background = element_rect(fill = NA, 
#                                         colour = "black", 
#                                         size = 0.7 ),
#         axis.text = element_text(size = 14), 
#         axis.text.x = element_blank(), 
#         axis.ticks.x = element_blank(),
#         axis.title.y = element_blank())