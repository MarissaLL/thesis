library(tidyverse)

get_cpt_files <- function(chr, ncycles){
    read_delim(paste0("/media/drive_6tb/projects/kakapo-genomics/output/25_tidy_alphapeel/changepoint_files/changepoints_", chr, "_ncycles", ncycles, ".txt"), delim = "\t") %>%
    select(indiv, changepoints, focal_parent) %>%
    rename(!!paste0("changepoints_", chr) := changepoints)
}



chrs <- c("S1","S2","S4", "S5", "S6", "S7", "S8", "S9", "S10", 
          "S11", "S12",  "S14", "S15", "S16", "S17", "S18", "S19", 
          "S20", "S21", "S22", "S23", "S24", "S25", "S26")

### n cpts data
ncycles20 <- lapply(chrs, get_cpt_files, ncycles = 20) %>%
  reduce(full_join, by = c("indiv", "focal_parent")) %>%
  pivot_longer(cols = c(-indiv, -focal_parent),
               names_to = "chromosome", values_to = "n_cpts20") %>% 
  mutate(n_cpts20 = as.numeric(n_cpts20))

ncycles100 <- lapply(chrs, get_cpt_files, ncycles = 100) %>%
  reduce(full_join, by = c("indiv", "focal_parent")) %>%
  pivot_longer(cols = c(-indiv, -focal_parent),
               names_to = "chromosome", values_to = "n_cpts100") %>% 
  mutate(n_cpts100 = as.numeric(n_cpts100))

ncycles200 <- lapply(chrs, get_cpt_files, ncycles = 200) %>%
  reduce(full_join, by = c("indiv", "focal_parent")) %>%
  pivot_longer(cols = c(-indiv, -focal_parent),
             names_to = "chromosome", values_to = "n_cpts200") %>% 
  mutate(n_cpts200 = as.numeric(n_cpts200))

ncycles1000 <- lapply(chrs, get_cpt_files, ncycles = 1000) %>%
  reduce(full_join, by = c("indiv", "focal_parent")) %>%
  pivot_longer(cols = c(-indiv, -focal_parent),
               names_to = "chromosome", values_to = "n_cpts1000") %>% 
  mutate(n_cpts1000 = as.numeric(n_cpts1000))


### Not used in this script currently
# ncycles20 <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/25_tidy_alphapeel/recombination_positions/all_recombination_ncycles_20.txt", delim = ' ')
# ncycles100 <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/25_tidy_alphapeel/recombination_positions/all_recombination_ncycles_100.txt", delim = ' ')
# ncycles200 <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/25_tidy_alphapeel/recombination_positions/all_recombination_ncycles_200.txt", delim = ' ')
# ncycles1000 <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/25_tidy_alphapeel/recombination_positions/all_recombination_ncycles_1000.txt", delim = ' ')

######## Need to fix this up to include the fourth case
all_ncycles <- full_join(ncycles20, ncycles100) %>% 
  full_join(ncycles200) %>% 
  full_join(ncycles1000) %>% 
  mutate(concordant = case_when(n_cpts20 == n_cpts100 & n_cpts20 == n_cpts200 ~ "all",
                                n_cpts20 == n_cpts100 & n_cpts20 != n_cpts200 ~ "20_100",
                                n_cpts20 != n_cpts100 & n_cpts20 == n_cpts200 ~ "20_200",
                                n_cpts20 != n_cpts100 & n_cpts20 != n_cpts200 & n_cpts100 == n_cpts200 ~ "100_200",
                                n_cpts20 != n_cpts100 & n_cpts20 != n_cpts200 & n_cpts100 != n_cpts200 ~ "none",
                                TRUE ~ "poor_phasing_na"))

### filtering data
get_kept_indivs <- function(chr, ncycles, parent){
  read_delim(paste0("/media/drive_6tb/projects/kakapo-genomics/output/25_tidy_alphapeel/phasing_filter/indivs_to_keep_", chr, "_ncycles_", ncycles, "_", parent, ".txt"), delim = "\t") %>% 
  rename(!!paste0("kept_", ncycles) := chick) %>% 
  mutate(chromosome = !!paste0("changepoints_", chr))
}

kept_by_ncycles <- function(chrs, ncycles, cpts_df, parent){
ncycles_kept <- lapply(chrs, get_kept_indivs, ncycles = ncycles, parent = parent) %>% 
  reduce(bind_rows) %>% 
  mutate(kept_ncycles = "kept")

cpts_df %>% 
  filter(focal_parent == paste0(parent, "_only")) %>% 
  left_join(ncycles_kept , by = c(indiv = paste0("kept_", ncycles), chromosome = "chromosome")) %>% 
  mutate(kept_ncycles = case_when(is.na(kept_ncycles) ~ "discarded",
                          TRUE ~ kept_ncycles)) %>% 
  rename(!!paste0("kept_ncycles", ncycles) := kept_ncycles)
}



# ## This section still needs condensing into a function - done. see condense
# kept_ncycles20_mat <- kept_by_ncycles(chrs = chrs, ncycles = 20, cpts_df = ncycles20, parent = "mat") %>% 
#     mutate(n_cpts20 = case_when(kept_ncycles20 == "discarded" ~ -55,
#                               kept_ncycles20 == "kept" ~ n_cpts20))
# kept_ncycles20_pat <- kept_by_ncycles(chrs = chrs, ncycles = 20, cpts_df = ncycles20, parent = "pat") %>% 
#   mutate(n_cpts20 = case_when(kept_ncycles20 == "discarded" ~ -55,
#                               kept_ncycles20 == "kept" ~ n_cpts20))
# 
# 
# kept_ncycles100_mat <- kept_by_ncycles(chrs = chrs, ncycles = 100, cpts_df = ncycles100, parent = "mat") %>% 
#   mutate(n_cpts100 = case_when(kept_ncycles100 == "discarded" ~ -55,
#                                      kept_ncycles100 == "kept" ~ n_cpts100))
# kept_ncycles100_pat <- kept_by_ncycles(chrs = chrs, ncycles = 100, cpts_df = ncycles100, parent = "pat") %>% 
#   mutate(n_cpts100 = case_when(kept_ncycles100 == "discarded" ~ -55,
#                                kept_ncycles100 == "kept" ~ n_cpts100))
# 
# 
# kept_ncycles200_mat <- kept_by_ncycles(chrs = chrs, ncycles = 200, cpts_df = ncycles200, parent = "mat") %>% 
#   mutate(n_cpts200 = case_when(kept_ncycles200 == "discarded" ~ -55,
#                                      kept_ncycles200 == "kept" ~ n_cpts200))
# kept_ncycles200_pat <- kept_by_ncycles(chrs = chrs, ncycles = 200, cpts_df = ncycles200, parent = "pat") %>% 
#   mutate(n_cpts200 = case_when(kept_ncycles200 == "discarded" ~ -55,
#                                kept_ncycles200 == "kept" ~ n_cpts200))
# 
# kept_ncycles1000_mat <- kept_by_ncycles(chrs = chrs, ncycles = 1000, cpts_df = ncycles1000, parent = "mat") %>% 
#   mutate(n_cpts1000 = case_when(kept_ncycles1000 == "discarded" ~ -55,
#                                kept_ncycles1000 == "kept" ~ n_cpts1000))
# kept_ncycles1000_pat <- kept_by_ncycles(chrs = chrs, ncycles = 1000, cpts_df = ncycles1000, parent = "pat") %>% 
#   mutate(n_cpts1000 = case_when(kept_ncycles1000 == "discarded" ~ -55,
#                                 kept_ncycles1000 == "kept" ~ n_cpts1000))
# kept_ncycles1000 <-  bind_rows(kept_ncycles1000_mat, kept_ncycles1000_pat)


### Need to remember the right format for the rhs of this 
condense <- function(chrs, ncycles, cpts_df){
  mat <- kept_by_ncycles(chrs = chrs, ncycles = ncycles, cpts_df = cpts_df, parent = "mat") %>% 
    mutate(!!paste0("n_cpts", ncycles) := case_when(!!sym(paste0("kept_ncycles", ncycles)) == "discarded" ~ -55, 
                                                    !!sym(paste0("kept_ncycles", ncycles)) == "kept" ~ !!sym(paste0("n_cpts", ncycles))))
  pat <- kept_by_ncycles(chrs = chrs, ncycles = ncycles, cpts_df = cpts_df, parent = "pat") %>% 
    mutate(!!paste0("n_cpts", ncycles) := case_when(!!sym(paste0("kept_ncycles", ncycles)) == "discarded" ~ -55, 
                                                    !!sym(paste0("kept_ncycles", ncycles)) == "kept" ~ !!sym(paste0("n_cpts", ncycles))))
  bind_rows(mat, pat)
}

kept_ncycles20  <- condense(chrs = chrs, 20, ncycles20)
kept_ncycles100  <- condense(chrs = chrs, 100, ncycles100)
kept_ncycles200  <- condense(chrs = chrs, 200, ncycles200)
kept_ncycles1000  <- condense(chrs = chrs, 1000, ncycles1000)

all_ncycles_kept_info <- full_join(kept_ncycles20, kept_ncycles100) %>% 
  full_join(kept_ncycles200) %>% 
  full_join(kept_ncycles1000) 

all_ncycles_na <- all_ncycles_kept_info %>% 
  replace((. == -55), NA_real_) %>% 
  mutate(concordant = case_when(n_cpts20 == n_cpts100 & n_cpts20 == n_cpts200 ~ "all",
                                n_cpts20 == n_cpts100 & n_cpts20 != n_cpts200 ~ "20_100",
                                n_cpts20 != n_cpts100 & n_cpts20 == n_cpts200 ~ "20_200",
                                n_cpts20 != n_cpts100 & n_cpts20 != n_cpts200 & n_cpts100 == n_cpts200 ~ "100_200",
                                n_cpts20 != n_cpts100 & n_cpts20 != n_cpts200 & n_cpts100 != n_cpts200 ~ "none",
                                TRUE ~ "poor_phasing_na")) %>% 
  mutate(actual_na_status = case_when(concordant != "poor_phasing_na" ~ concordant,
                                      concordant == "poor_phasing_na" & n_cpts20 == n_cpts100 ~ "20_100",
                                      concordant == "poor_phasing_na" & n_cpts20 == n_cpts200 ~ "20_200",
                                      concordant == "poor_phasing_na" & n_cpts100 == n_cpts200 ~ "100_200",
                                      TRUE ~ "unreliable")) %>% 
  mutate(chr = as.numeric(str_remove(chromosome, "changepoints_S"))) 
  


##### Heatmaps

for_heatmap <- all_ncycles_kept_info %>% 
  mutate(chr = as.numeric(str_remove(chromosome, "changepoints_S"))) %>% 
  select(indiv, chr, n_cpts20, n_cpts100, n_cpts200, n_cpts1000) %>% 
  replace(is.na(.), -99) 


#### Filtered vs not
ggplot(for_heatmap, aes(x = chr, y = reorder(indiv, desc(indiv)))) +
  geom_tile(aes(fill = n_cpts20)) +
  scale_x_continuous(breaks = c(1:26))+
  labs(title = "Ncycles 20", x = "Chromosome",y = "Individual")+
  theme(legend.position = 'none')
  

ggplot(for_heatmap, aes(x = chr, y = reorder(indiv, desc(indiv)))) +
  geom_tile(aes(fill = n_cpts100)) +
  scale_x_continuous(breaks = c(1:26))+
  labs(title = "Ncycles 100", x = "Chromosome", y = "Individual")+
  theme(legend.position = 'none')

ggplot(for_heatmap, aes(x = chr, y = reorder(indiv, desc(indiv)))) +
  geom_tile(aes(fill = n_cpts200)) +
  scale_x_continuous(breaks = c(1:26))+
  labs(title = "Ncycles 200",x = "Chromosome", y = "Individual")+
  theme(legend.position = 'none')

ggplot(for_heatmap, aes(x = chr, y = reorder(indiv, desc(indiv)))) +
  geom_tile(aes(fill = n_cpts1000)) +
  scale_x_continuous(breaks = c(1:26))+
  labs(title = "Ncycles 1000", x = "Chromosome",y = "Individual")+
  theme(legend.position = 'none')


####  Indiv and chromosome scores (Just ncycles 20 for now)
library(grid)
library(gtable)

a <-  ggplot(for_heatmap, aes(x = chr, y = reorder(indiv, desc(indiv)))) +
  geom_tile(aes(fill = n_cpts1000)) +
  scale_x_continuous(breaks = c(1:26))+
  # labs(title = "Ncycles 20", y = "Individual")+
  theme_void()+
  theme(legend.position = 'none')

indiv_scores <- all_ncycles_na %>% 
  group_by(indiv) %>% 
  summarise(bad = sum(is.na(n_cpts1000))) %>% 
  mutate(overall_score = 48 - bad)

chr_scores <- all_ncycles_na %>% 
  group_by(chr) %>% 
  summarise(bad = sum(is.na(n_cpts1000))) %>% 
  mutate(overall_score = 248 - bad)

b <- ggplot(indiv_scores, aes(x = 1, y = reorder(indiv, desc(indiv)), fill = overall_score))+ 
  geom_tile()+
  scale_fill_gradient(low="darkred", high="white", limits = c(0,48))+
  theme_void()+
  theme(legend.position = 'none')


c <- ggplot(chr_scores, aes(x = chr, y = 1, fill = overall_score))+ 
  geom_tile()+
  scale_fill_gradient(low="darkred", high="white", limits = c(0,248))+
  theme_void()+
  theme(legend.position = 'none')


aa <- ggplotGrob(a)
bb <- ggplotGrob(b)
cc <- ggplotGrob(c)
blank_one <- rectGrob(gp = gpar(fill = "white",col = "white"))

g <- gtable_matrix(name = "ncycles20_plot",
                   grobs = matrix(list(aa,bb,
                                       blank_one, blank_one,
                                       cc, blank_one), byrow = TRUE, nrow =3, ncol = 2),
                   widths = unit(c(20,1), "cm"),
                   heights = unit(c(20,0.7, 1), "cm"))

grid.newpage()
grid.draw(g)



#### Recombo when unfiltered

for_na_heatmap <- for_heatmap %>% 
  replace((. == -55), NA_real_) %>% 
  replace((. == -99), NA_real_)


ggplot(for_na_heatmap, aes(x = chr, y = reorder(indiv, desc(indiv)))) +
  geom_tile(aes(fill = n_cpts20)) +
  scale_x_continuous(breaks = c(4:26))+
  scale_fill_continuous(na.value = 'grey')

ggplot(for_na_heatmap, aes(x = chr, y = reorder(indiv, desc(indiv)))) +
  geom_tile(aes(fill = n_cpts100)) +
  scale_x_continuous(breaks = c(4:26))+
  scale_fill_continuous(na.value = 'grey')

ggplot(for_na_heatmap, aes(x = chr, y = reorder(indiv, desc(indiv)))) +
  geom_tile(aes(fill = n_cpts200)) +
  scale_x_continuous(breaks = c(4:26))+
  scale_fill_continuous(na.value = 'grey')

ggplot(for_na_heatmap, aes(x = chr, y = reorder(indiv, desc(indiv)))) +
  geom_tile(aes(fill = n_cpts1000)) +
  scale_x_continuous(breaks = c(4:26))+
  scale_fill_continuous(na.value = 'grey')





#### Overall n_cpts per indiv and n_cpts per individual

missingness <- all_ncycles_na %>%
  group_by(indiv) %>%
  summarise(missingness20 = sum(is.na(n_cpts20)== TRUE),
            ncpts20 = sum(n_cpts20,na.rm = TRUE),
            missingness100 = sum(is.na(n_cpts100)== TRUE),
            ncpts100 = sum(n_cpts100,na.rm = TRUE),
            missingness200 = sum(is.na(n_cpts200)== TRUE),
            ncpts200 = sum(n_cpts200,na.rm = TRUE),
            missingness1000 = sum(is.na(n_cpts1000)== TRUE),
            ncpts1000 = sum(n_cpts1000,na.rm = TRUE)) %>% 
  pivot_longer(cols = -indiv) %>% 
  mutate(type = str_remove(name, "[:digit:]+"),
         ncycles = as.numeric(str_remove(name, "[:alpha:]+"))) %>% 
  select(-name) %>% 
  mutate(ncycles = as.factor(ncycles))



ggplot(subset(missingness, type == "missingness"), aes(x = value, fill = ncycles))+
  facet_wrap(~ ncycles)+
  labs(x = "Chromosomes excluded per individual due to poor phasing", y = "Number of individuals") +
  geom_histogram(binwidth = 1)+
  scale_fill_manual(values = c( "#CED38C", "#DCC949", "#BCA888", "#7D9D33")) +
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)))


ggplot(subset(missingness, type == "ncpts"), aes(x = value, fill = ncycles))+
  facet_wrap(~ ncycles)+
  labs(x = "Total predicted recombinations per individual", y = "Number of individuals") +
  geom_histogram(binwidth = 2)+
  scale_fill_manual(values = c( "#CED38C", "#DCC949", "#BCA888", "#7D9D33")) +
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)))



ncpts_per_chr <- missingness %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  mutate(presentness = 48 - missingness,
         ncpts_chr = ncpts/presentness)


ggplot(ncpts_per_chr, aes(x = ncpts_chr, fill = ncycles))+
  facet_wrap(~ncycles)+
  labs(x = "Mean recombinations per chromosome", y = "Number of individuals") +
  geom_histogram(binwidth = 0.25)+
  scale_fill_manual(values = c( "#CED38C", "#DCC949", "#BCA888", "#7D9D33")) +
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)))
