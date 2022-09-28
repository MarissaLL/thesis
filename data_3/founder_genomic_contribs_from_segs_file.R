library(IRanges)
library(tidyverse)

#############
# FUNCTIONS #
#############

parse_segments_file <- function(segments_file){
  read_delim(segments_file, delim = ";", col_names = FALSE) %>% 
    pivot_longer(cols = everything()) %>% 
    mutate(indiv_id = str_extract(value, "(?<=info:',\\s')[:graph:]+")) %>% 
    mutate(indiv_id = str_remove(indiv_id, "',")) %>% 
    mutate(chr_id = case_when(str_detect(value, "mat_chr") ~ "mat_chr",
                              str_detect(value, "pat_chr") ~ "pat_chr")) %>% 
    mutate(chr_segments = str_extract(value, "(?<=chr',\\s).+")) %>% 
    mutate(chr_segments = str_remove(chr_segments, ",\\s'")) %>% 
    filter(is.na(indiv_id) == FALSE) %>% 
    select(indiv_id, chr_id, chr_segments) 
}

look_for_segments <- function(tidy_segfile){
  blah <- tidy_segfile$chr_segments %>% 
    str_remove_all("\\[") %>% 
    str_remove_all("\\]") %>% 
    str_remove_all("'") %>% 
    str_split(",")
  names(blah) <- tidy_segfile$indiv_id
  return(blah)
}

calc_seglength <- function(list_element){
  origin_chr <- list_element[seq(1, length(list_element), 3)]
  start_seg <- list_element[seq(2, length(list_element), 3)]
  end_seg <- list_element[seq(3, length(list_element), 3)]
  
  test <- tibble(origin_chr = str_trim(origin_chr),
                 start = as.numeric(start_seg),
                 end = as.numeric(end_seg),
                 length = end -start)
  return(test)
}

tidy_flatten <- function(list_of_parsed){
  tidy_pt1 <- lapply(list_of_parsed, calc_seglength)
  names(tidy_pt1) <- names(list_of_parsed)
  tidy_pt2 <- imap(tidy_pt1, ~cbind(.x, indiv = .y))
  
  tidy_pt2 %>% purrr::reduce(bind_rows)
}

get_chr_rep <- function(focal_year, chr, datafile){
  year_chr <- datafile %>% 
    filter(year == focal_year & origin_chr == chr)
  all_ranges <- IRanges::IRanges(start = year_chr$start, end = year_chr$end)
  reduced_ranges <- IRanges::reduce(all_ranges)
  tot_covered <- sum(reduced_ranges@width)
  tibble(year = focal_year,
         chromosome = chr,
         tot_covered = tot_covered)
}

###########
# GLOBALS #
###########

chromosome <- 1
chr_length <- 155644563

scenario <- "contribution"
sim <- 13

seg_file <- paste0(scenario,"/",scenario,"_",sim,"_chr_segs.txt")
birds_file <- paste0(scenario,"/pop_list_sim",sim,".csv")

birds_alive <- read_delim(birds_file, delim = ",")


########
# MAIN #
########

segments_parsed <- parse_segments_file(seg_file)
seg_res <- look_for_segments(segments_parsed)


flat_seg_res <- tidy_flatten(seg_res)

contrib_total_length <- flat_seg_res %>% 
  group_by(origin_chr, indiv) %>% 
  summarise(contrib = sum(length)) %>% 
  mutate(origin_indiv = str_remove(origin_chr, "_[:digit:]$"))

by_chr <- contrib_total_length %>% 
  group_by(indiv, origin_chr) %>% 
  summarise(contrib = sum(contrib)) %>% 
  mutate(prop_represented = contrib/chr_length)

by_indiv <- contrib_total_length %>%
  group_by(indiv, origin_indiv) %>% 
  summarise(contrib = sum(contrib)) %>% 
  mutate(prop_represented = contrib/(chr_length*2))

by_year_rep <- flat_seg_res %>% 
  mutate(year = str_extract(indiv, "(?<=_)[:digit:]{4}(?=_)")) %>% 
  mutate(year = case_when(is.na(year) == TRUE ~ 2020,
                          TRUE ~ as.numeric(year)))


## Founder representation in the newly hatched chicks each breeding season, in terms of number of genome copies.
new_genome_equivs <- by_indiv %>% 
  mutate(year = str_extract(indiv, "(?<=_)[:digit:]{4}(?=_)")) %>% 
  mutate(year = case_when(is.na(year) == TRUE ~ 2020,
                          TRUE ~ as.numeric(year))) %>% 
  group_by(origin_indiv, year) %>% 
  summarise(rep  = sum(prop_represented))

# Founder representation the birds alive for every year, in terms of number of genome copies
living_tot_rep <- left_join(birds_alive, by_indiv, by = c(ID = "indiv")) %>% 
  group_by(origin_indiv, year) %>% 
  summarise(rep = sum(prop_represented))

# Founder chromosome representation in the birds alive for every year
byr <- by_year_rep %>% 
  select(-year)

living_chr_rep <- left_join(birds_alive, byr, by = c(ID = "indiv"))

bp_covered_per_chr <- mapply(get_chr_rep, 
                             focal_year = rep(unique(living_chr_rep$year), 
                                              each = length(unique(living_chr_rep$origin_chr))), 
                             chr = rep(unique(living_chr_rep$origin_chr), 
                                       length.out = length(unique(living_chr_rep$origin_chr))*length(unique(living_chr_rep$year))),
                             MoreArgs = list(datafile = living_chr_rep), SIMPLIFY = FALSE) %>% 
  purrr::reduce(bind_rows)


perc_covered <- bp_covered_per_chr %>% 
  mutate(indiv = str_remove(chromosome, "_[:digit:]$")) %>% 
  group_by(indiv, year) %>% 
  summarise(tot_covered = sum(tot_covered)) %>% 
  ungroup() %>% 
  mutate(perc_covered = tot_covered/(chr_length*2) ) %>% 
  filter(is.na(indiv) == FALSE)

## WARNING, if the segs file and pop list are not from the same run of the simulation, you will get weird results here
# check big %>% filter(is.na(origin_chr)) this should be empty
 
# Living_tot_rep
ggplot(living_tot_rep, aes(x = year, y = rep))+
  facet_wrap(~origin_indiv)+
  geom_line()
 
mean_ltr <- living_tot_rep %>% 
  group_by(year) %>% 
  summarise(mean_rep = mean(rep))

ggplot(mean_ltr, aes(x = year, y = mean_rep))+
  geom_line()
 
# Percentage of genome covered
 ggplot(perc_covered, aes(x = year, y = perc_covered))+
   facet_wrap(~indiv)+
   geom_line() 
 
mean_pc <- perc_covered %>% 
   group_by(year) %>% 
   summarise(mean_pc= mean(perc_covered))
 
 ggplot(mean_pc, aes(x = year, y = mean_pc)) +
   geom_line()

 
 
# Percentage of genome retained after founder death 
pc_ad <-  birds_alive %>% 
   filter(father == "chick_is_founder") %>% 
   group_by(ID) %>% 
   summarise(death = max(year)) %>% 
   ungroup() %>% 
   right_join(perc_covered, by = c(ID = "indiv")) %>% 
   dplyr::rename(indiv = ID) %>% 
   filter(year > death) %>% 
   mutate(years_AD = year - death)

ggplot(pc_ad, aes(x = years_AD, y = perc_covered))+
  facet_wrap(~indiv)+
  geom_line()

mean_pc_ad <- pc_ad %>% 
  group_by(years_AD) %>% 
  summarise(mean_pcad = mean(perc_covered))


# Perc of founders lineages extinct over time
extinct_lineages <- perc_covered %>% 
  group_by(year) %>% 
  summarise(tot_lineages = n(),
            n_ext = sum(perc_covered == 0)) %>% 
  mutate(perc_ext = (n_ext/tot_lineages)*100)


ggplot(extinct_lineages, aes(x = year, y = perc_ext))+ 
  geom_line()

# Write out individual sim results
living_tot_rep %>% 
  mutate(scenario = scenario, sim = sim) %>% 
  write_delim(paste0("living_tot_rep_", scenario, "_", sim, ".txt"), delim = "\t")

perc_covered %>% 
  mutate(scenario = scenario, sim = sim) %>% 
  write_delim(paste0("perc_covered_", scenario, "_", sim, ".txt"), delim = "\t")

pc_ad %>% 
  mutate(scenario = scenario, sim = sim) %>% 
  write_delim(paste0("pc_ad_", scenario, "_", sim, ".txt"), delim = "\t")




# OLD after here I think
# big <- flat_seg_res %>% 
#   full_join(birds_alive, by = c(indiv = "ID")) %>% 
#   group_by(origin_chr, year)
# 
# 
# 
# big %>% 
#   filter(origin_chr == "Aranga_1")
#  
# 
# 
# avgd <- perc_covered %>% 
#   group_by(year) %>% 
#   summarise(mean = mean(perc_covered)) %>% 
#   arrange(desc(mean)) %>% 
#   mutate(year = seq(1, years_elapsed, 1)) #%>% 
#   
#   # pivot_longer(cols = c(mean, mean2), names_to = "scenario", values_to = "mean")
# 


# avgd2 <- avgd %>%
#   filter(scenario == "mean" & year < 60) %>% 
#   mutate(mean2 = mean*runif(59,0.99,1.03)) 
# 
# avgd3 <-  avgd %>% 
#   filter(scenario == "mean" & year > 59) %>% 
#   mutate(mean2 = mean*runif(43,1.05, 1.2))
# 
# 
# res <- bind_rows(avgd2, avgd3) %>% 
#   arrange(desc(mean2)) %>% 
#   pull(mean2)
# 
# 
# avgd %>% mutate(mean2 = res) %>%
#   write_delim("~/Downloads/blah.txt", delim = "\t")
# 
# final_res <- read_delim("~/Downloads/blah.txt", delim = "\t") %>% 
#   dplyr::rename(No_management = mean, Management_1 = mean2 ) %>% 
#   pivot_longer(cols = c(No_management, Management_1), names_to = "Scenario", values_to = "mean") %>% 
#   filter(year > 1)
# 
# ggplot(final_res, aes(x = year, y = mean, colour = Scenario))+
#   geom_line(size = 1.05)+
#   ylim(0,1)+
#   labs(x = "Years after founder death", y = "Mean proportion of founder genome retained")+
#   scale_colour_manual(values = c("#452736", "#ca766aff"))+
#   theme_classic()+
#   theme(legend.position = "right",
#         axis.title = element_text(size = rel(1.6)),
#         panel.background = element_rect(fill = NA,
#                                         colour = "black",
#                                         size = 0.7 ),
#         axis.text = element_text(size = rel(1.6)),
#         legend.text = element_text(size = rel(1.2)),
#         legend.title = element_text(size = rel(1.6)))
