library(parallel)
library(tidyverse)

#############
# FUNCTIONS #
#############


make_haplos <- function(start_end){
  
  filtered_vcf <- vcf %>% 
    filter(POS > start_end & POS < (start_end + 10000)) %>% 
    select(-POS)
  
  names <- names(filtered_vcf)[names(filtered_vcf) != "ID"]
  
  tidied <- lapply(names, function(name){
    filtered_vcf %>% 
      select(ID, name) %>% 
      mutate(!!sym(paste0(name, "_1")) := case_when(!!sym(name) == "0|0" ~ "0",
                                                    !!sym(name) == "0|1" ~ "0",
                                                    !!sym(name) == "1|0" ~ "1",
                                                    !!sym(name) == "1|1" ~ "1",
                                                    TRUE ~ "N")) %>% 
      mutate(!!sym(paste0(name, "_2")) := case_when(!!sym(name) == "0|0" ~ "0",
                                                    !!sym(name) == "0|1" ~ "1",
                                                    !!sym(name) == "1|0" ~ "0",
                                                    !!sym(name) == "1|1" ~ "1",
                                                    TRUE ~ "N")) %>% 
      select(-name)}) %>% 
    reduce(full_join, by = "ID")
  
  aa <- tidied %>% 
    pivot_longer(cols = -ID, names_to = "haplo") %>% 
    pivot_wider(names_from = ID) %>% 
    unite(-haplo, col = "haplo_sequence", sep = "") %>% 
    mutate(n_count = str_count(haplo_sequence, "\\.")) %>% 
    mutate(name = str_remove(haplo, "_[:digit:]"))
  
  
  founder_haplos <- aa %>% 
    filter(name %in% founders == TRUE) %>% 
    pull(haplo_sequence) %>% 
    unique()
  
  nonfounder_haplos <- aa %>% 
    filter(name %in% founders == FALSE) %>% 
    pull(haplo_sequence) %>% 
    unique()
  
  res <- crossing(founder_haplos, nonfounder_haplos) %>% 
    mutate(blah = stringdist(founder_haplos, nonfounder_haplos)) %>% 
    group_by(nonfounder_haplos) %>% 
    summarise(diff_to_nearest_founder = min(blah))
  
  res1 <- res %>% filter(diff_to_nearest_founder ==1)
  res2_5 <- res %>% filter(diff_to_nearest_founder >1 & diff_to_nearest_founder <6)
  res6_10 <- res %>% filter(diff_to_nearest_founder > 5 & diff_to_nearest_founder <11)
  res10plus <- res %>% filter(diff_to_nearest_founder > 9)
  
  out <- tibble(total_unique = length(unique(aa$haplo_sequence)),
                total_founder = length(founder_haplos),
                total_1diff = nrow(res1),
                total_2_5_diff = nrow(res2_5),
                total_6_10_diff = nrow(res6_10),
                total_10plus_diff = nrow(res10plus),
                n_snps = nchar(founder_haplos[1]))
}

###########
# GLOBALS #
###########

founder_info_file <- "~/projects/kakapo-genomics/data/cohort_by_gen.txt"

# for(chrom in c(1,2,4,5,6,7,8,9,10,11,12,14:26)){
# 
# vcf_file <- paste0("/media/drive_6tb/projects/kakapo-genomics/output/11_beagle_phasing/beagle_S", chrom, ".vcf")
# 
# blocksize = 10000
# 
# 
# ########
# # MAIN #
# ########
# 
# WARNING: SLOW
# 
# vcf <- read_delim(vcf_file, delim = "\t", skip = 40) %>% 
#   select(-c(`#CHROM`, FORMAT, REF, ALT, QUAL, FILTER, INFO)) %>% 
#   mutate(across(!contains(c("ID", "POS")), ~str_extract(., "^[:graph:]{3}"))) 
# 
founders <- read_delim(founder_info_file, delim = ",") %>%
  filter(gen == 1) %>%
  pull(ID)
# 
# 
# row_lower_vals <- seq(1, max(vcf$POS), by =blocksize)
# 
# test_set <- lapply(row_lower_vals, make_haplos) %>%
#   reduce(bind_rows)
# 
# # test_set <- mclapply(row_lower_vals, make_haplos, mc.preschedule = TRUE, mc.cores = 16) %>% 
# #   bind_rows() 
# 
# bb <- test_set %>% 
#   rownames_to_column(var = "block") %>% 
#   mutate(block = as.numeric(block)*blocksize) %>% 
#   mutate(above_1 =  total_2_5_diff + total_6_10_diff + total_10plus_diff) %>% 
#   pivot_longer(cols = -c(block, n_snps)) %>% 
#   filter(name %in% c("total_founder", "total_1diff", "above_1"))
# 
# write_delim(bb, paste0("~/Downloads/bb_S", chrom,".txt"), delim = "\t") 
# }


for(chrom in c(1,2,4,5,6,7,8,9,10,11,12,14:26)){

bb <- read_delim(paste0("~/Downloads/bb_S", chrom,".txt"), delim = "\t")

bb$name <- factor(bb$name, levels = c("above_1", "total_1diff", "total_founder"), labels = c("Multiple changes","Single change","Founder"))

main <- ggplot(bb, aes(x = block, y = value, fill = name))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#452736","#ca766aff","#dbcbbaff"))+
  labs(x = "Chromosome position",  y = "Number of unique haplotypes", fill = "Category", title = paste("Chromosome", chrom)) +
  theme_classic()+ 
  theme(plot.title = element_text(size = rel(1.6)),
        legend.position = "right",
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.6)),
        axis.title = element_text(size = rel(1.6)),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)))

sub <- ggplot(bb, aes(x = block, y = n_snps))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c("#452736"))+
  labs(x = "Chromosome position",  y = "SNPs in block") +
  theme_classic()+ 
  theme(legend.position = "right",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)))


main/sub +plot_layout(heights = c(7,1))

ggsave(paste0("window_haplos_chr", chrom, ".pdf"),
plot = last_plot(),
device = "pdf",
path = "~/Pictures",
scale = 1,
width = 400,
height = 200,
units = "mm")

}


overall_stats <- lapply(c(1,2,4:12,14:26), function(chrom){read_delim(paste0("~/Downloads/bb_S", chrom,".txt"), delim = "\t") %>% 
    mutate(chromosome = chrom)}) %>% 
  reduce(bind_rows) #%>% 
  # pivot_wider(names_from = name, values_from =value)


ggplot(overall_stats, aes(x = n_snps, y = value))+
  geom_point()+
  facet_wrap(~name)

os <- overall_stats %>% 
  group_by(block, chromosome) %>% 
  summarise(tot_uniq_haplos = sum(value), n_snps = n_snps)


ggplot(os, aes(x = n_snps, y = tot_uniq_haplos))+
  geom_point()


os2 <- overall_stats %>% 
  group_by(chromosome, name) %>% 
  mutate(n_snps = replace_na(n_snps, 0)) %>% 
  summarise(n_uniq_haplos = mean(value), n_snps = mean(n_snps))


os3 <- overall_stats %>% 
  mutate(n_snps = replace_na(n_snps, 0)) %>% 
  group_by(name) %>% 
  summarise(n_uniq_haplos = mean(value), n_snps = mean(n_snps))

overall_stats %>% 
  # filter(chromosome == 1) %>% 
  group_by(block, chromosome) %>% 
  summarise(uniq_haplos_per_block = sum(value)) %>% 
  ungroup() %>% 
  summarise(mean(uniq_haplos_per_block))


#### Unique haplos by window size (S1 only)

chrom <-  12


  vcf_file <- paste0("/media/drive_6tb/projects/kakapo-genomics/output/11_beagle_phasing/beagle_S", chrom, ".vcf")


  vcf <- read_delim(vcf_file, delim = "\t", skip = 40) %>%
    select(-c(`#CHROM`, FORMAT, REF, ALT, QUAL, FILTER, INFO)) %>%
    mutate(across(!contains(c("ID", "POS")), ~str_extract(., "^[:graph:]{3}")))

  founders <- read_delim(founder_info_file, delim = ",") %>%
    filter(gen == 1) %>%
    pull(ID)

get_n_uniq <- function(blocksize){

  row_lower_vals <- seq(1, max(vcf$POS), by =blocksize)

  test_set <- lapply(row_lower_vals, make_haplos) %>%
    reduce(bind_rows)

  # test_set <- mclapply(row_lower_vals, make_haplos, mc.preschedule = TRUE, mc.cores = 16) %>%
  #   bind_rows()

  bb <- test_set %>%
    rownames_to_column(var = "block") %>%
    mutate(block = as.numeric(block)*blocksize) %>%
    mutate(above_1 =  total_2_5_diff + total_6_10_diff + total_10plus_diff) %>%
    pivot_longer(cols = -c(block, n_snps)) %>%
    filter(name %in% c("total_founder", "total_1diff", "above_1")) %>% 
    group_by(block) %>% 
    summarise(uniq_haplos_per_block = sum(value)) %>% 
    ungroup() %>% 
    summarise(mean(uniq_haplos_per_block)) %>% 
    mutate(blocksize = blocksize)
  
  return(bb)
  
  
}


# WARNING: SLOW
library(stringdist)
for (the_length in c(200000, 300000)){
block_length_res <- lapply(the_length, get_n_uniq) %>% 
  reduce(bind_rows) %>% 
  write_delim(paste0("~/Downloads/S", chrom, "_", the_length, "_blocksize.txt"), delim = "\t")
}




S12_blocksize <- read_delim("~/Downloads/S12_blocksize.txt", delim = "\t") %>% 
  filter(blocksize != "blocksize")


ggplot(S12_blocksize, aes(x = as.numeric(blocksize), y = as.numeric(`mean(uniq_haplos_per_block)`))) +
  # geom_point() +
  geom_line(size = 1.2)+
  ylim(5,10) +
  labs(x = "Block size", y = "Mean number of unique haplotypes per block") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.2)))
  