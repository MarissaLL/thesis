library(tidyverse)
library(patchwork)

founders <- c("Alice", "Arab", "Barnard", "Basil", "Bella", "Ben", "Bill", 
              "Blades", "Bonus", "Boss", "Cyndy", "Felix", "Flossie", "Fuchsia", 
              "Gumboots", "Gunner", "Heather", "Jean", "Jimmy", "Joe", 
              "Lee", "Lionel", "Lisa", "Luke", "Maggie", "Margaret-Maree", 
              "Merv", "Nog", "Nora", "Ox", "Piripi", "Rangi", "Richard_Henry", 
              "Ruth", "Sandra", "Sass", "Smoko", "Solstice", "Stumpy", "Sue", 
              "Suzanne", "Waynebo", "Whiskas")

haplos <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/31_haploblocker/beagle_coverage_90/coverage_90_haploblocks.vcf", delim = "\t", skip = 31) %>% 
  mutate(blockid = paste0("S", `#CHROM`, "_", ID)) %>% 
  select(blockid, all_of(founders)) %>% 
  mutate(across(!blockid, ~ case_when(. == "0/0" ~ "0",
                                                TRUE ~ "1")))





test_all <- c()
res <- tibble(name = "test", num = 0)

for(name in founders){
  test <- haplos %>% 
    filter((blockid %in% test_all) == FALSE) %>% 
    filter(!!sym(name) == 1) %>% 
    pull(blockid)
  test_all <- c(test_all, test)
  res <- res %>% 
    add_row(name = name, num = length(test_all))
}


res <- res %>% 
  filter(name != "test")

ggplot(res, aes(x = name, y = num, group = "a"))+
  geom_line()+
  labs(x = "Founder",
       y = "Number of haplotype blocks in population") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



occurrance <- haplos %>% 
  pivot_longer(cols = -blockid) %>% 
  group_by(blockid) %>% 
  summarise(occurrance = sum(as.numeric(value))) 


a <– ggplot(occurrance, aes(x = occurrance, fill = "a"))+
  geom_histogram(binwidth = 1) +
  labs(x = "Number of founders haplotype block occurs in",
       y = "Number of haplotype blocks") +
  scale_fill_viridis_d()+
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))

unique_blocknames <- occurrance %>% 
  filter(occurrance == 1) %>% 
  pull(blockid)

unique_blocks <- haplos %>% 
  pivot_longer(cols = -blockid) %>% 
  filter(blockid %in% unique_blocknames & value == "1")

n_unique_blocks <- unique_blocks %>% 
  group_by(name) %>% 
  summarise(n = n())

all_founders_n_unique <- enframe(founders,name = NULL, value = "name") %>% 
  full_join(n_unique_blocks) %>% 
  mutate(n = case_when(is.na(n) ~ 0,
                       TRUE ~ as.numeric(n)))


ggplot(all_founders_n_unique, aes(x = name, y = n, fill = "a"))+
  geom_histogram(stat = "identity")+
  labs(x = "Founder",
       y = "Number of unique haplotype blocks") +
  scale_fill_viridis_d()+
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Split by unique Rakiura vs RH
unique_by_location <- haplos %>% 
  pivot_longer(cols = -blockid) %>% 
  mutate(loc = case_when(name == "Richard_Henry" ~ "Fiordland",
                         name != "Richard_Henry" ~ "Rakiura")) %>% 
  group_by(blockid, loc) %>% 
  summarise(occurrances = sum(as.numeric(value))) %>% 
  pivot_wider(names_from = loc, values_from = occurrances) %>% 
  mutate(category = case_when(Fiordland == 0 & Rakiura > 0 ~ "Rakiura only",
                              Fiordland > 0 & Rakiura > 0 ~ "Both",
                              Fiordland > 0 & Rakiura == 0 ~ "Fiordland only",
                              Fiordland == 0 & Rakiura == 0 ~ "Non-founder",
                              TRUE ~ "error case")) 

summary_unique_by_location <-  unique_by_location %>% 
  group_by(category) %>% 
  summarise(n = n())

ggplot(subset(summary_unique_by_location, category != "Non-founder"), aes(x = category, y = n, fill = "a"))+
  geom_histogram(stat = "identity")+
  labs(x = "Founder population present in",
       y = "Number of haplotype blocks") +
  scale_fill_viridis_d()+
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))

