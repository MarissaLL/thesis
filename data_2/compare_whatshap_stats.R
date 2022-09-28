library(tidyverse)

methods <- c("alphapeel", "beagle", "fimpute", "whatshap")

all <- lapply (methods, function(method){read_delim(paste0("/media/drive_6tb/projects/kakapo-genomics/output/29_whatshap_compare/whatshap_stats/", method, "_all_chr.tsv"), delim = "\t") %>% 
    mutate(method = method)}) %>% 
  reduce(bind_rows)


b <- all %>% 
  filter(chromosome == "ALL") %>% 
  mutate(phased/variants)
        