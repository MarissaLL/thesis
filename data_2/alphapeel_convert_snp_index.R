library(tidyverse)


args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
ncycles <- args[2]
parent <- args[3]




convert_snp_index <- function(chr, ncycles, parent){
cpts <- read_delim(paste0("cpt_pos_S", chr, "_ncycles_", ncycles, "_", parent,".txt"), 
                   delim = "\t", 
                   col_names = FALSE) %>% 
  select(-X1) %>% 
  rename(snp_index = X2, 
         individual = X3, 
         focal_parent = X5,
         chromosome = X4)
vcf_pos <- read_delim(paste0("S", chr, "_pos.txt"), delim = "\t", col_names = FALSE) %>% 
  rownames_to_column()%>% 
  mutate(rowname = as.numeric(rowname)) %>% 
  rename(chromosome = X1,
         chr_position = X2)

cpts_vcf <- left_join(cpts, vcf_pos, by = c(snp_index = "rowname", "chromosome")) 
return(cpts_vcf)
}


aaa <- convert_snp_index(chr, ncycles, parent)

write_delim(aaa, paste0("recombination_S", chr, "_ncycles_", ncycles, "_", parent, ".txt"))

# all_mat <- lapply(chrs, convert_snp_index, ncycles = 20, parent = "mat") %>% 
#   reduce(bind_rows) 


# write_delim(all_mat, "/media/drive_6tb/projects/kakapo-genomics/output/21_alphapeel/recombination_mat_chr.txt", delim = '\t')





# 
# S8_vcf <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/21_alphapeel/S8_pos.txt", delim = "\t", col_names = FALSE) %>% 
#   rownames_to_column() %>% 
#   mutate(rowname = as.numeric(rowname))
# S9_vcf <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/21_alphapeel/S9_pos.txt", delim = "\t", col_names = FALSE) %>% 
#   rownames_to_column()%>% 
#   mutate(rowname = as.numeric(rowname))
# 
# S8_cpts <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/21_alphapeel/positions/S8_recombos_mat_pat.txt", delim = "\t", col_names = FALSE)
# S9_cpts <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/21_alphapeel/positions/S9_recombos_mat_pat.txt", delim = "\t", col_names = FALSE)
# 
# 
# S8_cpts_vcf <- left_join(S8_cpts, S8_vcf, by = c(X2 = "rowname", X4 = "X1")) %>% 
#   rename(n_cpts = X1, 
#          snp_index = X2, 
#          individual = X3, 
#          focal_parent = X5,
#          chromosome = X4, 
#          chr_position = X2.y) %>% 
#   select(n_cpts, snp_index, individual, focal_parent, chromosome, chr_position) %>% 
#   write_delim("/media/drive_6tb/projects/kakapo-genomics/output/21_alphapeel/positions/S8_recombos_chrpos.txt", delim = "\t")
# 
# 
# S9_cpts_vcf <- left_join(S9_cpts, S9_vcf, by = c(X2 = "rowname", X4 = "X1")) %>% 
#   rename(n_cpts = X1, 
#          snp_index = X2, 
#          individual = X3, 
#          focal_parent = X5,
#          chromosome = X4, 
#          chr_position = X2.y) %>% 
#   select(n_cpts, snp_index, individual, focal_parent, chromosome, chr_position) %>% 
#   write_delim("/media/drive_6tb/projects/kakapo-genomics/output/21_alphapeel/positions/S9_recombos_chrpos.txt", delim = "\t")
