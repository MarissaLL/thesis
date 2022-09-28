library(HaploBlocker)
library(tidyverse)

###########
# GLOBALS #
###########
target_coverage <- 95
# phased_file <- "/media/drive_6tb/projects/kakapo-genomics/output/25_tidy_alphapeel/haplotypes/alphapeel_phased_S"
phased_file <- "/media/drive_6tb/projects/kakapo-genomics/output/11_beagle_phasing/beagle_S"
outdir <- paste0("/media/drive_6tb/projects/kakapo-genomics/output/31_haploblocker/beagle_coverage_", target_coverage, "/")

################################
# Determine coverage vs nblock #
################################
# WARNING THIS IS VERY SLOW TO RUN

get_target_cov_info <- function(coverage){
  a <- block_calculation(dhm = "/media/drive_6tb/projects/kakapo-genomics/output/11_beagle_phasing/beagle2.vcf.gz",
                         target_coverage = coverage)
  tibble(coverage = coverage,
         n_blocks = length(a))
}

test_cov2 <- lapply(seq(from = 0.05, to = 1, by = 0.05), get_target_cov_info)

test_cov_for_plot <- test_cov2 %>% 
  reduce(bind_rows) %>% 
  add_row(coverage = 0, n_blocks = 0) %>% 
  mutate(perc_phased = coverage*100)

ggplot(test_cov_for_plot, aes(x = perc_phased, y = n_blocks, colour = "a"))+ 
  geom_line(size = 2)+
  scale_colour_viridis_d()+
  geom_vline(xintercept = 90)+
  xlim(0,100) +
  labs(x = "% of genome phased",
       y = "N haplotype blocks") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        panel.background = element_rect(fill = NA, 
                                        colour = "black", 
                                        size = 0.7 ),
        axis.text = element_text(size = 14))

# write_delim(test_cov_for_plot, "/media/drive_6tb/projects/kakapo-genomics/output/31_haploblocker/coverage_x_nblock.txt", delim = "\t")


#############################
# Generate haplotype blocks #
#############################


for (chrom in c(1,2,4:12, 14:26)) {

vcf_file <- paste0(phased_file, chrom, ".vcf")


blocklist <-  block_calculation(dhm = vcf_file, target_coverage = (target_coverage/100))

preprocessing_grep_command <- sprintf("grep %s %s", 
                                      paste0("POS"),
                                      paste(vcf_file, collapse = " "))


vcf_names <- data.table::fread(cmd = preprocessing_grep_command) %>% 
  select(-`#CHROM`, -POS, -REF, -ALT, -ID, -INFO, -QUAL, -FILTER, -FORMAT) %>% 
  names() 

haplo_names <- paste0(rep(vcf_names, each = 2), "_", rep(c(1,2), each = 1))

bb <-  block_matrix_construction(blocklist = blocklist) 

colnames(bb) <- haplo_names
write.table(bb, file=paste0(outdir, "haplotype_blocks_S", chrom, ".txt"), row.names=TRUE, col.names=TRUE, quote = FALSE)

pdf(file = paste0("~/Pictures/haplos_S", chrom, "_cov", target_coverage, ".pdf"), width = 8, height = 3.5)
plot_block(blocklist, type = "snp", orientation = "mid")
dev.off()

startend <- blocklist_startend(blocklist, type = "bp")

rownames(startend) <- str_replace_all(rownames(startend), " ", "")
write.table(startend, file = paste0(outdir, "haplotype_startend_S", chrom, ".txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)


}


####################
# Old notes/trials #
####################

# blocklist_plot_xsixe()
# testcoverage <- coverage_test(blocklist_alphapeel)


# Can't run for haplotypes from only one individual.
# block_calculation(dhm = "/media/drive_6tb/projects/kakapo-genomics/output/13_whatshap_phasing/Adelaide.phased.vcf")

# blocklist_whatshap <- block_calculation(dhm = "/media/drive_6tb/projects/kakapo-genomics/whatshap_nesi/S8.phased.vcf")
# 
# plot_block(blocklist_whatshap, type = "snp", orientation = "snp")
# startend <- blocklist_startend(blocklist_whatshap)
# coverage <- coverage_test(blocklist_whatshap)



#### NOTES
# Also to try extended block identification as mentioned in the paper. Requires t and also adaptive mode
# File S4 includes the R-code used to generate an exemplary MAGIC population for the section on recovering founder haplotypes. - Look at this too
# incorporate basepairs too

# Try block_windowdataset if I can assign haplotypes back to individuals because then I could make my plots of haplos in founders vs not on the window basis

# 
# raw <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/21_alphapeel/test_line.vcf", delim = '\t') 
# phased <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/25_tidy_alphapeel/haplotypes/test_line.vcf", delim = '\t') 
# 
# 
# blah <- full_join(raw, phased)



###### Adjusting target coverage
# blocklist_alphapeel_tc <-  block_calculation(dhm = "/media/drive_6tb/projects/kakapo-genomics/output/25_tidy_alphapeel/haplotypes/alphapeel_phased_S8.vcf",
#                                           target_coverage = 0.1)


# 
# aa <- test_cov2 %>% unlist() %>% 
#   enframe() %>%
#   mutate(iteration = rep(1:10,each = 2)) %>% 
#   pivot_wider(id_cols = iteration, names_from = name, values_from = value) %>% 
#   mutate(coverage = seq(from = 0.1, to = 1, by = 0.1))
# 
# ggplot(aa, aes(x = coverage, y = n_blocks))+
#   geom_point()+
#   geom_line()+
#   theme_classic()+
#   labs(x = "Target coverage", y = "Number of haplotype blocks")+
#   theme(legend.position = "none",
#       axis.title = element_text(size = 18),
#       panel.background = element_rect(fill = NA, 
#                                       colour = "black", 
#                                       size = 0.7 ),
#       axis.text = element_text(size = 14))
# 
# 


### Try big_output

# blocklist_alphapeel_big <-  block_calculation(dhm = "/media/drive_6tb/projects/kakapo-genomics/output/25_tidy_alphapeel/haplotypes/S8/alphapeel_phased_S8.vcf", big_output = TRUE)
# 
# 
# write_json(blocklist_alphapeel_big, "haploblocker_S8_alphapeel.json")


#### With bp_map

## 
# write_vcf(blocklist_alphapeel, path = "haploblocker_S8_alphapeel", hetero = TRUE)
# 
# vcf_res <- read_delim("haploblocker_S8_alphapeel.vcf", delim = "\t", skip = 4)
