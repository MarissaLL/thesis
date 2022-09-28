library(sequoia)
library(tidyverse)


#########
# INPUT #
#########

files <- commandArgs(trailingOnly = TRUE)

ds_geno_file <- files[1]
birds_file <- files[2]

##Dev
# ds_geno_file <- "/media/drive_6tb/projects/kakapo-genomics/output/37_sequoia_pt2/sample_400_deepvariant.r1.filtered.raw"
# birds_file <- "/media/drive_6tb/projects/kakapo-genomics/output/37_sequoia_pt2/birds.txt"

########
# MAIN #
########

converted <- GenoConvert(ds_geno_file)
birds_data <- read.table(birds_file,  sep = " ",  skip = 1)
results <- sequoia(GenoM = converted, LifeHistData = birds_data, Tfilter = -500)

outname <- str_extract(ds_geno_file, "sample_[:digit:]+_rep[:digit:]+_[:alpha:]+_?[:alpha:]+")
write_delim(results$Pedigree, paste0("sequoia_", outname, ".txt"), delim = "\t")
 
 