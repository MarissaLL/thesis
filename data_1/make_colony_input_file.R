library(tidyverse)


#############
# FUNCTIONS #
#############

get_colony_formatted_gts <- function(name, dataset){
  gts <- dataset %>% 
    filter(IID == name) %>% 
    select(-IID) %>% 
    pivot_longer(cols = everything(), names_to = "SNP", values_to = "gt") %>% 
    mutate(new_gt = case_when(gt == 0 ~ "3 3",
                              gt == 1 ~ "3 4",
                              gt == 2 ~ "4 4",
                              is.na(gt) == TRUE ~ "0 0")) %>% 
    pull(new_gt)
  gts_again <- paste(gts, collapse = " ")
  return(paste(name,gts_again, sep = " "))
}

###########
# GLOBALS #
###########
cmdargs <- commandArgs(trailingOnly = TRUE)


project_name <- cmdargs[1]
n_snp <- cmdargs[2]
rep <- cmdargs[3]
dataset <- cmdargs[4]
# project_name <- "kakapo_colony_pedigree_test"
# n_snp <- 200
output_name <- paste0("colony_out_sample", n_snp, "_rep",rep)
# gts_filename <- paste0("/media/drive_6tb/projects/kakapo-genomics/output/37_sequoia_pt2/sample_1000_rep2_gwas_set.raw")
gts_filename <- paste0("sample_", n_snp,"_rep", rep, "_", dataset, ".raw")

########
# MAIN #
########

gts_file <- read_delim(gts_filename, delim = " ", trim_ws = TRUE) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE, -last_col())

marker_ids <- names(gts_file)[names(gts_file) %in% "IID" == FALSE]%>% 
  str_remove('022045') %>% 
  str_remove('ctg1_') %>% 
  str_remove('_[:alpha:]_[:alpha:]_[:alpha:]$')

birds <- read_delim("bird_info_rships_nov20.csv", delim = ",")

# offspring_names <- c("Hoki", "Te_Kingi")
# mother_candidate_names <- c("Zephyr", "Bella", "Alice", "Jean")
# father_candidate_names <- c("Arab", "Felix", "Sass", "Waynebo")

offspring_names <- birds %>% filter(is.na(father) == FALSE) %>% pull(ID)
mother_candidate_names <- birds %>%  filter(sex == "F") %>%  pull(ID)
father_candidate_names <- birds %>%  filter(sex == "M") %>%  pull(ID)

message("checkpoint 1")

offspring_id_gt <- lapply(offspring_names, get_colony_formatted_gts, dataset = gts_file) %>% 
  unlist()
mother_id_gt <- lapply(mother_candidate_names, get_colony_formatted_gts, dataset = gts_file) %>% 
  unlist()
father_id_gt <- lapply(father_candidate_names, get_colony_formatted_gts, dataset = gts_file) %>% 
  unlist()

message("checkpoint 2")

colony_file <- paste(paste0(project_name," !Project name"),
      paste0(output_name, " !Output name"),
      paste0(length(offspring_names)," ! Number of offspring in the OFS sample (Integer)"),
      paste0(length(marker_ids), " ! Number of loci (Integer)"),
      "276453 ! Seed for random number generator",
      "0 ! Update allele frequency????????????????????????????????????????",
      "2 ! Dioecious",
      "1 ! Inbreeding",
      "0 ! Diploid",
      "0 0 ! Polygamy",
      "0 ! No clone inference",
      "0 ! Don't scale full sibship",
      "0 ! No sibship prior",
      "0 ! Unknown population allele frequency",
      "10 ! Number of runs",
      "1 ! Short run",
      "0 ! 0/1=Monitor method by Iterate#/Time in second",
      "10000 ! Monitor interval in Iterate# / in seconds",
      "0 ! 0/1=No/Yes for run with Windows GUI",
      "1 ! 0/1/2=PairLikelihood score/Fulllikelihood/FPLS",
      "2 ! 0/1/2/3=Low/Medium/High/Very high precision with Fulllikelihoo",
      paste0(paste(marker_ids, collapse = " "), " !Marker Ids"),
      "0@ !Marker types, 0/1=Codominant/Dominant. @ indicates all are the same",
      "0.0000@ !Allelic dropout rate at each locus",
      "0.0001@    !Other typing error rate at each locus",
      paste0(paste(offspring_id_gt, collapse = "\n")," !Offspring ID and genotypes at locus"),
      "0.9  0.9     !probabilities that the father and mother of an offspring are included in candidates",
      paste0(length(father_candidate_names), " ", length(mother_candidate_names), " !Numbers of candidate males and females"),
      paste0(paste(father_id_gt, collapse = "\n"), " !Candidate male ID and genotypes at locus"),
      paste0(paste(mother_id_gt, collapse = "\n"), " !Candidate female ID and genotypes at locus"),
      "0 0 !Number of offspring with known paternity, exclusion threshold",
      "0 0 !Number of offspring with known maternity, exclusion threshold",
      "0 !Number of known paternal sibship",
      "0 !Number of known maternal sibship",
      "0 !Number of offspring with known excluded paternity",
      "0 !Offspring ID, number of excluded males, the IDs of excluded males",
      "0 !Numberof offspring with known excluded maternity",
      "0 !Offspring ID, number of excluded females, the IDs of excluded females",
      "0 !Number of offspring with known excluded paternal sibships",
      "0 !Number of offspring with known excluded maternal sibships",
  sep= "\n")


write_lines(colony_file, paste0("colony_sample_", n_snp,"_rep",rep, ".dat"))
