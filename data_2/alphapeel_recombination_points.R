#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(changepoint))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(tidyverse,
                                       warn.conflicts = FALSE))





# args <- commandArgs(trailingOnly = TRUE)
# chromosome <- args[1]
# ncycles <- args[2]
# inherited_chr <- args[3]


## Dev
chromosome <- "S1"
ncycles <- 20
inherited_chr <- "pat_only"
added_haplos <- read_delim(paste0("/media/drive_6tb/projects/kakapo-genomics/output/25_tidy_alphapeel/added_haplos/added_haplos_", chromosome, "_ncycles", ncycles, ".txt"), delim = "\t")
# added_haplos <- read_delim(paste0("/media/drive_6tb/projects/kakapo-genomics/output/21_alphapeel/long_but_broken/added_haplos_", chromosome, "_ncycles", ncycles, ".txt"), delim = "\t")
# recombo_birds <- "Adelaide"
# recombo_birds <- read_delim("/media/drive_6tb/projects/kakapo-genomics/output/25_tidy_alphapeel/phasing_filter/indivs_to_keep_S1_ncycles_20_pat.txt", delim = '\t') %>% 
#   pull()
# recombo_birds <- unique(added_haplos$chick)
recombo_birds <- c("Egilsay", "Konini", "Sinbad")

parent <- str_remove(string = inherited_chr, pattern = "_only")
message(paste("Processing chromosome", chromosome, "ncycles", ncycles, parent))


# added_haplos <- read_delim(paste0("added_haplos_", chromosome, "_ncycles", ncycles, ".txt"), delim = "\t")
# 
# # Birds that passed the phasing quality filter for this chromosome/ncycles/parent
# recombo_birds <- read_delim(paste0("indivs_to_keep_", chromosome, "_ncycles_", ncycles, "_", parent, ".txt"), delim = "\t") %>%
#   pull()


phased_only <- added_haplos %>% 
  filter(chick %in% recombo_birds == TRUE)

only_095 <- phased_only %>% 
  filter(sum_val >= 0.95)



################## Find changepoints and filter low quality ones #############

filter_cpts <- function(indiv, make_plot, return_cpt_pos, focal_parent){
  perc_opp_allowed <- 0.2  # This doesn't really do anything other than affect the display of the second panel of the plot (makes res easier to interpret).
  sw_size <-  50 # This is only for the display of the second panel in the plot
  
  wotp <- only_095 %>% 
    filter(chick == !!(indiv)) %>% 
    filter(focal_parent == !!(focal_parent)) %>% 
    mutate(inherited = case_when(inherited == "mat" ~ "maternal", # To look nicer on the plots
                                 inherited == "pat" ~ "paternal"))
  
  newest <- wotp %>%
    mutate(inherited2 = case_when(inherited == "maternal" ~ 1,
                                  inherited == "paternal" ~ 0)) %>% 
    arrange(pos) %>%
    mutate(lagged_mean = rollmean(x = inherited2, sw_size, align = "right", fill = NA)) %>% 
    mutate(good_rgn_lagged = case_when(lagged_mean <= perc_opp_allowed ~ "paternal",
                                       lagged_mean >= (1-perc_opp_allowed) ~ "maternal",
                                       TRUE ~ "uncertain")) %>% 
    mutate(leading_mean = rollmean(x = inherited2, sw_size, align = "left", fill = NA)) %>% 
    mutate(good_rgn_leading = case_when(leading_mean < perc_opp_allowed ~ "paternal",
                                        leading_mean > (1-perc_opp_allowed) ~ "maternal",
                                        TRUE ~ "uncertain")) %>% 
    filter(is.na(lagged_mean) == FALSE) %>%
    rownames_to_column(var = "name") %>% 
    mutate(name = as.numeric(name))
  
  
  
  # Find changepoints
  cpt_mean_res <- cpt.mean(data = newest$inherited2, method = "PELT")
  
  
  
  
  # Rearrange and filter the changepoint output
  # Convert changepoint locations to pos
  pts <- cpt_mean_res@cpts %>%
    head(-1) # last changepoint always appears to be on the last datapoint
  
  pts_pos <- newest %>%
    filter(newest$name %in% pts == TRUE) %>%
    select(pos) %>%
    pull()
  
  prev_means <- cpt_mean_res@param.est$mean %>% 
    head(n = -1)
  
  next_means <- cpt_mean_res@param.est$mean %>% 
    tail(n = -1)
  
  
  ###  Remove segments shorter than 5000 snps long, keep only the first changepoint and update the mean behind it to be the one of the next sufficiently long segment.
  pos_mean_cpts <- pts_pos %>% 
    enframe() %>% 
    mutate(prev_mean = prev_means,
           next_mean = next_means,
           prev_value = lag(value, default = 0),
           next_value = lead(value, default = max(newest$pos))) %>% 
    mutate(prev_seg_size = value - prev_value) %>% 
    mutate(next_seg_size = next_value - value) %>% 
    mutate(step = "unfiltered")
  
  
  filteredpos_cpts <- pos_mean_cpts %>% 
    filter(prev_seg_size > 2500) %>% 
    mutate(update_nxt = lead(prev_mean)) %>% 
    mutate(update_nxt = case_when(is.na(update_nxt)== TRUE ~ last(pos_mean_cpts$next_mean),
                                  TRUE ~ update_nxt)) %>% 
    mutate(step = "pos_filter") %>% 
    filter((value == max(pos_mean_cpts$value) & next_seg_size < 2500) == FALSE)
  
  
  filtered_cpts <- filteredpos_cpts %>% 
    filter(prev_mean > 0.75 & update_nxt < 0.25 |
             prev_mean < 0.25 & update_nxt > 0.75) %>% 
    select(value) %>% 
    pull()
  
  
  filtered5000 <- pos_mean_cpts$value[pos_mean_cpts$value %in% filteredpos_cpts$value == FALSE]
  filtered_mean <- pts_pos[pts_pos %in% c(filtered_cpts, filtered5000) == FALSE]

  
  cpt_current_res <- tibble(indiv = indiv,
                            changepoints = length(filtered_cpts),
                            approach = "filter",
                            chromosome = chromosome,
                            focal_parent = focal_parent)
  
  if(return_cpt_pos == TRUE){
    cpt_current_pos <- filtered_cpts %>% 
      enframe(value = "cpt_pos") %>% 
      mutate(indiv = indiv,
             chromosome = chromosome,
             focal_parent = focal_parent) %>% 
      write_delim(paste0("cpt_pos_", indiv, "_", chromosome, "_ncycles", ncycles, "_", parent, ".txt"), delim = "\t")
  }
  
  
  
  # cpt_res <- bind_rows(cpt_current_res, cpt_res)
  message(paste0(length(filtered_cpts), " changepoints detected in ", indiv, " for chromosome ", chromosome," ncycles ", ncycles, "."))
  
  if(make_plot == TRUE){
    
    raw_origin_colours <- c(maternal = "#452736", paternal = "#dbcbbaff")
    p1 <- ggplot(wotp, aes(x = pos, y = sum_val, colour = inherited))+
      labs(x = "SNP number", y = "Summed probability", colour = "SNP origin")+
      geom_point() +
      theme_classic() +
      theme(legend.position = "right",
            axis.title.y = element_text(size = rel(1.4)),
            axis.title.x = element_blank(),
            panel.background = element_rect(fill = NA, 
                                            colour = "black", 
                                            size = 0.7 ),
            axis.text.y = element_text(size = rel(1.2)),
            # axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+ 
      scale_colour_manual(values = raw_origin_colours)
    
    
    haplo_colours <- c(maternal = "#452736", paternal = "#dbcbbaff", uncertain ="#ca766aff")
    p2 <- ggplot(newest, aes(x = pos, y = sum_val, colour = good_rgn_lagged))+
      geom_point()+
      labs(x = "SNP number", y = "Summed probability", colour = "SNP origin (window)")+
      theme_classic() +
      theme(legend.position = "right",
            axis.title.y = element_text(size = rel(1.4)),
            axis.title.x = element_blank(),
            panel.background = element_rect(fill = NA, 
                                            colour = "black", 
                                            size = 0.7 ),
            axis.text.y = element_text(size = rel(1.2)),
            # axis.text.x = element_blank(),
            axis.ticks.x = element_blank())+
      scale_colour_manual(values = haplo_colours)
    
    p7 <- ggplot(newest, aes(x = pos, y = lagged_mean))+
      labs(x = "SNP number", y = "Window mean") +
      geom_point() +
      ylim(0,1) +
      theme_classic() +
      theme(axis.title = element_text(size = rel(1.4)),
            panel.background = element_rect(fill = NA, 
                                            colour = "black", 
                                            size = 0.7 ),
            axis.text = element_text(size = rel(1.2))) +
      geom_vline(xintercept = filtered_cpts, colour = "red", size = 0.8)
    
    p1 / p2 / p7
    
    ggsave(paste0("~/Pictures/sw_cpts_pat2/", indiv, "_sw", sw_size, "_cpts_pat_filtered_", chromosome, "_ncycles", ncycles, ".eps"),
           width = 12, 
           height = 7,
           units = "in")
  }
    return(cpt_current_res)
}

#### Don't stop if an indiv doesn't have enough SNPs to detect changepoints
filter_cpts_ignore_unphased <- function(indiv, make_plot, return_cpt_pos, focal_parent) {
  tryCatch({result <- filter_cpts(indiv, make_plot, return_cpt_pos, focal_parent);}, 
           error = function(e) {result <<- NULL
                                message(paste0("Error in ", indiv, ", not enough phased SNPs to detect recombination for chromosome ", chromosome, " ncycles ", ncycles))});
  return(result)
  }


bbb <- lapply(recombo_birds, filter_cpts_ignore_unphased, make_plot = TRUE, return_cpt_pos = TRUE, focal_parent = inherited_chr ) %>% 
  reduce(bind_rows)


# write_delim(bbb, paste0("changepoints_", chromosome, "_ncycles", ncycles, "_", parent, ".txt"), delim = "\t")
