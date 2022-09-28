library(tidyverse)
id_name_mapping <- read_delim("/media/drive_6tb/projects/novaseq_kakapo/names_ids_mapping.csv", delim = ",", col_names = FALSE)
adapter_metrics <- read_delim("/media/drive_6tb/projects/novaseq_kakapo/Reports/Adapter_Metrics.csv", delim = ",")
demux_stats <- read_delim("/media/drive_6tb/projects/novaseq_kakapo/Reports/Demultiplex_Stats.csv", delim = ",")
fastq_list <- read_delim("/media/drive_6tb/projects/novaseq_kakapo/Reports/fastq_list.csv", delim = ",")
index_hopping_counts <- read_delim("/media/drive_6tb/projects/novaseq_kakapo/Reports/Index_Hopping_Counts.csv", delim = ",")
unknown_barcodes <- read_delim("/media/drive_6tb/projects/novaseq_kakapo/Reports/Top_Unknown_Barcodes.csv", delim = ",")

# Perc of bases that were assignable to an individual (and not adapter)
a <- adapter_metrics %>% 
  mutate(tot_seq = R1_SampleBases + R2_SampleBases, 
         tot_adapt =  R1_AdapterBases + R2_AdapterBases)

b <- sum(a$tot_seq)+ sum(a$tot_adapt)

c <- a %>% mutate(perc = (tot_seq/b)*100)
sum(c$perc)- sum(c$perc[c$Sample_ID == "Undetermined"])

# Perc of reads that were assignable
adapter_metrics %>% mutate(category = case_when(Sample_ID == "Undetermined" ~ "Undetermined",
                           TRUE ~ "Assigned")) %>% 
  group_by(category) %>% 
  summarise(reads = sum(`# Reads`)) %>% 
  mutate(perc = reads/sum(reads))


# Yield per lane
d <- adapter_metrics %>% mutate(category = case_when(Sample_ID == "Undetermined" ~ "Undetermined",
                                                TRUE ~ "Assigned")) %>% 
  mutate(tot_seq = R1_SampleBases + R2_SampleBases) %>% 
  group_by(Lane,category) %>% 
  summarise(aa = sum(tot_seq)) 

### Demux stats
dmx <- demux_stats %>% 
  mutate(category = case_when(SampleID == "Undetermined" ~ "Undetermined",
                              TRUE ~ "Assigned"))

dmx %>% 
  group_by(category) %>% 
  summarise(mean(`Mean Quality Score (PF)`))


aaaa <- dmx %>% 
  filter(category == "Assigned") %>%
  group_by(SampleID) %>% 
  summarise(perf_idx = sum(`# Perfect Index Reads`), 
            tot_reads = sum(`# Reads`),
            one_mm_idx = sum(`# One Mismatch Index Reads`)) %>% 
  mutate(perc_perf =perf_idx/tot_reads) %>% 
  mutate(perc_oneoff = one_mm_idx/tot_reads) %>% 
  pivot_longer(cols = c(perc_oneoff, perc_perf)) %>% 
  left_join(id_name_mapping, by = c("SampleID" = "X1"))
  # summarise(mean(perc_perf), min(perc_perf), max(perc_perf))

  ggplot(aaaa, aes(x = X2, y = value, fill = name))+
    geom_histogram(stat = "identity", position = "fill")+
    scale_fill_manual(values =  c("#452736","#b69a96ff")) +
    labs(x = "Kakapo Name",  y = "Proportion of assigned reads", fill = "") +
    theme_classic()+ 
    theme(legend.position = "top",
          axis.title = element_text(size = rel(1.6)),
          panel.background = element_rect(fill = NA,
                                          colour = "black",
                                          size = 0.7 ),
          axis.text = element_text(size = rel(1.0)),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.ticks.x = element_blank())
bbbb <- aaaa %>% 
  filter(name == "perc_oneoff") %>%
  select(-name, -value) %>% 
  pivot_longer(cols = c(perf_idx, one_mm_idx)) %>% 
  arrange(tot_reads) %>% 
  mutate(value = value/1000000) %>% 
  mutate(name = case_when(name == "perf_idx" ~ "Exact index match",
                          name == "one_mm_idx" ~ "Single index mismatch"))

ggplot(bbbb, aes(x = reorder(X2, tot_reads), y = value, fill = name ))+
  coord_flip()+
  geom_histogram(stat = "identity", position = "stack")+
  theme_classic()+ 
  scale_fill_manual(values =  c("#b69a96ff", "#452736")) +
  labs(x = "Kakapo Name",  y = "Assigned reads (millions)", fill = "") +
  theme(legend.position = "top",
        axis.title = element_text(size = rel(1.6)),
        panel.background = element_rect(fill = NA,
                                        colour = "black",
                                        size = 0.7 ),
        axis.text = element_text(size = rel(1.0)),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x = element_blank())
  
# What are index hopping counts??
  
  
# Unknown barcodes
  sum(unknown_barcodes$`# Reads`)
  # There are far more undetermined reads in the demux stats than there are in the unknown barcodes file. Presumably these are variants that only occur a small number of times.
  demux_stats %>% filter(SampleID == "Undetermined") %>% 
    summarise(sum(`# Reads`))

    