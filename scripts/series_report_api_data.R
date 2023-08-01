## ---------------------------
## Script name: series_report_api_data.R
##
## Purpose of script: Pull report CSVs directly from REDCap for filtering and subsetting
##
## Author: Nancy Scott
##
## Date Created: 2022-12-08
##
## Email: scot0854@umn.edu
## ---------------------------
## Notes:
## ---------------------------
## load packages
library(tidyverse)
library(lubridate)
library(magrittr)
library(writexl)

# redcap report IDs
avail_seq_data <- '52998'
samples <- '41047'
mic_results <- '41628'
chef_data <- '52263'
spot_plates <- '56394'
genes <- '56393'

# function to import report from redcap
import_report <- function(report_number) {
  token <- ""
  url <- "https://redcap.ahc.umn.edu/api/"
  formData <- list("token"=token,
                   content='report',
                   format='csv',
                   report_id=report_number,
                   csvDelimiter='',
                   rawOrLabel='label',
                   rawOrLabelHeaders='raw',
                   exportCheckboxLabel='true',
                   returnFormat='csv'
  )
  response <- httr::POST(url, body = formData, encode = "form")
  result <- httr::content(response, show_col_types = FALSE)
}

# Sample ID, species, series and cluster IDs
sample_info <- import_report(samples) %>%
    select(-c(starts_with('redcap_repeat')))

# snp data
gene_vars <- import_report(genes) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gene, protein_change, alt_freq)  %>%
    pivot_wider(names_from = "redcap_repeat_instance",
                values_from = c("gene", "protein_change", "alt_freq"))

# CHEF gel results
chef_done <- import_report(chef_data) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gel_date, blot_prepared) %>%
    pivot_wider(names_from = "redcap_repeat_instance",
                values_from = c("gel_date", "blot_prepared"))

# Missing CHEF results
todo <- sample_info %>%
    filter(isolate_type == "clinical") %>%
    left_join(chef_done, by="primary_id") %>%
    filter(is.na(redcap_repeat_instrument))

# export chef-todo
write_xlsx(todo,paste0(Sys.Date(),"_CHEF_todo.xlsx"))

# MSI location if sequencing data exists
seq_info <- import_report(avail_seq_data) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, msi_path_r1, msi_path_r2, msi_nanopore_path) %>%
    pivot_wider(names_from = "redcap_repeat_instance",
                values_from = c("msi_path_r1", "msi_path_r2", "msi_nanopore_path"))

# MIC and SMG results
mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, drug, start_date, mic50, smg) %>%
    pivot_wider(names_from = "redcap_repeat_instance", 
                values_from = c("drug","start_date","mic50","smg"))

# all_data_merged
all_samples <- left_join(sample_info, gene_vars) %>%
    left_join(mic_info) %>%
    left_join(chef_done) %>%
    filter(primary_id != "MEC103") 

write_xlsx(all_samples,paste0(Sys.Date(),"merged_Candida_data.xlsx"))

all_samples$cluster_id <- as.factor(all_samples$cluster_id)
all_samples$series_id <- as.factor(all_samples$series_id)

# calculate timeframes and isolate numbers for each series or cluster
# and append to all_samples tibble
# turn this into a function with series or cluster passed as parameter
serial_timespans <- all_samples %>%
    filter(series_id !="NA") %>%
    group_by(series_id) %>%
    mutate(series_span = max(collection_date) - min(collection_date) + 1)

count_serial_samples <- serial_timespans %>%
    group_by(series_id) %>%
    filter(series_span >0) %>%
    summarize(series_count = n())

serial_timespans <- serial_timespans %>%
    inner_join(count_serial_samples, by = 'series_id')
    
serial_calculations <- serial_timespans %>%
    select(primary_id, series_id, series_span, series_count) 

cluster_timespans <- all_samples %>%
    filter(cluster_id != "NA") %>%
    group_by(cluster_id) %>%
    mutate(cluster_span = max(collection_date) - min(collection_date) + 1)

cluster_count <- cluster_timespans %>%
    group_by(cluster_id) %>%
    summarize(cluster_count=n())

cluster_timespans <- cluster_timespans %>%
    inner_join(cluster_count, by="cluster_id")

cluster_calcs <- cluster_timespans %>%
    select(primary_id, cluster_id, cluster_span, cluster_count)

all_samples <- all_samples %>%
    full_join(serial_calculations, by = c("primary_id", "series_id")) %>%
    full_join(cluster_calcs, by = c("primary_id", "cluster_id")) 

# priorities: rare isolates, longer than one day, recurrent or clusters
priority_series <- all_samples %>%
    filter(series_id !="NA") %>%
    #filter(genus_species %in% c("Candida krusei (Pichia kudriavzevii)",
     #                   "Clavispora (Candida) lusitaniae",
      #                  "Candida utilis (Cyberlindnera jadinii)",
       #                 "Candida kefyr (Kluyveromyces marxianus)") |
    filter(series_span > 1, series_span < 30, series_count >2) %>%
    group_by(series_id) %>%
    slice(2L) %>%
    ungroup() %>%
    select(genus_species, series_id, series_span, series_count) %>%
    arrange(genus_species, desc(series_span))

recurrent <- all_samples %>%
    filter(series_id != "NA" & series_span >30) %>%
    group_by(series_id) %>%
    slice(2L) %>%
    ungroup() %>%
    select(genus_species, series_id, series_span, series_count) %>%
    arrange(genus_species, desc(series_span))

# cluster info
clusters <- all_samples %>%
    filter(cluster_id != "NA") %>%
    group_by(cluster_id) %>%
    slice(1L) %>%
    ungroup %>%
    select(genus_species, cluster_id, cluster_span, cluster_count) %>%
    arrange(genus_species, desc(cluster_span))


## little plot of samples by species
species_count <- all_samples %>%
    group_by(genus_species) %>%
    summarize(species_count=n()) %>%
    arrange(desc(species_count))

species_colors <- c("#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
                     "#882255","#BBBBBB", "#AA4499", "#DDCC77", "black")

sp_plot <- ggplot(species_count, aes(x=genus_species, y=species_count)) +
    geom_col(fill=species_colors) +
    scale_x_discrete(limits=species_count$genus_species,
                     labels = c("C. albicans", 
                                "C. glabrata",
                                "C. parapsilosis",
                                "C. lusitaniae",
                                "C. krusei",
                                "C. kefyr",
                                "C. orthopsilosis",
                                "C. tropicalis",
                                "C. utilis",
                                "C. dubliniensis",
                                "C. nivariensis")) +
    ylab("Number of Isolates") +
    xlab("Species") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 30, 
                                     hjust = 1, 
                                     vjust = 1,
                                     color = "black",
                                     size = 14,
                                     face = "italic")) +
    theme(axis.title = element_text(color = "black", size = 16))

ggsave("images/2022_Candida_MEC_spp_distribution.tiff", sp_plot, width = 7, height = 4, units = "in")

# export all data to csv
all_samples %>% write.csv("redcap_all_samples.csv")

# count number of patients per spp.

inds <- all_samples %>%
    group_by(genus_species) %>%
    count(series_id) %>%
    mutate(total=sum(n))

# create vectors for venn diagram of cluster-series-na pts.

neither <- all_samples %>%
    filter(is.na(series_id) & is.na(cluster_id) & genus_species == "Candida albicans") %>%
    pull(primary_id)

series <- all_samples %>%
    filter(series_id != "NA" & is.na(cluster_id) & genus_species == "Candida albicans") %>%
    pull(primary_id)

cluster <- all_samples %>%
    filter(is.na(series_id) & cluster_id != "NA" & genus_species == "Candida albicans") %>%
    pull(primary_id)

both <- all_samples %>%
    filter(series_id != "NA" & cluster_id != "NA" & genus_species == "Candida albicans") %>%
    pull(primary_id)

overlaps <- all_samples %>%
    mutate(series_sample = (series_id != "NA" & is.na(cluster_id))) %>%
    mutate(cluster_sample = (cluster_id != "NA" & is.na(series_id))) %>%
    mutate(both_sample = (cluster_id != "NA" & series_id != "NA")) %>%
    mutate(neither_sample = is.na(cluster_id) & is.na(series_id)) %>%
    select(primary_id, genus_species, series_id, cluster_id, series_sample, cluster_sample, both_sample, neither_sample) %>%
    replace_na(list(series_sample=FALSE, cluster_sample=FALSE, both_sample=FALSE, neither_sample=FALSE))

# generate stacked bar plot per spp. of cluster and series sets
summary_overlaps <- overlaps %>%
    group_by(genus_species) %>%
    summarize(neither_sum=sum(neither_sample==TRUE), both_sum=sum(both_sample==TRUE),
              cluster_sum=sum(cluster_sample==TRUE),series_sum=sum(series_sample==TRUE))
              #total=sum(series_sum, cluster_sum, both_sum, neither_sum))

library(reshape2)

reshape_overlaps <- melt(summary_overlaps, id.var="genus_species")
viz_overlaps <- ggplot(reshape_overlaps, aes(x=genus_species, y=value, fill=variable)) + 
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("#EECC66", "#EE99AA", "#6699CC", "#004488"),
                      labels = c("One-off isolates", "Cluster and series", "Cluster isolate", "Series isolate")) +
    xlab("Species") +
    ylab("Count") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, 
                                     hjust = 1, 
                                     vjust = 1,
                                     color = "black",
                                     size = 10,
                                     face = "italic")) +
    theme(axis.title = element_text(color = "black", size = 16))
    
ggsave("images/2022_Candida_patient_categories.tiff", viz_overlaps, width = 11.5, height = 8, units = "in")
