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
avail_seq_data='40639'
samples='41047'
mic_results='41628'
chef_data <- '52263'

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

# CHEF gel results
chef_done <- import_report(chef_data) %>%
    filter(redcap_repeat_instrument != "NA")

# Missing CHEF results
todo <- sample_info %>%
    left_join(chef_done, by="mec_id") %>%
    filter(is.na(redcap_repeat_instrument))

# export chef-todo
write_xlsx(todo,paste0(Sys.Date(),"_CHEF_todo.xlsx"))

# MSI location if sequencing data exists
seq_info <- import_report(avail_seq_data) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(mec_id, msi_path_r1, msi_path_r2)

# MIC and SMG results
mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(mec_id, redcap_repeat_instance, drug, start_date, mic50, smg) %>%
    pivot_wider(names_from = "redcap_repeat_instance", 
                values_from = c("drug","start_date","mic50","smg"))

# all_data_merged
all_samples <- full_join(sample_info, seq_info) %>%
    full_join(mic_info) %>%
    filter(mec_id != "MEC103") #exclude one known mixed culture; we have 103-1 and 103-2

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
    select(mec_id, series_id, series_span, series_count) 

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
    select(mec_id, cluster_id, cluster_span, cluster_count)

all_samples <- all_samples %>%
    left_join(serial_calculations, by = c("mec_id", "series_id")) %>%
    left_join(cluster_calcs, by = c("mec_id", "cluster_id")) %T>%
    write.csv("clin_isolate_data.csv")

# priorities: rare isolates, longer than one day, recurrent or clusters
priority_series <- all_samples %>%
    filter(series_id !="NA") %>%
    #filter(genus %in% c("Candida krusei (Pichia kudriavzevii)",
     #                   "Clavispora (Candida) lusitaniae",
      #                  "Candida utilis (Cyberlindnera jadinii)",
       #                 "Candida kefyr (Kluyveromyces marxianus)") |
    filter(series_span > 3, series_span < 30, series_count >2) %>%
    group_by(series_id) %>%
    slice(n=2) %>%
    ungroup() %>%
    select(genus, series_id, series_span, series_count) %>%
    arrange(genus, desc(series_span))

recurrent <- all_samples %>%
    filter(series_id != "NA" & series_span >30) %>%
    group_by(series_id) %>%
    slice_head(n=1) %>%
    ungroup() %>%
    select(genus, series_id, series_span, series_count) %>%
    arrange(genus, desc(series_span))

# cluster info
clusters <- all_samples %>%
    filter(cluster_id != "NA") %>%
    group_by(cluster_id) %>%
    slice_head(n=1) %>%
    ungroup %>%
    select(genus, cluster_id, cluster_span, cluster_count) %>%
    arrange(genus, desc(cluster_span))

# remaining sequencing
series_to_do <- serial_timespans %>%
    filter(genus %in% c("Candida krusei (Pichia kudriavzevii)",
                        "Clavispora (Candida) lusitaniae",
                        "Candida utilis (Cyberlindnera jadinii)",
                        "Candida kefyr (Kluyveromyces marxianus)") |
               (series_span > 3) & is.na(msi_path_r1)) %>%
    ungroup() %>%
    select(mec_id)

glab_to_do <- all_samples %>%
    filter(genus=="Candida glabrata" & is.na(msi_path_r1)) %>%
    select(mec_id)

to_do <- series_to_do %>%
    full_join(series_to_do) %>%
    full_join(glab_to_do) %>%
    pull(mec_id)

seq_remain <- all_samples %>%
    filter(mec_id %in% to_do & is.na(msi_path_r1)) %>%
    select(mec_id, genus, collection_date) %T>%
    write.csv("to_sequences.csv")

## little plot of samples by species
species_count <- all_samples %>%
    group_by(genus) %>%
    summarize(species_count=n()) %>%
    arrange(desc(species_count))

sp_plot <- ggplot(species_count, aes(x=genus, y=species_count)) +
    geom_col() +
    scale_x_discrete(limits=species_count$genus,
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
                                     size = 14)) +
    theme(axis.title = element_text(color = "black", size = 16))

ggsave("species_plot.svg", sp_plot, width = 7, height = 4, units = "in")

# export all data to csv
all_samples %>% write.csv("redcap_all_samples.csv")

# count number of patients per spp.

inds <- all_samples %>%
    group_by(genus) %>%
    count(series_id) %>%
    mutate(total=sum(n))

# create vectors for venn diagram of cluster-series-na pts.

neither <- all_samples %>%
    filter(is.na(series_id) & is.na(cluster_id) & genus == "Candida albicans") %>%
    pull(mec_id)

series <- all_samples %>%
    filter(series_id != "NA" & is.na(cluster_id) & genus == "Candida albicans") %>%
    pull(mec_id)

cluster <- all_samples %>%
    filter(is.na(series_id) & cluster_id != "NA" & genus == "Candida albicans") %>%
    pull(mec_id)

both <- all_samples %>%
    filter(series_id != "NA" & cluster_id != "NA" & genus == "Candida albicans") %>%
    pull(mec_id)

overlaps <- all_samples %>%
    mutate(series_sample = (series_id != "NA" & is.na(cluster_id))) %>%
    mutate(cluster_sample = (cluster_id != "NA" & is.na(series_id))) %>%
    mutate(both_sample = (cluster_id != "NA" & series_id != "NA")) %>%
    mutate(neither_sample = is.na(cluster_id) & is.na(series_id)) %>%
    select(mec_id, genus, series_id, cluster_id, series_sample, cluster_sample, both_sample, neither_sample) #%>%
    #replace_na(list(series_sample=FALSE, cluster_sample=FALSE, both_sample=FALSE, neither_sample=FALSE))

# god-awful effort to get stacked bar plot per spp. of cluster and series set operations
try_again <- overlaps %>%
group_by(genus) %>%
    summarize(neither_sum=sum(neither_sample), both_sum=sum(both_sample),
              cluster_sum=sum(cluster_sample),series_sum=sum(series_sample))
              #total=sum(series_sample, cluster_sample, both_sample, neither_sample))
library(reshape2)
library(viridis)
library(RColorBrewer)

another_try <- melt(try_again, id.var="genus")
ggplot(another_try, aes(x=genus, y=value, fill=variable)) + 
    geom_bar(stat="identity") +
    scale_fill_viridis(discrete=TRUE, option = "D")
    
ggplot(overlaps, aes(interaction(series_sample, cluster_sample), 
                     fill=series_id)) + 
    geom_bar() + facet_grid(. ~ genus) #+
    #scale_fill_viridis(discrete = TRUE)

