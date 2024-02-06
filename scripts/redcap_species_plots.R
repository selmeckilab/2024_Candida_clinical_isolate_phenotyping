## ---------------------------
## Script name: redcap_species_plots.R
##
## Purpose of script: MEC isolate distributions by species and category
##
## Author: Nancy Scott
##
## Date Created: 2023-11-28
##
## Email: scot0854@umn.edu
## ---------------------------
## Notes:
## ---------------------------
## load packages
library(tidyverse)
library(lubridate)
library(magrittr)
library(reshape2)

# redcap report ID
samples <- '58043'

# api token
api_token <- ""

# function to import report from redcap
import_report <- function(report_number) {
    url <- "https://redcap.ahc.umn.edu/api/"
  formData <- list("token"=api_token,
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

# Pull sample data
sample_info <- import_report(samples) %>%
    select(-c(starts_with('redcap_repeat'))) %>% 
    filter(primary_id != "MEC103") %>%
    filter(isolate_type == "clinical")

# Calculate timeframes and isolate numbers for each series or cluster
serial_timespans <- sample_info %>%
    filter(series_id !="NA") %>%
    group_by(series_id) %>%
    mutate(series_span = max(relative_days))

count_serial_samples <- serial_timespans %>%
    group_by(series_id) %>%
    filter(series_span >0) %>%
    summarize(series_count = n())

serial_timespans <- serial_timespans %>%
    inner_join(count_serial_samples, by = 'series_id')
    
serial_calculations <- serial_timespans %>%
    select(primary_id, series_id, series_span, series_count) 

count_cluster_samples <- sample_info %>%
    filter(cluster_id !="NA") %>% 
    group_by(cluster_id) %>%
    summarize(cluster_count=n()) 

sample_info <- sample_info %>%
    full_join(serial_calculations, by = c("primary_id", "series_id")) %>%
    full_join(count_cluster_samples, by = "cluster_id") 

# Tables of priorities: longer spans/more samples, recurrent infection or clusters
priority_series <- sample_info %>%
    filter(series_id !="NA") %>%
    filter(series_span > 1, series_span < 30, series_count >2) %>%
    group_by(series_id) %>%
    slice(2L) %>%
    ungroup() %>%
    select(genus_species, series_id, series_span, series_count) %>%
    arrange(desc(series_span))

recurrent <- sample_info %>%
    filter(series_id != "NA" & series_span >30) %>%
    group_by(series_id) %>%
    slice(2L) %>%
    ungroup() %>%
    select(genus_species, series_id, series_span, series_count) %>%
    arrange(genus_species, desc(series_span))

## Plot of samples by species
species_count <- sample_info %>%
    group_by(genus_species) %>%
    summarize(species_count=n(), patients=length(unique(patient_code))) %>%
    arrange(desc(species_count))

species_colors <- c("#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
                     "#882255","#BBBBBB", "#AA4499", "#DDCC77", "black")

species_colors <- species_colors %>% 
    set_names(species_count$genus_species)

sp_plot <- ggplot(species_count, aes(x=genus_species, y=species_count)) +
    geom_col(fill=species_colors) +
    scale_x_discrete(limits=species_count$genus_species) +
    ylab("Number of Isolates") +
    xlab("Species") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 35, 
                                     hjust = 1, 
                                     vjust = 1,
                                     color = "black",
                                     size = 10,
                                     face = "italic")) +
    theme(axis.title = element_text(color = "black", size = 16))

ggsave("images/2022_Candida_MEC_spp_distribution.png", sp_plot, 
       device = png, dpi=300, bg="white",
       width = 7, height = 5, units = "in")

pt_plot <- ggplot(species_count, aes(x=genus_species, y=patients)) +
    geom_col(fill=species_colors) +
    scale_x_discrete(limits=species_count$genus_species) +
    scale_y_continuous(limits = c(0,100))+
    ylab("Number of Patients") +
    xlab("Species") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 35, 
                                     hjust = 1, 
                                     vjust = 1,
                                     color = "black",
                                     size = 10,
                                     face = "italic")) +
    theme(axis.title = element_text(color = "black", size = 16))

ggsave("images/2022_Candida_MEC_spp_pt_count.png", pt_plot, 
       device = png, dpi=300, bg="white",
       width = 7, height = 5, units = "in")

# Create vectors for plot comparing cluster-series-na counts
category_counts <- sample_info %>%
    mutate(series_sample = (series_id != "NA" & is.na(cluster_id))) %>%
    mutate(cluster_sample = (cluster_id != "NA" & is.na(series_id))) %>%
    mutate(both_sample = (cluster_id != "NA" & series_id != "NA")) %>%
    mutate(neither_sample = is.na(cluster_id) & is.na(series_id)) %>%
    select(primary_id, genus_species, series_id, cluster_id, series_sample, cluster_sample, both_sample, neither_sample) %>%
    replace_na(list(series_sample=FALSE, cluster_sample=FALSE, both_sample=FALSE, neither_sample=FALSE)) %>% 
    group_by(genus_species) %>%
    summarize(neither_sum=sum(neither_sample==TRUE), 
              both_sum=sum(both_sample==TRUE),
              cluster_sum=sum(cluster_sample==TRUE),
              series_sum=sum(series_sample==TRUE), 
              species_sum = sum(neither_sum, both_sum, cluster_sum, series_sum)) %>%
    arrange(desc(species_sum)) %>% # use species counts to force order as above
    select(!(species_sum)) %>% 
    pivot_longer(names_to = "Category", values_to = "counts", -c(genus_species))
    
# generate stacked bar plot per spp. of cluster, series and individuals
# dataset is reshaped from wide to long using melt
viz_categories <- ggplot(category_counts, aes(x=genus_species, 
                                              y=counts, 
                                              fill=factor(Category, levels = c("neither_sum", "cluster_sum", "both_sum", "series_sum")))) + 
    geom_bar(stat="identity") +
    scale_x_discrete(limits=unique(category_counts$genus_species)) +
    scale_fill_manual(values = c("#EECC66", "#EE99AA", "#6699CC", "#004488"),
                      labels = c("Single", "Demographic Cluster", "Cluster + series", "Series")) +
    xlab("Species") +
    ylab("Count of isolates") +
    labs(fill="Category") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 35, 
                                     hjust = 1, 
                                     vjust = 1,
                                     color = "black",
                                     size = 10,
                                     face = "italic")) +
    theme(axis.title = element_text(color = "black", size = 16)) +
    theme(legend.text = element_text(size = 8)) +
    theme(legend.box.spacing = unit(0, "pt"))
    
ggsave("images/2022_Candida_patient_categories.png", viz_categories,
       device = png, dpi=300, bg="white",
       width = 7.2, height = 5, units = "in")
