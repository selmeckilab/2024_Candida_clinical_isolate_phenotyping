## ---------------------------
## Purpose: MEC isolate distributions by species and category
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
## Notes:
## ---------------------------
## Load packages
library(tidyverse)
library(magrittr)
library(reshape2)
library(patchwork)

# Redcap report ID
samples <- '58043'

# API token
api_token <- ""

species_colors <- c("#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
                    "#882255","#BBBBBB", "#AA4499", "#DDCC77", "black")

# Function to import report from redcap
import_report <- function(report_number) {
    url <- "https://redcap.ahc.umn.edu/redcap/api/"
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
    filter(isolate_type == "clinical") 

sample_info$genus_species <- str_replace(sample_info$genus_species, "Candida", "C.")

# Calculate time spans and isolate numbers for each series or cluster
serial_timespans <- sample_info %>%
    filter(series_id !="NA") %>%
    group_by(series_id) %>%
    mutate(series_span = max(relative_days) +1) %>% 
    mutate(series_count = n())

serial_summary <- serial_timespans %>%
    select(primary_id, genus_species, series_id, series_span, series_count) 

count_cluster_samples <- sample_info %>%
    filter(cluster_id !="NA") %>% 
    group_by(cluster_id) %>%
    summarise(patients_per_cluster=n_distinct(patient_code), isolates_per_cluster=n()) 

sample_info <- sample_info %>%
    full_join(serial_summary, by = c("primary_id", "series_id", "genus_species")) %>%
    full_join(count_cluster_samples, by = "cluster_id") 

# Data check - single isolates
singles <- sample_info %>% 
    filter(is.na(series_id), multiple_infections=="No", poly_fungal=="No") %>% 
    arrange(patient_code)

# Summary tables
priority_series <- sample_info %>%
    filter(series_id !="NA") %>%
    filter(series_span <= 30, series_count >1) %>%
    group_by(series_id) %>%
    slice(2L) %>%
    ungroup() %>%
    select(genus_species, series_id, series_span, series_count) %>%
    arrange(desc(series_count))

series_tally <- priority_series %>% 
    group_by(genus_species) %>% 
    summarise(serial_sample_totals = sum(series_count)) 


recurrent <- sample_info %>%
    filter(series_id != "NA" & series_span >30) %>%
    group_by(series_id) %>%
    slice(n()) %>%
    ungroup() %>%
    select(genus_species, series_id, series_span, series_count) %>%
    arrange(genus_species, desc(series_span))

polyfungal <- sample_info %>% 
    filter(poly_fungal=="Yes") %>% 
    arrange(patient_code, relative_days) %>% 
    group_by(patient_code, genus_species) %>% 
    summarise(sp_count=n(), collection_span=max(relative_days) + 1)

series_avg <- priority_series %>% 
    group_by(genus_species) %>%
    summarise(num_series =n(),
              mean_series = round(mean(series_count, na.rm =TRUE), digits = 0),
              min_series = min(series_count),
              max_series = max(series_count),
              mean_length=round(mean(series_span, na.rm = TRUE), digits = 1),
              min_length = min(series_span),
              max_length = max(series_span),)

species_count <- sample_info %>%
    group_by(genus_species) %>%
    summarize(species_count=n(), patients=length(unique(patient_code))) %>%
    arrange(desc(species_count)) %>% 
    left_join(series_tally, by=join_by("genus_species")) 

species_count <- species_count %>% 
    group_by(genus_species) %>% 
    mutate(percentage_series=round(serial_sample_totals/species_count, digits = 3))

species_colors <- species_colors %>% 
    set_names(species_count$genus_species)

# Keep species sorting consistent
sorted_species <- which(species_count$genus_species %in% priority_series$genus_species)

################################################################################
# Plot of samples by species

# x-axis commented out for combining with pt plot
sp_plot <- ggplot(species_count, aes(x=genus_species, y=species_count)) +
    geom_col(fill=species_colors, colour = "grey26") +
    scale_x_discrete(limits=species_count$genus_species) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100))+
    ylab("Number of Isolates") +
    xlab(NULL) +
    theme_bw() +
    theme(axis.ticks = element_blank())+
    #theme(axis.text.x = element_text(angle = 35, 
                                    # hjust = 1, 
                                     #vjust = 1,
                                     #color = "black",
                                     #size = 11,
                                    # face = "italic")) +
    theme(axis.text.x = element_blank())+
    theme(axis.title = element_text(color = "black", size = 11))

#ggsave("images/2022_Candida_MEC_spp_distribution.png", sp_plot, 
       #device = png, dpi=300, bg="white",
       #width = 6, height = 4, units = "in")

# Plot of patients by species
pt_plot <- ggplot(species_count, aes(x=genus_species, y=patients)) +
    geom_col(fill=species_colors, colour = "grey26") +
    scale_x_discrete(limits=species_count$genus_species) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100))+
    ylab("Number of Patients") +
    xlab(NULL) +
    theme_bw() +
    theme(axis.ticks = element_blank())+
    theme(axis.text.x = element_text(angle = 35, 
                                     hjust = 1, 
                                     vjust = 1,
                                     color = "black",
                                     size = 11,
                                     face = "italic")) +
    theme(axis.title = element_text(color = "black", size = 11))

#ggsave("images/2022_Candida_MEC_spp_pt_count.png", pt_plot, 
#       device = png, dpi=300, bg="white",
#       width = 6, height = 4, units = "in")

# Combine sp and pt in single figure
both_counts_plot <- sp_plot/pt_plot + plot_annotation(tag_levels = "A")

ggsave("images/2022_Candida_MEC_combined_counts.png", both_counts_plot,
       device = png, dpi = 300, bg="white")

# Stacked bar plot summing patient series
series_plot <- priority_series %>% 
    ggplot( aes(x = genus_species, y = series_count, fill=genus_species))+
    scale_fill_manual(values = species_colors)+
    scale_x_discrete(limits=species_count$genus_species[sorted_species]) +
    scale_y_continuous(limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100))+
    geom_col(colour = "grey26")+ 
    ylab("Number of isolates, grouped by patient") +
    xlab(NULL) +
    theme_bw() +
    theme(axis.ticks = element_blank())+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 35, 
                                     hjust = 1, 
                                     vjust = 1,
                                     color = "black",
                                     size = 11,
                                     face = "italic")) +
    theme(axis.title = element_text(color = "black", size = 11))

ggsave("images/2022_Candida_MEC_series_totals.png", series_plot, 
       device = png, dpi=300, bg="white",
       width = 6, height = 4, units = "in")


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
    theme_bw() +
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
       width = 6.1, height = 4, units = "in")