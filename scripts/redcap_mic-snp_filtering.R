## ---------------------------
## Purpose: Search for AMR SNPs in resistant isolates 
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

## Load packages
library(tidyverse)
library(readxl)
library(writexl)
library(patchwork)

# Local antifungal data
af_spreadsheet <- "data/metadata/MEC_antifungal_history.xlsx"
breakpoint_file <- "data/metadata/Candida_eucast_breakpoints.xlsx"

# Redcap report IDs
samples <- '58043'
mic_results <- '58044'
growth_curves <- '58045'
genes <- '58048'
rpmi_od <- '64098'

token <- ''

# Function to import report from redcap
import_report <- function(report_number) {
  url <- "https://redcap.ahc.umn.edu/redcap/api/"
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
################################################################################
# Read in reports.

# Sample ID, species, series and cluster IDs
sample_info <- import_report(samples) %>%
    select(-c(starts_with('redcap_repeat'))) %>%
    filter(isolate_type == "clinical")

sample_info$genus_species <- str_replace(sample_info$genus_species, "Candida", "C.")

# OD data
od_vals <- import_report(rpmi_od) %>%
    select(-c(starts_with('redcap_repeat'))) 

# SNP data
gene_vars <- import_report(genes) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gene, protein_change, alt_freq) %>%
    left_join((sample_info %>% select(primary_id, genus_species))) %>% 
    unite("gene_pos", gene:protein_change)

# MIC and SMG results
mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, drug, mic_media,mic_date, mic50, eucast_breakpoint, smg, qc_ok) %>% 
    left_join((sample_info %>% 
                   select(primary_id, genus_species, series_id)), 
              by=join_by(primary_id))

# Subset to RPMI data with valid control results
# Slice head removes repeated assays (for now, filtering isn't working)
mic_info <- mic_info %>% 
    filter(mic_media=="RPMI", !(primary_id %in% c("AMS5122","AMS5123")), qc_ok=="Yes") %>% 
    inner_join(sample_info %>% select(primary_id, patient_code)) %>% 
    group_by(primary_id, drug) %>% 
    slice_head()

mic_info <- mic_info %>% left_join(od_vals, by=c("primary_id", "mic_date"="rpmi_date"))

# For ordering species and drug levels
species_count <- sample_info %>%
    group_by(genus_species) %>%
    summarize(species_count=n()) %>%
    arrange(desc(species_count))

mic_info$mic50 <- factor(mic_info$mic50, levels=c("0.016", "0.032", "0.064", "0.128",
                                                  "0.256", "0.5", "1", ">1", "2", 
                                                  "4", "8", "16", "32", ">32"))

mic_info$genus_species <- factor(mic_info$genus_species, levels = species_count$genus_species)

################################################################################
# Drug exposure data

# Process de-identified antifungal relative timeframe data
antifungals <- read_xlsx(af_spreadsheet)

remote_exposure <- antifungals %>% 
    filter(relative_start_days <= 0, relative_collection_day == 0)

no_results <- sample_info %>% filter(!(primary_id %in% antifungals$primary_id)) %>% 
    select(primary_id, genus_species, relative_days, patient_code, series_id, secondary_id)

################################################################################
# Summary tables for resistant isolates

# All isolates done for each drug?
isolate_counts <- mic_info %>%
    group_by(drug) %>%
    count(genus_species)

# All resistant isolates
resistant_isolates <- mic_info %>% 
    filter(eucast_breakpoint=="R")

resistant_series <- mic_info %>% 
    group_by(patient_code, drug) %>% 
    filter(any(eucast_breakpoint=="R"))

patient_counts <- mic_info %>%
    group_by(drug, genus_species) %>%
    summarize(patients=length(unique(patient_code)))

patient_resistance <- resistant_isolates %>% 
    group_by(drug, genus_species) %>%
    summarize(patient_res=length(unique(patient_code)))

patient_differences <- resistant_series %>% 
    group_by(patient_code, series_id, drug, eucast_breakpoint) %>%
    summarise(patient_diffs = n())

# Frequency of resistant isolates among patients
patient_freq <- patient_counts %>% 
    inner_join(patient_resistance, by=join_by(drug, genus_species)) %>% 
    mutate(patient_res_freq = patient_res/patients) %>% 
    arrange(genus_species, drug)

################################################################################
# Candidate AMR gene SNPs
# Variants found in sensitive isolates 

non_driver_vars <- gene_vars %>% 
    filter(!primary_id %in% resistant_isolates$primary_id) %>% 
    left_join(sample_info %>% select(primary_id, series_id))

# Variants found only in resistant isolates
possible_drivers <- gene_vars %>% 
    filter(!gene_pos %in% non_driver_vars$gene_pos) %>% 
    left_join(sample_info %>% select(primary_id, series_id))

################################################################################
# How does low OD in no-drug control relate to resistance status?
