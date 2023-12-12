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
library(writexl)

# redcap report IDs
samples <- '58043'
mic_results <- '58044'
growth_curves <- '58045'
chef_data <- '58046'
spot_plates <- '58047'
genes <- '58048'
avail_seq_data <- '58050'

token <- "" # don't forget to delete before gh

# function to import report from redcap
import_report <- function(report_number) {
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
    select(-c(starts_with('redcap_repeat'))) %>%
    filter(!primary_id %in% c("MEC103", "MEC113")) %>%
    filter(isolate_type == "clinical")

# CHEF gel results
chef_done <- import_report(chef_data) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gel_date, blot_prepared) %>%
    pivot_wider(names_from = "redcap_repeat_instance",
                values_from = c("gel_date", "blot_prepared"),
                names_vary = "slowest")

# No CHEF results yet
todo <- sample_info %>%
    filter(isolate_type == "clinical") %>%
    anti_join(chef_done)

# export chef-todo
#write_xlsx(todo,paste0(Sys.Date(),"_CHEF_todo.xlsx"))

# MIC and SMG results
mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, drug, mic_date, mic50, eucast_breakpoint, smg) %>% 
    left_join((sample_info %>% select(primary_id, genus_species, series_id)), by=join_by(primary_id)) 

# MIC to-do
mic_to_do <- sample_info %>% 
    filter(isolate_type == "clinical") %>%
    anti_join(mic_info)

# growth curve results
gc <- import_report(growth_curves) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gc_date, drug_used, gc_temp, gc_time, k, r, t_gen, auc_l)

# MSI location if sequencing data exists
seq_info <- import_report(avail_seq_data) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, msi_path_r1, msi_path_r2, msi_long_read_path) %>%
    pivot_wider(names_from = "redcap_repeat_instance",
                values_from = c("msi_path_r1", "msi_path_r2", "msi_long_read_path"),
                names_vary = "slowest")

# snp data
gene_vars <- import_report(genes) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gene, protein_change, alt_freq) %>%
    left_join((sample_info %>% select(primary_id, genus_species)))
