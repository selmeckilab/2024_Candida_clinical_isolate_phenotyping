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
library(reshape2)

# redcap report IDs
samples <- '41047'
mic_results <- '41628'
growth_curves <- '53245'
chef_data <- '52263'
spot_plates <- '56394'
genes <- '56393'
avail_seq_data <- '52998'

# function to import report from redcap
import_report <- function(report_number) {
  token <- token # API token (string) here; I've got a placeholder bc of github
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
    select(primary_id, redcap_repeat_instance, gene, protein_change, alt_freq) %>%
    pivot_wider(names_from = "redcap_repeat_instance",
                values_from = c("gene", "protein_change", "alt_freq"),
                names_vary = "slowest")

# CHEF gel results
chef_done <- import_report(chef_data) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gel_date, blot_prepared) %>%
    pivot_wider(names_from = "redcap_repeat_instance",
                values_from = c("gel_date", "blot_prepared"),
                names_vary = "slowest")

gc <- import_report(growth_curves) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gc_date, drug_used, gc_temp, gc_time, k, t_gen, auc_l) %>%
    pivot_wider(names_from = "redcap_repeat_instance",
                values_from = c("gc_date", "drug_used", "gc_temp", "gc_time", "k", "t_gen", "auc_l"),
                names_vary = "slowest")

# No CHEF results yet
todo <- sample_info %>%
    filter(isolate_type == "clinical") %>%
    anti_join(chef_done)

# export chef-todo
#write_xlsx(todo,paste0(Sys.Date(),"_CHEF_todo.xlsx"))

# MSI location if sequencing data exists
seq_info <- import_report(avail_seq_data) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, msi_path_r1, msi_path_r2, msi_nanopore_path) %>%
    pivot_wider(names_from = "redcap_repeat_instance",
                values_from = c("msi_path_r1", "msi_path_r2", "msi_nanopore_path"),
                names_vary = "slowest")

# MIC and SMG results
mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, drug, mic_date, mic50, smg) %>%
    pivot_wider(names_from = "redcap_repeat_instance", 
                values_from = c("drug","mic50","smg", "mic_date"), 
                names_vary = "slowest")

# all_data_merged
all_samples <- left_join(sample_info, gene_vars) %>%
   left_join(mic_info, by = join_by("primary_id")) %>%
   left_join(chef_done, by = join_by("primary_id")) %>%
   left_join(gc, by = join_by("primary_id")) %>%
   filter(!primary_id %in% c("MEC103", "MEC113")) 

write_xlsx(all_samples,paste0(Sys.Date(),"merged_Candida_data.xlsx"))

