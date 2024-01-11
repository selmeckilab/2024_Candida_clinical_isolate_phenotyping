## ---------------------------
## Script name: redcap_reports.R
##
## Purpose of script: Pull, filter, subset REDCap reports
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
library(writexl)

# redcap report IDs
samples <- '58043'
mic_results <- '58044'
growth_curves <- '58045'
chef_data <- '58046'
spot_plates <- '58047'
genes <- '58048'
avail_seq_data <- '58050'
calb_mlst <- '58053'
cglab_mlst <- '58052'

token <- '' # no gh

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

# mlst types
albicans_sts <- import_report(calb_mlst) %>% 
    select(primary_id, st, aat1a_exact_match, acc1_exact_match, adp1_exact_match,
           mpib_exact_match, sya1_exact_match, vps13_exact_match, zwf1b_exact_match) %>% 
    mutate(across(st:zwf1b_exact_match, as.character)) %>% 
    mutate(concat_alleles=paste(aat1a_exact_match, acc1_exact_match, 
                                adp1_exact_match, mpib_exact_match, 
                                sya1_exact_match, vps13_exact_match, 
                                zwf1b_exact_match, sep = "")) %>% 
    mutate(st = case_when(is.na(st) ~ concat_alleles, .default = st))

