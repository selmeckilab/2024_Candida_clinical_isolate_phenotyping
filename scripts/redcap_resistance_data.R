## ---------------------------
## Script name: redcap_resistance_data.R
##
## Purpose of script: Pull REDCap reports to analyze MIC and gene variants
##
## Author: Nancy Scott
##
## Date Created: 2023-12-28
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
genes <- '58048'
calb_mlst <- '58053'
cglab_mlst <- '58052'

token <- '' 

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

# MIC and SMG results
mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, drug, mic_date, mic50, eucast_breakpoint, smg) %>% 
    left_join((sample_info %>% select(primary_id, genus_species, series_id)), by=join_by(primary_id)) 

# MIC to-do
mic_to_do <- sample_info %>% 
    filter(isolate_type == "clinical") %>%
    anti_join(mic_info)

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
    mutate(concat_alleles=paste(aat1a_exact_match, acc1_exact_match, adp1_exact_match, mpib_exact_match, sya1_exact_match, vps13_exact_match, zwf1b_exact_match, sep = ""))

flc_smg <- mic_info %>% filter(drug=="fluconazole", smg > 0.3)

mcf_res <- mic_info %>% filter(drug=="micafungin", eucast_breakpoint=="R")
mcf_sensitive <- mic_info %>% filter(drug=="micafungin", eucast_breakpoint=="S", genus_species=="Candida glabrata")
glabrata_vars <- gene_vars %>% filter(genus_species=="Candida glabrata")
glab_mcf <- mic_info %>% filter(genus_species=="Candida glabrata", drug=="micafungin")
mcf_res_vars <- glabrata_vars %>% filter(primary_id %in% mcf_res$primary_id) %>% select(gene, protein_change)
mcf_sens_vars <- glabrata_vars %>% filter(primary_id %in% mcf_sensitive$primary_id) %>% select(gene, protein_change)
possible_mcf <- mcf_res_vars %>% anti_join(mcf_sens_vars)


glab_flc <- mic_info %>% filter(genus_species=="Candida glabrata", drug=="fluconazole", !(is.na(eucast_breakpoint)))
flc_res <- glab_flc %>% filter(eucast_breakpoint=="R")
flc_int <- glab_flc %>% filter(eucast_breakpoint=="I")
flc_sens_vars <- glabrata_vars %>% filter(primary_id %in% flc_int$primary_id) %>% select(gene, protein_change)
flc_res_vars <- glabrata_vars %>% filter(primary_id %in% flc_res$primary_id) %>% select(gene, protein_change)
possible_flc <- flc_res_vars %>% anti_join(flc_sens_vars)

possible_mcf <- unique(mcf_res_vars %>% anti_join(mcf_sens_vars))
possible_flc <- unique(flc_res_vars %>% anti_join(flc_sens_vars))

multi_res_vars <- glabrata_vars %>% filter(primary_id %in% mcf_res$primary_id, primary_id %in% flc_res$primary_id) %>% select(gene, protein_change)
possible_multi <- unique(multi_res_vars %>% anti_join(mcf_sens_vars)) %>% anti_join(flc_sens_vars)
