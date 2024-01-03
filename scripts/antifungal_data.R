## ---------------------------
## Script name: antifungal_data.R
##
## Purpose of script: Analyze redcap and antifungal exposure data jointly
##
## Author: Nancy Scott
##
## Date Created: 2024-01-03
##
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999) 

## load packages
library(tidyverse)
library(writexl)

# local antifungal data
af_spreadsheet <- "antifungal_relative_days_by_sample.xlsx"

# redcap report IDs
samples <- '58043'
mic_results <- '58044'
growth_curves <- '58045'
genes <- '58048'

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

# snp data
gene_vars <- import_report(genes) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gene, protein_change, alt_freq) %>%
    left_join((sample_info %>% select(primary_id, genus_species)))

# MIC and SMG results
mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, drug, mic_date, mic50, eucast_breakpoint, smg) %>% 
    left_join((sample_info %>% 
                   select(primary_id, genus_species, series_id)), 
              by=join_by(primary_id))

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


# de-identified antifungal data from ahc
antifungals <- read_xlsx(af_spreadsheet)
antifungals <- antifungals %>% 
    inner_join((sample_info %>% 
                    select(primary_id, genus_species, series_id, relative_days)), 
               by=join_by(sample==primary_id))

remote_exposure <- antifungals %>% 
    group_by(series_id) %>% 
    filter(relative_start < 0, relative_days == 0)

no_results <- sample_info %>% filter(!(primary_id %in% antifungals$sample)) %>% 
    select(primary_id, genus_species, relative_days, series_id, secondary_id)
