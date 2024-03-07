## ---------------------------
## Purpose: Upload gene variant records by API
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
## Notes:
## ---------------------------
## Load packages
library(tidyverse)
library(jsonlite)

# Spreadsheet of variants
vars <- read_csv("Cglabrata_MEC_bwa_filtered_annotated.csv", show_col_types = FALSE)
vars$POS <- as.character(vars$POS)
vars <- vars %>% 
    mutate(across(where(is.numeric), ~num(., digits = 2)))
    

amr_genes <- read_csv("redcap_gene_ids.csv")
ref <- "CBS138_s05m03r02"

# API info
api_token <- ""
api_url <-  "https://redcap.ahc.umn.edu/redcap/api/"

# Redcap report IDs
sample_report <- '58043'
gene_report <- '58048'

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

# General sample info 
sample_info <- import_report(sample_report) %>%
    select(-c(starts_with('redcap_repeat'))) %>%
    filter(!primary_id %in% c("MEC103", "MEC113")) %>%
    filter(isolate_type == "clinical")

# Variants in redcap
gene_vars <- import_report(gene_report) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gene, protein_change, alt_freq)

## New data processing
# Add redcap gene codes
vars <- vars %>% 
    left_join(amr_genes, by=join_by(GENE==gene)) %>% 
    filter(!(is.na(redcap_code))) %>% 
    filter(Impact != "LOW")

# Extract sample IDs
samples <- colnames(vars)[str_detect(colnames(vars), "AMS|MEC")] %>%
  str_extract("AMS\\d+|MEC\\d+") %>% unique()

# Non-sample colnames
gene_data <- colnames(vars)[str_detect(colnames(vars), "AMS|MEC", negate = TRUE)]

# Data upload 
# For each strain, create new record and send form
for(j in 1:length(samples)){
  single_sample <- vars[,c(gene_data, paste0(samples[j], ".VAF"), paste0(samples[j], ".GT"))] %>% 
    filter(get(paste0(samples[j],".VAF")) > 0.4)

  sample_name <- samples[j]

# primary vs secondary name
  if(sample_name %in% sample_info$secondary_id){
    primary_id <- as.character(sample_info %>% 
                                   filter(secondary_id==sample_name) %>% 
                                   select(primary_id))
  }else{
    primary_id <- sample_name
  }


  for(i in 1:length(single_sample$POS)){
    
    record <- c(
        primary_id = primary_id,
        redcap_repeat_instrument = "variants_of_interest",
        redcap_repeat_instance = "new",
        gene = single_sample$redcap_code[i],
        snp_ref_genome = ref,
        locus_tag = single_sample$LOCUS[i],
        snp_chromosome = single_sample$CHROM[i],
        snp_pos = single_sample$POS[i],
        snp_ref = single_sample$REF[i],
        snp_alt = single_sample$ALT[i],
        alt_freq = substr(as.character(single_sample[i, paste0(sample_name,".VAF")]),1,4),
        predicted_change = single_sample$Impact[i],
        coding_change = single_sample$Coding_Change[i],
        protein_change = single_sample$AA_change[i],
        alt_gt = as.character(single_sample[i, paste0(sample_name,".GT")]),
        variants_of_interest_complete = 1
       
    )
    
    result_data <- toJSON(list(as.list(record)), auto_unbox=TRUE)
    
    formData <- list("token"=api_token,
                     content='record',
                     action='import',
                     format='json',
                     type='flat',
                     overwriteBehavior='normal',
                     forceAutoNumber='false',
                     data=result_data,
                     returnContent='count',
                     returnFormat='json'
    )
    response <- httr::POST(api_url, body = formData, encode = "form")
    result <- httr::content(response)
    print(result)
  }
}
