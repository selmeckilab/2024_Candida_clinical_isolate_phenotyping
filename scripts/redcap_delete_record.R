## ---------------------------
## Script name: redcap_delete_record.R
##
## Purpose of script: Delete specified records by API
##
## Author: Nancy Scott
##
## Date Created: 2023-12-04
##
## Email: scot0854@umn.edu

library(jsonlite)

# redcap report IDs
samples <- '58043'
genes <- '58048'

# api
token <- ""
url <- "https://redcap.ahc.umn.edu/api/"

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

# get report of interest
gene_vars <- import_report(genes) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gene, protein_change, alt_freq) %>%
    left_join((sample_info %>% select(primary_id, genus_species)))

# filter to what needs to be removed
glab_problem_variants <- gene_vars %>% 
    filter(genus_species == "Candida glabrata", gene == "HOS2")

for(i in 2:length(glab_problem_variants$primary_id)){
  # remove repeat instances
  formData <- list("token"=token,
                   action='delete',
                   content='record',
                   'records[0]'= glab_problem_variants$primary_id[i],
                   instrument='variants_of_interest',
                   repeat_instance = as.character(glab_problem_variants$redcap_repeat_instance[i]),
                   returnContent='count'
  )
  response <- httr::POST(url, body = formData, encode = "form")
  #result <- httr::content(response, as="text")
  #print(result)
}
