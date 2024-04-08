## ---------------------------
## Purpose: Delete specified records by API
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
library(jsonlite)

# Redcap report IDs
samples <- '58043'
genes <- '58048'

# API
api_url <- "https://redcap.ahc.umn.edu/redcap/api/"

# Function to import report from redcap
import_report <- function(report_number) {
    url <- api_url 
    formData <- list("token"=Sys.getenv("redcap_api_key"),
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

# Get report of interest
gene_vars <- import_report(genes) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gene, protein_change, alt_freq) %>%
    left_join((sample_info %>% select(primary_id, genus_species)))

# Filter to what needs to be removed
glab_problem_variants <- gene_vars %>% 
    filter(genus_species == "Candida glabrata", gene == "HOS2")

for(i in 2:length(glab_problem_variants$primary_id)){
  # remove repeat instances
  formData <- list("token"=Sys.getenv("redcap_api_key"),
                   action='delete',
                   content='record',
                   'records[0]'= glab_problem_variants$primary_id[i],
                   instrument='variants_of_interest',
                   repeat_instance = as.character(glab_problem_variants$redcap_repeat_instance[i]),
                   returnContent='count'
  )
  response <- httr::POST(api_url, body = formData, encode = "form")
  #result <- httr::content(response, as="text")
  #print(result)
}
