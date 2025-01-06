## ---------------------------
## Purpose: Get SRA data and sequencing paths
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
## Load packages
library(tidyverse)
library(writexl)

# Redcap report IDs
samples <- "58043"
avail_seq_data <- "58050"
sra_data <- "63448"

# Function to import report from redcap
import_report <- function(report_number) {
  url <- "https://redcap.ahc.umn.edu/redcap/api/"
  formData <- list(
    "token" = Sys.getenv("redcap_api_key"),
    content = "report",
    format = "csv",
    report_id = report_number,
    csvDelimiter = "",
    rawOrLabel = "label",
    rawOrLabelHeaders = "raw",
    exportCheckboxLabel = "true",
    returnFormat = "csv"
  )
  response <- httr::POST(url, body = formData, encode = "form")
  result <- httr::content(response, show_col_types = FALSE)
}

# Sample ID, species, series and cluster IDs
sample_info <- import_report(samples) %>%
  select(-c(starts_with("redcap_repeat"))) %>%
  filter(isolate_type == "clinical")

# SRA data
sra <- import_report(sra_data) %>%
  filter(redcap_repeat_instrument != "NA") %>%
  pivot_wider(
    names_from = "redcap_repeat_instance",
    values_from = c("bioproject", "biosample", "sra_1", "sra_2", "sra_3"),
    names_vary = "slowest"
  )

# MSI location if sequencing data exists
seq_info <- import_report(avail_seq_data) %>%
  filter(redcap_repeat_instrument != "NA") %>%
  select(primary_id, redcap_repeat_instance, msi_path_r1, msi_path_r2, msi_long_read_path) %>%
  pivot_wider(
    names_from = "redcap_repeat_instance",
    values_from = c("msi_path_r1", "msi_path_r2", "msi_long_read_path"),
    names_vary = "slowest"
  )
