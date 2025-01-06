## ---------------------------
## Purpose: CHEF gel summary
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
## Load packages
library(tidyverse)
library(writexl)

# Redcap report IDs
samples <- "58043"
chef_data <- "58046"

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

# CHEF gel results
chef_done <- import_report(chef_data) %>%
  filter(redcap_repeat_instrument != "NA") %>%
  select(primary_id, redcap_repeat_instance, gel_date, chef_file_name, chef_lane, blot_prepared) %>%
  pivot_wider(
    names_from = "redcap_repeat_instance",
    values_from = c("gel_date", "blot_prepared", "chef_file_name", "chef_lane"),
    names_vary = "slowest"
  )

# No CHEF results yet
todo <- sample_info %>%
  filter(isolate_type == "clinical") %>%
  anti_join(chef_done)

# Export chef-todo list
# write_xlsx(todo,paste0(Sys.Date(),"_CHEF_todo.xlsx"))
