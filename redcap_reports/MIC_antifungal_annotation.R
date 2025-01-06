## ---------------------------
## Purpose: Join de-identified antifungal history with sample info
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

source("MIC_data_summary.R")
# Local antifungal data
af_spreadsheet <- "~/umn/data/metadata/MEC_all_isolates_antifungal_history.xlsx"

# Process de-identified antifungal relative timeframe data
antifungals <- read_xlsx(af_spreadsheet)

remote_exposure <- antifungals %>%
  filter(relative_start_days <= 0, relative_collection_day == 0)

no_results <- sample_info %>%
  filter(!(primary_id %in% antifungals$primary_id)) %>%
  select(primary_id, genus_species, relative_days, patient_code, series_id, secondary_id)
