## ---------------------------
## Purpose: Clean and summarise MIC and SMG data for candidemia isolates
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

## Load packages
library(tidyverse)
library(readxl)
library(writexl)

# Redcap report IDs
samples <- "58043"
mic_results <- "58044"

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
################################################################################
# Read in reports.

# Sample ID, species, series and cluster IDs
sample_info <- import_report(samples) %>%
  select(-c(starts_with("redcap_repeat"))) %>%
  filter(isolate_type == "clinical")

sample_info$genus_species <- str_replace(sample_info$genus_species, "Candida", "C.")

# MIC and SMG results
mic_info <- import_report(mic_results) %>%
  filter(redcap_repeat_instrument != "NA") %>%
  select(
    primary_id,
    redcap_repeat_instance,
    drug,
    mic_media,
    mic_date,
    mic50,
    mean_mic50_relative_growth,
    mic90,
    mean_mic90_relative_growth,
    eucast_breakpoint,
    mean_smg,
    sem_smg,
    mean_no_drug_24hr_od,
    sd_no_drug_24hr_od,
    qc_ok,
    use_mic_result
  ) %>%
  left_join(
    (sample_info %>%
      select(primary_id, genus_species, series_id, patient_code, relative_days)),
    by = join_by(primary_id)
  )

# Subset to RPMI data with valid control results
mic_info <- mic_info %>%
  filter(
    qc_ok == "Yes",
    !(primary_id %in% c("AMS5122", "AMS5123", "AMS2401")),
    use_mic_result == "Yes",
  ) %>%
  group_by(primary_id, drug) %>%
  arrange(desc(mic_date)) %>%
  slice_head() # newest valid result if there are repeats

# For ordering species and drug levels
species_count <- sample_info %>%
  group_by(genus_species) %>%
  summarize(species_count = n()) %>%
  arrange(desc(species_count))

mic_info$mic50 <- factor(mic_info$mic50, levels = c(
  "0.016", "0.032", "0.047", "0.064", "0.125",
  "0.256", "0.5", "1", ">1", "2",
  "4", "8", "16", "32", ">32", "64",
  "128", "160", "256", ">256"
), ordered = TRUE)

mic_info$mic90 <- factor(mic_info$mic90, levels = c(
  "0.016", "0.023", "0.032", "0.047", "0.064", "0.125",
  "0.256", "0.38", "0.5", "1", ">1", "2",
  "4", "8", "16", "32", ">32", "64",
  "96", "128", "160", "256", ">256"
), ordered = TRUE)

mic_info$genus_species <- factor(mic_info$genus_species, levels = species_count$genus_species)

# No-drug well OD vals
control_od <- mic_info %>%
  group_by(primary_id) %>%
  summarise(mean_stationary_k = mean(mean_no_drug_24hr_od, na.rm = TRUE)) %>%
  left_join(sample_info)

control_od$genus_species <- factor(control_od$genus_species, levels = species_count$genus_species)
control_od_summary <- control_od %>%
  select(genus_species, mean_stationary_k) %>%
  split(.$genus_species) %>%
  map(summary)

################################################################################
# For subsetting SMG to isolates with adequate 24-hour growth
carrying_cap <- mic_info %>%
  group_by(primary_id) %>%
  filter(mean_no_drug_24hr_od < 0.2)

flc_carrying_cap <- carrying_cap %>%
  filter(drug == "fluconazole")

# Summarise all FLC SMG by species
smg <- mic_info %>%
  filter(!primary_id %in% flc_carrying_cap$primary_id, drug == "fluconazole") %>%
  group_by(genus_species, drug) %>%
  filter(!is.na(mean_smg)) %>%
  summarise(
    number_smg = n(),
    min_smg = round(min(mean_smg), digits = 2),
    max_smg = round(max(mean_smg), digits = 2),
    overall_mean_smg = round(mean(mean_smg), digits = 2),
    median_smg = round(median(mean_smg), digits = 2),
    IQR_smg = round(IQR(mean_smg), digits = 2)
  )

# Summarise changes in FLC SMG within each series
serial_flc_smgs <- mic_info %>%
  filter(!primary_id %in% flc_carrying_cap$primary_id, drug == "fluconazole") %>%
  group_by(genus_species, drug, series_id) %>%
  filter(!is.na(series_id), !is.na(mean_smg)) %>%
  summarise(
    number_smg = n(),
    min_smg = round(min(mean_smg), digits = 2),
    max_smg = round(max(mean_smg), digits = 2),
    diff_smg = max(mean_smg) - min(mean_smg)
  )

################################################################################
# Summary tables for resistant isolates

# All isolates done for each drug? No duplicates?
isolate_counts <- mic_info %>%
  group_by(drug) %>%
  count(genus_species) %>%
  select(genus_species, n) %>%
  group_by(genus_species, n) %>%
  count()

# All resistant isolates based on breakpoints
resistant_isolates <- mic_info %>%
  filter(eucast_breakpoint == "R") %>%
  select(
    primary_id,
    drug,
    mic50,
    eucast_breakpoint,
    genus_species,
    series_id,
    patient_code,
    relative_days
  )

resistance_freqs <- mic_info %>%
  group_by(genus_species, drug) %>%
  summarize(
    all_count = length(unique(primary_id)),
    res_count = sum(eucast_breakpoint == "R", na.rm = TRUE),
    res_percent = round((sum(eucast_breakpoint == "R", na.rm = TRUE)) / length(unique(primary_id)), digits = 3) * 100
  ) %>%
  filter(res_percent != 0)

# Differences within series
change_in_res_categories <- mic_info %>%
  group_by(series_id, patient_code, drug) %>%
  filter(!is.na(series_id)) %>%
  filter(any(eucast_breakpoint == "R")) %>%
  summarize(
    totals = length(unique(primary_id)),
    res = sum(eucast_breakpoint == "R", na.rm = TRUE),
    int = sum(eucast_breakpoint == "I", na.rm = TRUE),
    sens = sum(eucast_breakpoint == "S", na.rm = TRUE)
  ) %>%
  filter(!(sens == 0 & int == 0))

serial_mic50 <- mic_info %>%
  filter(!is.na(series_id), drug %in% c("fluconazole", "micafungin")) %>%
  group_by(genus_species, series_id, drug, mic_media, mic50) %>%
  count()

serial_mic90 <- mic_info %>%
  filter(!is.na(series_id), drug %in% c("amphotericin B")) %>%
  group_by(genus_species, series_id, drug, mic_media, mic90) %>%
  count()

changed_amb_series <- serial_mic90 %>%
  filter(drug %in% c("amphotericin B")) %>%
  group_by(genus_species, series_id, mic_media) %>%
  filter(n() > 1)

changed_mic50_series <- serial_mic50 %>%
  filter(drug %in% c("fluconazole", "micafungin")) %>%
  group_by(genus_species, series_id, mic_media) %>%
  filter(n() > 1)

changes_by_drug <- changed_mic50_series %>%
  rbind(changed_amb_series) %>%
  group_by(genus_species, series_id, drug) %>%
  count() %>%
  filter(n > 1) %>%
  arrange(genus_species, series_id)

patient_counts <- mic_info %>%
  group_by(drug, genus_species) %>%
  summarize(patients = length(unique(patient_code)))

patient_resistance <- resistant_isolates %>%
  group_by(drug, genus_species) %>%
  summarize(patient_res = length(unique(patient_code)))

# Frequency of resistant isolates among patients
patient_freq <- patient_counts %>%
  inner_join(patient_resistance, by = join_by(drug, genus_species)) %>%
  mutate(patient_res_freq = round(patient_res / patients, digits = 3) * 100) %>%
  arrange(genus_species, drug)
