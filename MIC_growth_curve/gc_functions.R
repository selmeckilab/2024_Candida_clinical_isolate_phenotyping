## ---------------------------
## Purpose: use growthcurver for gc metric calculation, modify Mew's custom plotting script
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
## Notes:
## ---------------------------
options(scipen = 6, digits = 4) # To view outputs in non-scientific notation
## ---------------------------
## load packages
library(readxl)
library(tidyverse)
library(growthcurver)
library(lubridate)
library(gtools)
################################################################################
## Functions
# Convert excel time to hours and leave in wide format for GrowthCurver functions
clean_growthcurver <- function(plate_reader_file) {
  read_excel(plate_reader_file, sheet = 1) %>%
    mutate(time = hour(Time) * 60 + minute(Time)) %>%
    mutate(time = case_when(
      date(Time) == "1899-12-31" ~ time,
      date(Time) == "1900-01-01" ~ time + 1440,
      date(Time) == "1900-01-02" ~ time + 2880,
      date(Time) == "1900-01-03" ~ time + 4320
    )) %>%
    dplyr::select(!(Time)) %>%
    mutate(time = time / 60)
}

# Read in plate data, wrangle correct time to hours for plotting, and tidy/pivot longer
growth_curve <- function(plate_reader_file) {
  read_excel(plate_reader_file, sheet = 1) %>%
    mutate(time = hour(Time) * 60 + minute(Time)) %>%
    mutate(time = case_when(
      date(Time) == "1899-12-31" ~ time,
      date(Time) == "1900-01-01" ~ time + 1440,
      date(Time) == "1900-01-02" ~ time + 2880
    )) %>%
    dplyr::select(!(Time)) %>%
    mutate(time = time / 60) %>%
    pivot_longer(!time, names_to = "well", values_to = "OD") # keep time value per each well
}

# Read in metadata
samples <- function(plate_reader_file, col_types = NULL) {
  read_excel(plate_reader_file, col_types, sheet = 2)
  # ID <- plate_reader_file$sample
}

# Join tidied plate and metadata
link_metadata <- function(growth_curve, samples, plate) {
  growth_curve %>%
    inner_join(samples, by = "well") %>%
    add_column(plate)
}

# Use SummarizeGrowth function per sample/condition to fit logistic curve, retrieve relevant vals
gc_summary_vals <- function(final_gc) {
  final_gc %>%
    mutate(condition = as.factor(condition)) %>%
    filter(sample != ("blank")) %>%
    group_by(well, condition) %>% # time and OD will be passed as vectors per each well
    mutate(
      k = SummarizeGrowth(time, OD)$vals$k, # carrying capacity
      r = SummarizeGrowth(time, OD)$vals$r, # growth rate
      t_mid = SummarizeGrowth(time, OD)$vals$t_mid, # time at half of carrying capacity
      t_gen = SummarizeGrowth(time, OD)$vals$t_gen, # max doubling time
      auc_e = SummarizeGrowth(time, OD)$vals$auc_e,
      auc_l = SummarizeGrowth(time, OD)$vals$auc_l
    ) # auc
}
