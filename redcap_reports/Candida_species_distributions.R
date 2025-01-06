## ---------------------------
## Purpose: MEC isolate distributions by species and category
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
library(tidyverse)
library(magrittr)
library(reshape2)
library(patchwork)

species_colors <- c(
  "#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
  "#882255", "#BBBBBB", "#AA4499", "#DDCC77", "black"
)

# Redcap report ID
samples <- "58043"

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

# Pull sample data
sample_info <- import_report(samples) %>%
  select(-c(starts_with("redcap_repeat"))) %>%
  filter(isolate_type == "clinical")

sample_info$genus_species <- str_replace(sample_info$genus_species, "Candida", "C.")

# Calculate time spans and isolate numbers for each series or cluster
serial_timespans <- sample_info %>%
  filter(series_id != "NA") %>%
  group_by(series_id) %>%
  mutate(series_span = max(relative_days) + 1) %>%
  mutate(series_count = n())

serial_summary <- serial_timespans %>%
  select(primary_id, genus_species, series_id, series_span, series_count)

count_cluster_samples <- sample_info %>%
  filter(cluster_id != "NA") %>%
  group_by(cluster_id) %>%
  summarise(patients_per_cluster = n_distinct(patient_code), isolates_per_cluster = n())

sample_info <- sample_info %>%
  full_join(serial_summary, by = c("primary_id", "series_id", "genus_species")) %>%
  full_join(count_cluster_samples, by = "cluster_id")

# Data check - single isolates
singles <- sample_info %>%
  filter(is.na(series_id), multiple_infections == "No", poly_fungal == "No") %>%
  arrange(patient_code)

# Summary tables
priority_series <- sample_info %>%
  filter(series_id != "NA") %>%
  filter(series_span <= 30, series_count > 1) %>%
  group_by(series_id) %>%
  slice(2L) %>%
  ungroup() %>%
  select(genus_species, series_id, series_span, series_count) %>%
  arrange(desc(series_count))

series_tally <- priority_series %>%
  group_by(genus_species) %>%
  summarise(serial_sample_totals = sum(series_count))


recurrent <- sample_info %>%
  filter(series_id != "NA" & series_span > 30) %>%
  group_by(series_id) %>%
  slice(n()) %>%
  ungroup() %>%
  select(genus_species, series_id, series_span, series_count) %>%
  arrange(genus_species, desc(series_span))

polyfungal <- sample_info %>%
  filter(poly_fungal == "Yes") %>%
  arrange(patient_code, relative_days) %>%
  group_by(patient_code, genus_species) %>%
  summarise(sp_count = n(), collection_span = max(relative_days) + 1)

series_avg <- priority_series %>%
  group_by(genus_species) %>%
  summarise(
    num_series = n(),
    mean_series = round(mean(series_count, na.rm = TRUE), digits = 0),
    min_series = min(series_count),
    max_series = max(series_count),
    mean_length = round(mean(series_span, na.rm = TRUE), digits = 1),
    min_length = min(series_span),
    max_length = max(series_span),
  )

species_count <- sample_info %>%
  group_by(genus_species) %>%
  summarize(species_count = n(), patients = length(unique(patient_code))) %>%
  arrange(desc(species_count)) %>%
  left_join(series_tally, by = join_by("genus_species"))

species_count <- species_count %>%
  group_by(genus_species) %>%
  mutate(percentage_series = round(serial_sample_totals / species_count, digits = 3))

species_colors <- species_colors %>%
  set_names(species_count$genus_species)

# Keep species sorting consistent
sorted_species <- which(species_count$genus_species %in% priority_series$genus_species)

################################################################################
# Plot sample and pt counts per species

sp_plot <- ggplot(species_count, aes(x = genus_species, y = species_count)) +
  geom_col(fill = species_colors, colour = "grey26") +
  scale_x_discrete(limits = species_count$genus_species) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  ylab("Number of Isolates") +
  xlab(NULL) +
  theme_bw() +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title = element_text(color = "black", size = 11))

pt_plot <- ggplot(species_count, aes(x = genus_species, y = patients)) +
  geom_col(fill = species_colors, colour = "grey26") +
  scale_x_discrete(limits = species_count$genus_species) +
  scale_y_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) +
  ylab("Number of Patients") +
  xlab(NULL) +
  theme_bw() +
  theme(axis.ticks = element_blank()) +
  theme(axis.text.x = element_text(
    angle = 35,
    hjust = 1,
    vjust = 1,
    color = "black",
    size = 10,
    face = "italic"
  )) +
  theme(axis.title = element_text(color = "black", size = 11))

both_counts_plot <- sp_plot / pt_plot + plot_annotation(tag_levels = "A")

ggsave("Candida_MEC_combined_counts.pdf", both_counts_plot, width = 7, height = 5.2, units = "in")

# Stacked bar plot summing patient series
series_plot <- priority_series %>%
  ggplot(aes(x = genus_species, y = series_count, fill = genus_species)) +
  scale_fill_manual(values = species_colors) +
  scale_x_discrete(limits = species_count$genus_species[sorted_species]) +
  scale_y_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80)) +
  geom_col(colour = "grey26") +
  ylab("Number of serial isolates,\n grouped by case") +
  xlab(NULL) +
  theme_bw() +
  theme(axis.ticks = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(
    angle = 35,
    hjust = 1,
    vjust = 1,
    color = "black",
    size = 11,
    face = "italic"
  )) +
  theme(axis.title = element_text(color = "black", size = 12))
