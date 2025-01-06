## ---------------------------
## Purpose: Plot timeline serial and recurrent candidemia cases
## Author: Christopher Zajac
## Email: zajac039@umn.edu
## ---------------------------

source("manuscript_changes/MIC_data_summary.R")

species_colors <- c(
  "#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
  "#882255", "#BBBBBB", "#AA4499", "#DDCC77", "black"
)

species_colors <- species_colors %>%
  set_names(species_count$genus_species)

# Get number of isolates in each series
# Order dataframe by species, with polyfungal cases on bottom

redcap_series <- sample_info %>%
  group_by(patient_code) %>%
  reframe(poly_fungal = length(unique(genus_species)), nseries = n(), across()) %>%
  arrange(poly_fungal, genus_species, .by_group = TRUE)

# Factor dataframe by series correctly
redcap_series$patient_code <- factor(redcap_series$patient_code, levels = unique(redcap_series$patient_code))

# Categorize recurrent series
redcap_series$graph <- sapply(redcap_series$relative_days, function(x) ifelse(x <= 25, "Serial", "Recurrent"))

# Filter for serial cases
timeline_plot <- redcap_series %>%
  subset(nseries > 1 & !is.na(patient_code)) %>%
  ggplot() +
  geom_point(aes(
    x = (relative_days), y = patient_code,
    colour = genus_species
  ), alpha = 0.7, size = 1, position = position_jitter(width = 0.25, height = 0.15)) +
  scale_colour_manual(values = species_colors) +
  labs(x = "Relative Day of Isolate Collection", y = "Patient Code", colour = "Species") +
  theme_minimal() +
  theme(
    axis.text = element_text(colour = "black"),
    legend.text = element_text(
      colour = "black", face = "italic"
    )
  ) +
  scale_y_discrete(limits = rev) +
  facet_grid(~ factor(graph, levels = c("Serial", "Recurrent")), scales = "free_x")
ggsave("timeline.pdf", timeline_plot, width = 7, height = 9, units = "in")
ggsave("timeline.tiff", timeline_plot, width = 7, height = 9, units = "in", dpi = 320, bg = "white")
