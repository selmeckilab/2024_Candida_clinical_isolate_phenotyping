## ---------------------------
## Purpose: Plot growth curve summaries and test correlations for all species
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

source("manuscript_redcap_MIC_summary.R")

library(patchwork)
library(ggbeeswarm)
library(correlation)

# Redcap report ID
growth_curves <- "58045"

species_colors <- c(
  "#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
  "#882255", "#BBBBBB", "#AA4499", "#DDCC77", "black"
)

# MIC vals for factor setting
flc_levels <- c("0.5", "1", "2", "4", "8", "16", "32", "64", "128", "256")
mcf_levels <- c("0.016", "0.032", "0.064", "0.125", "0.256", "0.5", "1", ">1", "2")
amb_levels <- c("0.016", "0.032", "0.064", "0.125", "0.256", "0.5", "1")

# function to import report from redcap
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
# Read in GC report and join to MIC data

gc <- import_report(growth_curves) %>%
  filter(
    redcap_repeat_instrument != "NA",
    !primary_id %in% c("AMS5122", "AMS5123")
  ) %>%
  select(
    primary_id, redcap_repeat_instance,
    gc_date, drug_used,
    gc_temp, gc_time,
    k, sem_k,
    r, sem_r,
    t_gen, sem_tgen,
    auc_l, sem_auc_l,
    auc_e, sem_auc_e
  )

gc <- gc %>%
  group_by(primary_id) %>%
  filter(drug_used == "None", gc_date == max(gc_date), gc_time < 24.1)

gc <- gc %>%
  left_join(mic_info, by = join_by(primary_id))

# Quick summary statistics by species
gc_summary <- gc %>%
  filter(drug == "fluconazole") %>%
  select(primary_id, genus_species, k, r, auc_e) %>%
  split(.$genus_species) %>%
  map(summary)

# For ordering species and colors for consistency
species_count <- sample_info %>%
  group_by(genus_species) %>%
  summarize(species_count = n()) %>%
  arrange(desc(species_count))

species_colors <- species_colors %>%
  set_names(species_count$genus_species)

gc$genus_species <- factor(gc$genus_species, levels = species_count$genus_species)

################################################################################
# Growth rate plots

## Add count of samples to box plots at bottom of plot
give.n <- function(x) {
  return(c(y = 0, label = length(x)))
  # experiment with the multiplier to find the perfect position
}

## YPAD growth rate box plot
r_metric <- ggplot(gc, aes(x = genus_species, y = r, color = genus_species)) +
  # geom_beeswarm(size = 0.8, cex = 0.5) +
  geom_boxplot() +
  theme_bw() +
  scale_color_manual(values = species_colors, guide = "none") +
  theme(axis.title = element_text(
    color = "black",
    size = 12
  )) +
  theme(axis.text.x = element_text(
    angle = 30,
    face = "italic",
    hjust = 1,
    vjust = 1,
    color = "black",
    size = 11
  )) +
  xlab(NULL) +
  ylab(expression(paste("Growth rate, ", "hr"^-1)))

r_metric <- r_metric + stat_summary(fun.data = give.n, geom = "text", fun.y = median)

## YPAD growth rate points + error bars
r_summary_plot <- gc %>%
  filter(drug == "fluconazole") %>%
  ggplot(aes(x = genus_species, y = r, color = genus_species)) +
  geom_pointrange(aes(ymin = (r - sem_r), ymax = (r + sem_r)),
    position = position_jitter(width = 0.4),
    alpha = 0.7,
    linetype = "dotted"
  ) +
  theme_bw() +
  scale_color_manual(values = species_colors, guide = "none") +
  theme(axis.text.x = element_text(
    angle = 30,
    face = "italic",
    hjust = 1,
    vjust = 1,
    color = "black",
    size = 10
  )) +
  xlab(NULL) +
  ylab(expression(paste("Growth rate, ", "hr"^-1)))

ggsave("Candida_MEC_YPAD_growth_rate.pdf", r_metric, width = 7, height = 5.5, units = "in")
ggsave("Candida_MEC_growth_rate_summary_with_error_bars.pdf", r_summary_plot, width = 7, height = 5.5, units = "in")

################################################################################
# Test for correlation between growth rate and MIC values

## Subset to relevant species and drug concentrations (species with >1 MIC value)
flc_gc <- gc %>%
  filter(drug == "fluconazole", genus_species != "C. krusei")
flc_gc$genus_species <- factor(flc_gc$genus_species)
flc_gc <- flc_gc %>% mutate(mic50 = replace(mic50, mic50 == ">32", "64")) # Bin C. tropicalis (temp)
flc_gc <- flc_gc %>% mutate(mic50 = replace(mic50, mic50 == "160", "256")) # Bin C. utilis
flc_gc$mic50 <- factor(flc_gc$mic50, levels = flc_levels)
flc_gc <- flc_gc %>%
  mutate(mic_int = as.integer(mic50)) # ordered factor = ordinal variable

flc_mic_counts <- flc_gc %>%
  group_by(genus_species, mic50) %>%
  count() %>%
  group_by(genus_species) %>%
  count()
flc_mic_counts <- flc_mic_counts %>% filter(n > 1) # Get species list with > 1 MIC value

mcf_gc <- gc %>%
  filter(drug == "micafungin")
mcf_gc$mic50 <- factor(mcf_gc$mic50, levels = mcf_levels)
mcf_gc <- mcf_gc %>%
  mutate(mic_int = as.integer(mic50)) # ordered factor = ordinal variable

mcf_mic_counts <- mcf_gc %>%
  group_by(genus_species, mic50) %>%
  count() %>%
  group_by(genus_species) %>%
  count()
mcf_mic_counts <- mcf_mic_counts %>% filter(n > 1) # Get species list with > 1 MIC value

amb_gc <- gc %>%
  filter(drug == "amphotericin B")
amb_gc <- amb_gc %>% mutate(mic90 = replace(mic90, mic90 == "0.023", "0.032"))
amb_gc$mic90 <- factor(amb_gc$mic90, levels = amb_levels)
amb_gc <- amb_gc %>%
  mutate(mic_int = as.integer(mic90)) # ordered factor = ordinal variable

amb_mic_counts <- amb_gc %>%
  group_by(genus_species, mic90) %>%
  count() %>%
  group_by(genus_species) %>%
  count()
amb_mic_counts <- amb_mic_counts %>% filter(n > 1) # Get species list with > 1 MIC value

# Correlation for each drug
flc_gc_vs_mic <- flc_gc %>%
  filter(genus_species %in% flc_mic_counts$genus_species)
flc_gc_vs_mic$genus_species <- factor(flc_gc_vs_mic$genus_species) # Drop species w/o varying MIC
flc_gc_vs_mic_corr <- flc_gc_vs_mic %>%
  group_by(genus_species) %>%
  select(genus_species, primary_id, r, mic_int) %>%
  correlation(method = "spearman")

flc_smg_vs_mic_corr <- flc_gc_vs_mic %>%
  group_by(genus_species) %>%
  select(genus_species, primary_id, mean_smg, mic_int) %>%
  correlation(method = "spearman")

mcf_gc_vs_mic <- mcf_gc %>%
  filter(genus_species %in% mcf_mic_counts$genus_species)
mcf_gc_vs_mic$genus_species <- factor(mcf_gc_vs_mic$genus_species) # Drop species w/o varying MIC
mcf_gc_vs_mic_corr <- mcf_gc_vs_mic %>%
  group_by(genus_species) %>%
  select(genus_species, primary_id, r, mic_int) %>%
  correlation(method = "spearman")

amb_gc_vs_mic <- amb_gc %>%
  filter(genus_species %in% amb_mic_counts$genus_species)
amb_gc_vs_mic$genus_species <- factor(amb_gc_vs_mic$genus_species) # Drop species w/o varying MIC
amb_gc_vs_mic_corr <- amb_gc_vs_mic %>%
  group_by(genus_species) %>%
  select(genus_species, primary_id, r, mic_int) %>%
  correlation(method = "spearman")

corr_results <- list(flc_gc_vs_mic_corr, mcf_gc_vs_mic_corr, amb_gc_vs_mic_corr)
write_xlsx(corr_results, "Candida_MEC_growth_rate_MIC_correlations.xlsx")

################################################################################
# Plot faceted MIC vs growth rate
flc_r_scatter_plot <- flc_gc %>%
  filter(genus_species %in% flc_mic_counts$genus_species) %>%
  ggplot(aes(x = mic50, y = r, color = genus_species)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), ) +
  facet_wrap(~genus_species, nrow = 4) +
  scale_color_manual(values = species_colors, guide = "none") +
  scale_x_discrete(limits = flc_levels, drop = FALSE) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic", size = 9)) +
  theme(axis.text = element_text(size = 7, colour = "black")) +
  xlab(expression(paste("Fluconazole MIC50, ", mu, "g/mL"))) +
  ylab(expression(paste("Mean growth rate, ", "hr"^-1)))

# Binning MCF concentrations for plotting purposes
mcf_gc <- mcf_gc %>%
  mutate(mic50 = replace(mic50, mic50 == "2", ">1"))
mcf_levels <- c("0.016", "0.032", "0.064", "0.125", "0.256", "0.5", "1", ">1")
mcf_gc$mic50 <- factor(mcf_gc$mic50, levels = mcf_levels)

mcf_r_scatter_plot <- mcf_gc %>%
  filter(genus_species %in% mcf_mic_counts$genus_species) %>%
  ggplot(aes(x = mic50, y = r, color = genus_species)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), ) +
  facet_wrap(~genus_species, nrow = 4) +
  scale_color_manual(values = species_colors, guide = "none") +
  scale_x_discrete(
    limits = mcf_levels, drop = FALSE,
    labels = c("0.016", "0.03", "0.06", "0.125", "0.25", "0.5", "1", ">1")
  ) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic", size = 9)) +
  theme(axis.text = element_text(size = 7, colour = "black")) +
  xlab(expression(paste("Micafungin MIC50, ", mu, "g/mL"))) +
  ylab(expression(paste("Mean growth rate, ", "hr"^-1)))

amb_r_scatter_plot <- amb_gc %>%
  filter(genus_species %in% amb_mic_counts$genus_species) %>%
  ggplot(aes(x = mic90, y = r, color = genus_species)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), ) +
  facet_wrap(~genus_species, nrow = 4) +
  scale_color_manual(values = species_colors, guide = "none") +
  scale_x_discrete(
    limits = amb_levels, drop = FALSE,
    labels = c("0.016", "0.03", "0.06", "0.125", "0.25", "0.5", "1")
  ) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic", size = 9)) +
  theme(axis.text = element_text(size = 7, colour = "black")) +
  xlab(expression(paste("Amphotericin B MIC90, ", mu, "g/mL"))) +
  ylab(expression(paste("Mean growth rate, ", "hr"^-1)))


ggsave("Candida_MEC_FLC_r_scatter.pdf", flc_r_scatter_plot, width = 7, height = 5.6, units = "in")
ggsave("Candida_MEC_MCF_r_scatter.pdf", mcf_r_scatter_plot, width = 7, height = 5.6, units = "in")
ggsave("Candida_MEC_AMB_r_scatter.pdf", amb_r_scatter_plot, width = 7, height = 7.5, units = "in")
