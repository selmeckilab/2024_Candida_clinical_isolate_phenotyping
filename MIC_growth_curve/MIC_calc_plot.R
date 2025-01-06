## ---------------------------
## Purpose: Calculate, plot and save MIC and SMG data
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

library(readxl)
library(writexl)
library(tidyverse)
library(gtools)
library(patchwork)

# Enter drug and spreadsheet data
drug_used <- "flc"
mic_breakpoint <- 0.5
replicates <- 3

save_dir <- ""

mic_spreadsheet <- "~/umn/data/MIC/2024-07-18_EW_Cglabrata_MIC24_RPMI35.xlsx"
smg_spreadsheet <- "~/umn/data/MIC/2024-07-19_EW_Cglabrata_SMG48_RPMI35.xlsx"

od_tab <- 1 # Excel tab of OD values
metadata_tab <- 2 # Excel tab of metadata


# For redcap uploading, and file naming if desired
mic_date <- str_extract(mic_spreadsheet, "\\d+-\\d+-\\d+")

# Std error function
stderror <- function(x) sd(x) / sqrt(length(x))

################################################################################
# MIC: read metadata, convert everything to strings after to avoid weird floats
mic_plates <- read_excel(mic_spreadsheet,
  sheet = metadata_tab,
  col_names = TRUE,
  col_types = "list"
)

mic_plates <- mic_plates %>%
  mutate(across(1:ncol(mic_plates), as.character))


# Grab these values for any redcap uploading
media <- mic_plates$media[1]
temp <- mic_plates$temp[1]

# Read in OD vals for all drugs
for (a in 1:replicates) {
  mic_plates <-
    mic_plates %>%
    mutate(read_excel(mic_spreadsheet,
      sheet = od_tab,
      range = mic_plates[[drug_used]][a]
    ) %>%
      pivot_longer(everything(),
        names_to = NULL,
        values_to = "OD600"
      ))
  mic_plates <- mic_plates %>%
    rename("{drug_used}_{a}_OD600" := "OD600") # this rename gives replicate info
}


# Select columns, using regex to keep only drug of interest
mic_plates <- mic_plates %>%
  select(c(well, strain, grep(paste0(drug_used, "_"), colnames(mic_plates), value = TRUE)))

# Simplify column name
names(mic_plates)[names(mic_plates) == paste0(drug_used, "_concentration")] <- "concentration"

# Set drug concentration as factor, and specify the order
mic_plates$concentration <- factor(mic_plates$concentration, levels = mixedsort(unique(mic_plates$concentration)))
mic_plates <- mic_plates %>%
  pivot_longer(-c(well, strain, concentration), names_to = "replicate", values_to = "OD")

# Get average background, remove unusually high values
mic_blank <- mic_plates %>%
  filter(concentration == "blank") %>%
  filter(OD < 0.2) %>%
  summarise(blank_OD = mean(OD))

# Background correction
mic_plates <- mic_plates %>%
  filter(concentration != "blank") %>%
  group_by(replicate, strain, concentration) %>%
  mutate(corrected_OD = OD - mic_blank$blank_OD)

# Normalize by control well
mic_plates <- mic_plates %>%
  filter(strain != "NA") %>%
  group_by(replicate, strain) %>%
  mutate(relative_growth = round(corrected_OD / corrected_OD[concentration == "0" | concentration == "0.0"], digits = 3)) %>%
  droplevels() # this drops the blank level of the concentration factor

# Check for and correct negative relative growth
mic_plates <- mic_plates %>%
  mutate(relative_growth = ifelse(relative_growth < 0, 0, relative_growth))

# Add factor for MIC value above plate max
mic_plates$concentration <- fct_expand(mic_plates$concentration, paste0(">", levels(mic_plates$concentration)[length(levels(mic_plates$concentration))]))

# Summarize - average the replicates
summary_vals <- mic_plates %>%
  group_by(strain, concentration) %>%
  summarise(
    mean_corrected_OD = round(mean(corrected_OD), digits = 3),
    sd_corrected_OD = round(sd(corrected_OD), digits = 3),
    mean_relative_growth = round(mean(relative_growth), digits = 3),
    sd_relative_growth = round(sd(relative_growth), digits = 3)
  )

# Filter for MIC values per strain
mic_vals <- summary_vals %>%
  select(strain, concentration, mean_relative_growth, sd_relative_growth) %>%
  filter(mean_relative_growth < mic_breakpoint) %>%
  select(strain, concentration) %>%
  group_by(strain) %>%
  slice_head() %>%
  right_join(distinct(summary_vals, pick("strain")), by = "strain") # if any strains are > than max concentration, this adds the strain name back

# In case of strains above max concentration, assign the new MIC value that you created on line 99
mic_vals <- mic_vals %>%
  replace_na(list(concentration = levels(mic_vals$concentration)[length(levels(mic_vals$concentration))])) %>%
  arrange(strain) %>%
  left_join(summary_vals)

# Get wells for SMG
smg_wells <- summary_vals %>%
  group_by(strain) %>%
  filter(mean_relative_growth < mic_breakpoint) %>%
  select(strain, concentration)

################################################################################
# MIC heatmap
mic_plot <- summary_vals %>%
  group_by(strain) %>% # Can add a step to filter out strains or set a different order here
  ggplot(aes(x = concentration, y = strain, fill = mean_relative_growth)) +
  geom_raster(hjust = 1.0) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x = element_text(family = "Helvetica", color = "black", size = 10),
    axis.text.y = element_text(family = "Helvetica", color = "black", size = 12, vjust = 0.4),
    legend.text = element_text(family = "Helvetica", color = "black"),
    plot.title = element_text(family = "Helvetica", color = "black", size = 13)
  ) +
  guides(fill = guide_legend(title = "Relative \ngrowth")) +
  ylab(NULL) +
  xlab(paste(toupper(drug_used), "concentration")) +
  theme(
    plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm"),
    axis.title.y = element_text(hjust = 0.5, margin = margin(0, 10, 0, 0))
  ) +
  annotate("segment",
    x = mic_vals$concentration,
    xend = mic_vals$concentration,
    y = (seq(0.5, n_distinct(mic_vals$strain), by = 1)),
    yend = (seq(1.5, n_distinct(mic_vals$strain) + 1, by = 1)),
    color = "yellow",
    linewidth = 1.5
  )

################################################################################
# SMG values repeat as above: read metadata, convert everything to strings
smg_plates <- read_excel(smg_spreadsheet,
  sheet = metadata_tab,
  col_names = TRUE,
  col_types = "list"
)

smg_plates <- smg_plates %>%
  mutate(across(1:ncol(smg_plates), as.character))

# Read in OD vals for all drugs
for (a in 1:replicates) {
  smg_plates <-
    smg_plates %>%
    mutate(read_excel(smg_spreadsheet,
      sheet = od_tab,
      range = smg_plates[[drug_used]][a]
    ) %>%
      pivot_longer(everything(),
        names_to = NULL,
        values_to = "OD600"
      ))
  smg_plates <- smg_plates %>%
    rename("{drug_used}_{a}_OD600" := "OD600") # rename to include plate number
}


# Select columns w/regex to keep drug of interest, simplify column name, set drug levels as factors
smg_plates <- smg_plates %>%
  select(c(well, strain, grep(paste0(drug_used, "_"), colnames(smg_plates), value = TRUE)))

names(smg_plates)[names(smg_plates) == paste0(drug_used, "_concentration")] <- "concentration"

smg_plates$concentration <- factor(smg_plates$concentration, levels = mixedsort(unique(smg_plates$concentration)))

smg_plates <- smg_plates %>%
  pivot_longer(-c(well, strain, concentration), names_to = "replicate", values_to = "OD")

# Get average background, remove unusually high values
smg_blank <- smg_plates %>%
  filter(concentration == "blank") %>%
  filter(OD < 0.2) %>%
  summarise(blank_OD = mean(OD))

# Background correction
smg_plates <- smg_plates %>%
  filter(concentration != "blank") %>%
  group_by(replicate, strain, concentration) %>%
  mutate(corrected_OD = OD - smg_blank$blank_OD)

# Normalize by control well
smg_plates <- smg_plates %>%
  filter(strain != "NA") %>%
  group_by(replicate, strain) %>%
  mutate(relative_growth = round(corrected_OD / corrected_OD[concentration == "0" | concentration == "0.0"], digits = 3)) %>%
  droplevels()

smg_plates <- smg_plates %>%
  mutate(relative_growth = ifelse(relative_growth < 0, 0, relative_growth))

# Add factor for smg value above plate max
smg_plates$concentration <- fct_expand(smg_plates$concentration, paste0(">", levels(smg_plates$concentration)[length(levels(smg_plates$concentration))]))

# Summarize and subset to MIC and above
smg_summary_vals <- smg_plates %>%
  group_by(strain, concentration) %>%
  summarise(
    mean_corrected_OD = round(mean(corrected_OD), digits = 3),
    sd_corrected_OD = round(sd(corrected_OD), digits = 3),
    mean_relative_growth = round(mean(relative_growth), digits = 3),
    sd_relative_growth = round(sd(relative_growth), digits = 3)
  ) %>%
  right_join(smg_wells)

# Final dataframe of mean SMG (calculated from means of each well) and std error
# Any strains above plate max are added back in with NA value
final_smg <- smg_summary_vals %>%
  group_by(strain) %>%
  summarise(
    mean_smg = mean(mean_relative_growth),
    sem_smg = round(stderror(mean_relative_growth), digits = 3)
  ) %>%
  right_join(distinct(summary_vals, pick("strain")), by = "strain") %>%
  arrange(strain)


################################################################################
# SMG barplot
smg_plot <- final_smg %>%
  ggplot(aes(y = strain, x = mean_smg)) +
  geom_col(fill = "black") +
  theme_minimal() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(family = "Helvetica", color = "black", size = 10),
    axis.line.y = element_line(color = "black", linewidth = 0.2),
    axis.line.x = element_line(color = "black", linewidth = 0.3),
    plot.title = element_text(family = "Helvetica", color = "black", size = 13)
  ) +
  xlab("SMG") +
  ylab(NULL) +
  # expand_limits(y = 1.0) +
  coord_fixed(ratio = 0.8) + # ADJUST THE WIDTH HERE: smaller number means wider bar graph
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = c("0", "", "0.5", "", "1")
  )

# Pull together with patchwork
full_plot <- (mic_plot + smg_plot) + plot_layout(guides = "collect")

# Save final plot
ggsave(paste0(save_dir, mic_date, "_MEC_", toupper(drug_used), "_MIC.png"),
  full_plot,
  width = 7.6,
  height = 4.5,
  units = "in",
  device = png,
  bg = "white",
  dpi = 300
)

################################################################################
# Write out data to new spreadsheet with MIC date

# Only need to make a list if you want to save multiple dataframes as tabs
# I renamed these, which says something about my variable choices
outfiles <- list(
  summary_24hr = summary_vals,
  mic = mic_vals,
  summary_48hr = smg_summary_vals,
  smg = final_smg
)

write_xlsx(outfiles, path = sprintf("%s%s_%s_MICs.xlsx", save_dir, mic_date, toupper(drug_used)))
