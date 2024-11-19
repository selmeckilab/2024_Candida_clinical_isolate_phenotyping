## ---------------------------
## Purpose: Summarize and plot MEC isolate (EUCAST) MIC and SMG data by species 
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999) 

# Load packages
library(patchwork)
library(ggbeeswarm)
library(data.table)

source("MIC_data_summary.R")

species_colors <- c("#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
                    "#882255","#BBBBBB", "#AA4499", "#DDCC77", "black")

species_colors <- species_colors %>% 
    set_names(species_count$genus_species)

drugs <- as_labeller(c(fluconazole="Fluconazole", 
                       micafungin="Micafungin",
                       `amphotericin B` = "Amphotericin B"))

full_mic <- data.table::copy(mic_info)

# Bin Etest reults with broth MIC levels
full_mic <- full_mic %>% 
  mutate(mic50 = replace(mic50, mic50==">32", "64")) %>% 
  mutate(mic50 = replace(mic50, mic50=="160", "256"))

full_mic <- full_mic %>% 
  mutate(mic90=replace(mic90, mic90=="0.023", "0.032"))

full_mic <- full_mic %>% 
  mutate(mic50 = case_when((drug=="micafungin" & genus_species=="C. parapsilosis" & mic50=="2") ~ ">1", 
                         .default = mic50))

# Drug-specific levels for plotting
flc_levels <- c("0.5", "1", "2", "4", "8", "16", "32", "64", "128", "256")
mcf_amb_levels <- c("0.016", "0.032", "0.064", "0.125", "0.256", "0.5", "1", ">1")

# Drug-specific DFs
mic_flc <- full_mic %>% 
  filter(drug=="fluconazole")
mic_flc$mic50 <- factor(mic_flc$mic50, levels = flc_levels)

mic_mcf <- full_mic %>% 
  filter(drug=="micafungin")
mic_mcf$mic50 <- factor(mic_mcf$mic50, levels=mcf_amb_levels)

mic_amb <- full_mic %>% 
  filter(drug=="amphotericin B")
mic_amb$mic90 <- factor(mic_amb$mic90, levels = mcf_amb_levels)

# Find the series with > 2-fold changes
flc_step_diffs <- mic_flc %>% 
  filter(!is.na(series_id)) %>% 
  group_by(series_id) %>%
  summarize(diff = max(as.numeric(mic50), na.rm = T) - min(as.numeric(mic50), na.rm = T)) %>%
  filter(diff > 1)

mcf_step_diffs <- mic_mcf %>% 
  filter(!is.na(series_id)) %>% 
  group_by(series_id) %>%
  summarize(diff = max(as.numeric(mic50), na.rm = T) - min(as.numeric(mic50), na.rm = T)) %>%
  filter(diff > 1)

amb_step_diffs <- mic_amb %>% 
  filter(!is.na(series_id)) %>% 
  group_by(series_id) %>%
  summarize(diff = max(as.numeric(mic90), na.rm = T) - min(as.numeric(mic90), na.rm = T)) %>%
  filter(diff > 1)

# Plot per drug, facet by patient
flc_diff_plot <- mic_flc %>%  
  filter(series_id %in% flc_step_diffs$series_id) %>%
  mutate(relative_days = relative_days + 1) %>% 
  ggplot(aes(x = relative_days, y = mic50, color = genus_species)) +
  facet_wrap(.~patient_code) +
  geom_beeswarm(cex = 2) +
  scale_color_manual(values=species_colors, guide = "none") +
  scale_x_continuous(limits = c(1,12), breaks = c(1,3,5,7,9,11)) +
  theme_bw() +
  theme(strip.text = element_text(color = "black")) + 
  theme(axis.text = element_text(size = 8, color = "black")) +
  theme(axis.title = element_text(size=11, color = "black")) +
  ylab("Fluconazole MIC50") +
  xlab(NULL)

mcf_diff_plot <- mic_mcf %>% 
  filter(series_id %in% mcf_step_diffs$series_id) %>%
  mutate(relative_days = relative_days + 1) %>% 
  ggplot(aes(x = relative_days, y = mic50, color = genus_species)) +
  facet_wrap(.~patient_code) +
  geom_beeswarm(cex = 2, show.legend = TRUE) +
  scale_color_manual(values=species_colors, name = NULL, drop = FALSE) +
  scale_x_continuous(limits = c(1,12), breaks = c(1,3,5,7,9,11)) +
  theme_bw() +
  theme(strip.text = element_text(color = "black")) +
  theme(legend.text = element_text(color = "black", face = "italic")) +
  theme(axis.text = element_text(size = 8, color = "black")) +
  theme(axis.title = element_text(size=11, color = "black")) +
  ylab("Micafungin MIC50") +
  xlab(NULL)

amb_diff_plot <- mic_amb %>% 
  filter(series_id %in% amb_step_diffs$series_id) %>%
  mutate(relative_days = relative_days + 1) %>% 
  ggplot(aes(x = relative_days, y = mic90, color = genus_species)) +
  facet_wrap(.~patient_code) +
  geom_beeswarm(cex = 2) +
  scale_color_manual(values=species_colors, guide = "none") +
  scale_x_continuous(limits = c(1,12), breaks = c(1,3,5,7,9,11)) +
  theme_bw() +
  theme(strip.text = element_text(color = "black")) + 
  theme(axis.text = element_text(size = 8, color = "black")) +
  theme(axis.title = element_text(size=11, color = "black")) +
  ylab("Amphotericin B MIC90") +
  xlab("Collection day")

combined_diffs <- (flc_diff_plot/(mcf_diff_plot + plot_layout(ncol = 3, nrow = 1))/amb_diff_plot) +
  plot_layout(heights = c(2,1,2), guides = "collect") + plot_annotation(tag_levels = "A")

ggsave("MEC_within_patient_MIC_differences.pdf", combined_diffs, width = 7, height = 9, units = "in")