## ---------------------------
## Purpose: Summarize and plot MEC isolate (EUCAST) MIC and SMG data by species
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

source("MIC_data_summary.R")

# Load packages
library(patchwork)
library(ggbeeswarm)
library(data.table)

species_colors <- c(
  "#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
  "#882255", "#BBBBBB", "#AA4499", "#DDCC77", "black"
)

species_colors <- species_colors %>%
  set_names(species_count$genus_species)

drugs <- as_labeller(c(
  fluconazole = "Fluconazole",
  micafungin = "Micafungin",
  `amphotericin B` = "Amphotericin B"
))

# All values, MICs not binned
full_mic <- data.table::copy(mic_info)

# Bin Etest reults with broth MIC levels and bin extended FLC values with original
# screening concentrations
mic_info <- mic_info %>%
  mutate(mic90 = replace(mic90, mic90 == "0.023", "0.032"))

mic_info <- mic_info %>%
  mutate(mic50 = replace(mic50, mic50 == "0.047", "0.032")) %>%
  mutate(mic50 = replace(mic50, mic50 == "160", ">32")) %>%
  mutate(mic50 = replace(mic50, mic50 == "128", ">32")) %>%
  mutate(mic50 = replace(mic50, mic50 == "64", ">32")) %>%
  mutate(mic50 = replace(mic50, mic50 == "256", ">32"))

# Bin C. parapsilosis MCF values with original screening concentrations
mic_info <- mic_info %>%
  mutate(mic50 = case_when((drug == "micafungin" & genus_species == "C. parapsilosis" & mic50 == "2") ~ ">1",
    .default = mic50
  ))

# Drug-specific levels
flc_levels <- c("0.5", "1", "2", "4", "8", "16", "32", ">32")
mcf_amb_levels <- c("0.016", "0.032", "0.064", "0.125", "0.256", "0.5", "1", ">1")

# Drug-specific DFs, binned
mic_flc <- mic_info %>%
  filter(drug == "fluconazole")
mic_flc$mic50 <- factor(mic_flc$mic50, levels = flc_levels)

mic_mcf <- mic_info %>%
  filter(drug == "micafungin")
mic_mcf$mic50 <- factor(mic_mcf$mic50, levels = mcf_amb_levels)

mic_amb <- mic_info %>%
  filter(drug == "amphotericin B")
mic_amb$mic90 <- factor(mic_amb$mic90, levels = mcf_amb_levels)

################################################################################
# Bubble plots
flc_sum <- mic_flc %>%
  group_by(mic50, genus_species) %>%
  count()
flc_bubble <- ggplot(flc_sum, aes(x = genus_species, y = mic50, size = n, color = genus_species)) +
  geom_point() +
  theme_bw() +
  xlab(NULL) +
  ylab(expression(paste("Fluconazole MIC50, ", mu, "g/mL"))) +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.y = element_text(color = "black", size = 11)) +
  theme(axis.text.x = element_text(
    angle = 35,
    hjust = 1,
    vjust = 1,
    color = "black",
    size = 10,
    face = "italic"
  )) +
  scale_color_manual(values = species_colors, guide = "none") +
  scale_size_continuous(range = c(1, 11), breaks = c(1, 15, 30, 45, 60, 75, 90), limits = c(1, 100)) +
  guides(size = guide_legend(title = "Isolates (N)")) +
  annotate("segment", x = 0.5, xend = 1.3, y = 4.5, yend = 4.5, linetype = 2) +
  annotate("segment", x = 1.6, xend = 2.4, y = 6.5, yend = 6.5, linetype = 2) +
  annotate("segment", x = 2.6, xend = 4.3, y = 4.5, yend = 4.5, linetype = 2) +
  annotate("segment", x = 5.6, xend = 11.4, y = 4.5, yend = 4.5, linetype = 2)

mcf_sum <- mic_mcf %>%
  group_by(mic50, genus_species) %>%
  count()
mcf_bubble <- ggplot(mcf_sum, aes(x = genus_species, y = mic50, size = n, color = genus_species)) +
  geom_point() +
  theme_bw() +
  xlab(NULL) +
  ylab(expression(paste("Micafungin MIC50, ", mu, "g/mL"))) +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.y = element_text(color = "black", size = 11)) +
  theme(axis.text.x = element_text(
    angle = 35,
    hjust = 1,
    vjust = 1,
    color = "black",
    size = 10,
    face = "italic"
  )) +
  scale_color_manual(values = species_colors, guide = "none") +
  scale_size_continuous(range = c(1, 11), breaks = c(1, 15, 30, 45, 60, 75, 90), limits = c(1, 100)) +
  guides(size = guide_legend(title = "Isolates (N)")) +
  annotate("segment", x = 0.5, xend = 1.3, y = 1.5, yend = 1.5, linetype = 2) +
  annotate("segment", x = 1.6, xend = 2.4, y = 2.5, yend = 2.5, linetype = 2) +
  annotate("segment", x = 2.6, xend = 3.4, y = 8.5, yend = 8.5, linetype = 2)

amb_sum <- mic_amb %>%
  group_by(mic90, genus_species) %>%
  count()
amb_bubble <- ggplot(amb_sum, aes(x = genus_species, y = mic90, size = n, color = genus_species)) +
  geom_point() +
  theme_bw() +
  xlab(NULL) +
  ylab(expression(paste("Amphotericin B MIC90, ", mu, "g/mL"))) +
  theme(panel.grid.major = element_blank()) +
  theme(axis.text.y = element_text(color = "black", size = 11)) +
  theme(axis.text.x = element_text(
    angle = 35,
    hjust = 1,
    vjust = 1,
    color = "black",
    size = 10,
    face = "italic"
  )) +
  scale_color_manual(values = species_colors, guide = "none") +
  scale_size_continuous(range = c(1, 11), breaks = c(1, 6, 12, 18, 24, 30, 36), limits = c(1, 100)) +
  guides(size = guide_legend(title = "Isolates (N)")) +
  annotate("segment", x = 0.5, xend = 3.4, y = 7.5, yend = 7.5, linetype = 2) +
  annotate("segment", x = 4.6, xend = 5.4, y = 7.5, yend = 7.5, linetype = 2) +
  annotate("segment", x = 7.6, xend = 8.4, y = 7.5, yend = 7.5, linetype = 2) +
  annotate("segment", x = 9.6, xend = 10.4, y = 7.5, yend = 7.5, linetype = 2)

ggsave("MEC_FLC_summary.pdf", flc_bubble, width = 7, height = 5, units = "in", dpi = 320)
ggsave("MEC_MCF_summary.pdf", mcf_bubble, width = 7, height = 5, units = "in", dpi = 320)
ggsave("MEC_AMB_summary.pdf", amb_bubble, width = 7, height = 5, units = "in", dpi = 320)

################################################################################
# MIC barplot per drug, faceted
# Two y-axis scales combined with patchwork
flc_fig1 <- ggplot(mic_flc %>% filter(genus_species %in% c(
  "C. albicans",
  "C. glabrata",
  "C. parapsilosis",
  "C. lusitaniae"
)), aes(x = mic50)) +
  geom_bar(aes(fill = genus_species)) +
  facet_wrap(. ~ genus_species, nrow = 1) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic", color = "black")) +
  theme(axis.text.y = element_text(size = 8, color = "black")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  theme(panel.grid.major.x = element_blank()) +
  ylab("Isolates (N)") +
  xlab(NULL) +
  geom_vline(
    data = filter(mic_flc, genus_species == "C. glabrata"),
    aes(xintercept = 6.5),
    linetype = 2
  ) +
  geom_vline(
    data = filter(mic_flc, genus_species %in% c(
      "C. albicans",
      "C. lusitaniae",
      "C. parapsilosis"
    )),
    aes(xintercept = 4.5), linetype = 2
  )

flc_fig2 <- ggplot(mic_flc %>% filter(genus_species %in% c(
  "C. krusei",
  "C. orthopsilosis",
  "C. kefyr",
  "C. tropicalis",
  "C. utilis",
  "C. dubliniensis",
  "C. nivariensis"
)), aes(x = mic50)) +
  geom_bar(aes(fill = genus_species)) +
  facet_wrap(. ~ genus_species, nrow = 2) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic", color = "black")) +
  theme(axis.text.y = element_text(size = 8, color = "black")) +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 90)) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  theme(panel.grid.major.x = element_blank()) +
  ylab(NULL) +
  xlab(expression(paste("Fluconazole MIC50, ", mu, "g/mL"))) +
  geom_vline(
    data = filter(mic_flc, genus_species %in% c(
      "C. orthopsilosis",
      "C. kefyr",
      "C. tropicalis",
      "C. utilis",
      "C. dubliniensis",
      "C. nivariensis"
    )),
    aes(xintercept = 4.5), linetype = 2
  )


mcf_fig1 <- ggplot(mic_mcf %>% filter(genus_species %in% c(
  "C. albicans",
  "C. glabrata",
  "C. parapsilosis",
  "C. lusitaniae"
)), aes(x = mic50)) +
  geom_bar(aes(fill = genus_species)) +
  facet_wrap(. ~ genus_species, nrow = 1) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic", color = "black")) +
  theme(axis.text.y = element_text(size = 8, color = "black")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  theme(panel.grid.major.x = element_blank()) +
  ylab("Isolates (N)") +
  xlab(NULL) +
  geom_vline(
    data = filter(mic_mcf, genus_species == "C. albicans"),
    aes(xintercept = 1.5), linetype = 2
  ) +
  geom_vline(
    data = filter(mic_mcf, genus_species == "C. glabrata"),
    aes(xintercept = 2.5), linetype = 2
  ) +
  geom_vline(
    data = filter(mic_mcf, genus_species == "C. parapsilosis"),
    aes(xintercept = 8.5), linetype = 2
  )

mcf_fig2 <- ggplot(mic_mcf %>% filter(genus_species %in% c(
  "C. krusei",
  "C. orthopsilosis",
  "C. kefyr",
  "C. tropicalis",
  "C. utilis",
  "C. dubliniensis",
  "C. nivariensis"
)), aes(x = mic50)) +
  geom_bar(aes(fill = genus_species)) +
  facet_wrap(. ~ genus_species, nrow = 2) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_x_discrete(drop = FALSE, labels = c("<0.03", "0.03", "0.06", "0.12", "0.25", "0.5", "1", ">1")) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic", color = "black")) +
  theme(axis.text.y = element_text(size = 8, color = "black")) +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 90)) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  theme(panel.grid.major.x = element_blank()) +
  ylab(NULL) +
  xlab(expression(paste("Micafungin MIC50, ", mu, "g/mL")))

amb_fig1 <- ggplot(mic_amb %>% filter(genus_species %in% c(
  "C. albicans",
  "C. glabrata",
  "C. parapsilosis",
  "C. lusitaniae"
)), aes(x = mic90)) +
  geom_bar(aes(fill = genus_species)) +
  facet_wrap(. ~ genus_species, nrow = 1) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic", color = "black")) +
  theme(axis.text.y = element_text(size = 8, color = "black")) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  theme(panel.grid.major.x = element_blank()) +
  ylab("Isolates (N)") +
  xlab(NULL) +
  geom_vline(
    data = filter(mic_amb, genus_species %in% c(
      "C. albicans",
      "C. glabrata",
      "C. parapsilosis"
    )),
    aes(xintercept = 7.5), linetype = 2
  )

amb_fig2 <- ggplot(mic_amb %>% filter(genus_species %in% c(
  "C. krusei",
  "C. orthopsilosis",
  "C. kefyr",
  "C. tropicalis",
  "C. utilis",
  "C. dubliniensis",
  "C. nivariensis"
)), aes(x = mic90)) +
  geom_bar(aes(fill = genus_species)) +
  facet_wrap(. ~ genus_species, nrow = 2) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_x_discrete(drop = FALSE, labels = c("<0.03", "0.03", "0.06", "0.12", "0.25", "0.5", "1", ">1")) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic", color = "black")) +
  theme(axis.text.y = element_text(size = 8, color = "black")) +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 90)) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  theme(panel.grid.major.x = element_blank()) +
  ylab(NULL) +
  xlab(expression(paste("Amphotericin B MIC90, ", mu, "g/mL"))) +
  geom_vline(
    data = filter(mic_amb, genus_species %in% c(
      "C. krusei",
      "C. tropicalis",
      "C. dubliniensis"
    )),
    aes(xintercept = 7.5), linetype = 2
  )

flc_facet <- (flc_fig1 + theme(plot.margin = margin(b = 0))) / (flc_fig2) + plot_layout(heights = c(1, 2), axis_titles = "collect")

mcf_facet <- (mcf_fig1 + theme(plot.margin = margin(b = 0))) / (mcf_fig2) + plot_layout(heights = c(1, 2), axis_titles = "collect")

amb_facet <- (amb_fig1 + theme(plot.margin = margin(b = 0))) / (amb_fig2) + plot_layout(heights = c(1, 2), axis_titles = "collect")

ggsave("flc_facet_two_y_axis.pdf", flc_facet, width = 7, height = 7.5, units = "in")
ggsave("mcf_facet_two_y_axis.pdf", mcf_facet, width = 7, height = 7.5, units = "in")
ggsave("amb_facet_two_y_axis.pdf", amb_facet, width = 7, height = 7.5, units = "in")

################################################################################
# MIC barplots, facet by spp and combine all drugs with patchwork
flc1 <- ggplot(mic_flc %>% filter(genus_species %in% c(
  "C. albicans",
  "C. glabrata",
  "C. parapsilosis",
  "C. lusitaniae"
)), aes(x = mic50)) +
  geom_bar(aes(fill = genus_species), just = 1) +
  scale_fill_manual(values = species_colors, guide = "none") +
  facet_wrap(. ~ genus_species, nrow = 1) +
  theme_bw() +
  xlab(NULL) +
  ylab("Fluconazole\nIsolates (N)") +
  theme(strip.text = element_text(face = "italic", size = 9)) +
  theme(panel.grid.major.x = element_blank()) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = -0.1, size = 8, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  geom_vline(
    data = filter(mic_flc, genus_species %in% c(
      "C. albicans",
      "C. lusitaniae",
      "C. parapsilosis"
    )),
    aes(xintercept = "4"), linetype = 2
  ) +
  geom_vline(
    data = filter(mic_info, genus_species == "C. glabrata"),
    aes(xintercept = "16"), linetype = 2
  )

flc2 <- mic_flc %>%
  filter(genus_species %in% c(
    "C. krusei",
    "C. orthopsilosis",
    "C. kefyr",
    "C. tropicalis",
    "C. utilis",
    "C. dubliniensis",
    "C. nivariensis"
  )) %>%
  ggplot(aes(x = mic50)) +
  geom_bar(aes(fill = genus_species), just = 1) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_x_discrete(limits = flc_levels) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
  facet_wrap(. ~ genus_species, nrow = 1) +
  theme_bw() +
  xlab(NULL) +
  ylab(NULL) +
  # ylab("\n\nFluconazole") +
  theme(strip.text = element_text(face = "italic", size = 9)) +
  theme(panel.grid.major.x = element_blank()) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = -0.1, size = 8, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  geom_vline(
    data = filter(mic_flc, genus_species %in% c(
      "C. tropicalis",
      "C. kefyr",
      "C. orthopsilosis",
      "C. dubliniensis",
      "C. nivariensis",
      "C. utilis"
    )),
    aes(xintercept = "4"), linetype = 2
  )

mcf1 <- ggplot(mic_mcf %>% filter(genus_species %in% c(
  "C. albicans",
  "C. glabrata",
  "C. parapsilosis",
  "C. lusitaniae"
)), aes(x = mic50)) +
  geom_bar(aes(fill = genus_species), just = 1) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_x_discrete(limits = mcf_amb_levels) +
  facet_wrap(. ~ genus_species, nrow = 1, scales = "free_x") +
  theme_bw() +
  xlab(NULL) +
  ylab("Micafungin\nIsolates (N)") +
  theme(strip.text = element_blank()) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = -0.1, size = 8, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  theme(panel.grid.major.x = element_blank()) +
  geom_vline(
    data = filter(mic_mcf, genus_species == "C. albicans"),
    aes(xintercept = "0.016"), linetype = 2
  ) +
  geom_vline(
    data = filter(mic_mcf, genus_species == "C. glabrata"),
    aes(xintercept = "0.032"), linetype = 2
  ) +
  geom_vline(
    data = filter(mic_mcf, genus_species == "C. parapsilosis"),
    aes(xintercept = ">1"), linetype = 2
  )

mcf2 <- ggplot(mic_mcf %>% filter(genus_species %in% c(
  "C. krusei",
  "C. orthopsilosis",
  "C. kefyr",
  "C. tropicalis",
  "C. utilis",
  "C. dubliniensis",
  "C. nivariensis"
)), aes(x = mic50)) +
  geom_bar(aes(fill = genus_species), just = 1) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_x_discrete(limits = mcf_amb_levels) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
  facet_wrap(. ~ genus_species, nrow = 1) +
  theme_bw() +
  xlab(NULL) +
  ylab(NULL) +
  # ylab("Number of isolates\n\nMicafungin")+
  theme(strip.text = element_blank()) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = -0.1, size = 8, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  theme(panel.grid.major.x = element_blank())

amb1 <- ggplot(mic_amb %>% filter(genus_species %in% c(
  "C. albicans",
  "C. glabrata",
  "C. parapsilosis",
  "C. lusitaniae"
)), aes(x = mic90)) +
  geom_bar(aes(fill = genus_species), just = 1) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_x_discrete(limits = mcf_amb_levels) +
  facet_wrap(. ~ genus_species, nrow = 1) +
  theme_bw() +
  # xlab(expression(paste("MIC, ", mu,"g/mL"))) +
  xlab(NULL) +
  ylab("Amphotericin B\nIsolates (N)") +
  theme(strip.text = element_blank()) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = -0.1, size = 8, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  theme(panel.grid.major.x = element_blank()) +
  geom_vline(
    data = filter(mic_amb, genus_species %in% c(
      "C. albicans",
      "C. glabrata",
      "C. parapsilosis"
    )),
    aes(xintercept = "1"), linetype = 2
  )

amb2 <- ggplot(mic_amb %>% filter(genus_species %in% c(
  "C. krusei",
  "C. orthopsilosis",
  "C. kefyr",
  "C. tropicalis",
  "C. utilis",
  "C. dubliniensis",
  "C. nivariensis"
)), aes(x = mic90)) +
  geom_bar(aes(fill = genus_species), just = 1) +
  scale_fill_manual(values = species_colors, guide = "none") +
  scale_x_discrete(limits = mcf_amb_levels) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
  facet_wrap(. ~ genus_species, nrow = 1) +
  theme_bw() +
  xlab(expression(paste("MIC, ", mu, "g/mL"))) +
  ylab(NULL) +
  # ylab("\n\nAmphotericin B") +
  theme(strip.text = element_blank()) +
  theme(axis.title = element_text(size = 11, color = "black")) +
  theme(axis.title.x = element_text(vjust = -1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = -0.1, size = 8, color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
  theme(panel.grid.major.x = element_blank()) +
  geom_vline(
    data = filter(mic_amb, genus_species %in% c(
      "C. krusei",
      "C. tropicalis",
      "C. dubliniensis"
    )),
    aes(xintercept = "1"), linetype = 2
  )


full_plot <- (flc1 + flc2 + plot_layout(widths = c(1, 2))) / (mcf1 + mcf2 + plot_layout(widths = c(1, 2))) / (amb1 + amb2 + plot_layout(widths = c(1, 2)))

ggsave("MEC_all_drugs_all_spp.pdf", full_plot)

################################################################################
# SMG plotting

## Add count of samples to box plots at bottom of plot
give.n <- function(x) {
  return(c(y = 0, label = length(x)))
  # experiment with the multiplier to find the perfect position
}

## Box plot version
smg_boxplot <- mic_flc %>%
  filter(genus_species != "C. krusei") %>%
  filter(!primary_id %in% flc_carrying_cap$primary_id) %>%
  ggplot(aes(x = genus_species, y = mean_smg, color = genus_species)) +
  geom_boxplot() +
  # geom_violin() +
  # geom_beeswarm(size=1, cex=0.8) +
  # geom_point(data = filter(mic_info, genus_species=="C. nivariensis")) +
  scale_y_continuous(limits = c(0, 1)) +
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
  ylab("Fluconazole SMG")

smg_plot <- smg_plot + stat_summary(fun.data = give.n, geom = "text", fun.y = median)

## Point + error bar version
smg_pointrange_plot <- mic_flc %>%
  filter(genus_species != "C. krusei") %>%
  filter(!primary_id %in% flc_carrying_cap$primary_id) %>%
  ggplot(aes(x = genus_species, y = mean_smg, color = genus_species)) +
  geom_pointrange(aes(ymin = (mean_smg - sem_smg), ymax = (mean_smg + sem_smg)),
    position = position_jitter(width = 0.4),
    alpha = 0.7,
    linetype = "dotted"
  ) +
  theme_bw() +
  # scale_y_continuous(limits = c(0, 1.1)) +
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
  ylab("Fluconazole SMG")


ggsave("Candida_MEC_FLC_SMG_boxplot.pdf", smg_boxplot, width = 7, height = 5.5, units = "in")
ggsave("Candida_MEC_FLC_SMG_pointrange.pdf", smg_pointrange_plot, width = 7, height = 5.5, units = "in")

################################################################################
# SMG vs MIC scatterplot
smg_mic_plot <- mic_info %>%
  filter(!primary_id %in% flc_carrying_cap$primary_id, drug == "fluconazole", genus_species == "C. glabrata") %>%
  ggplot(aes(y = mean_smg, x = mic50, drop = FALSE)) +
  # geom_violin(aes(fill = genus_species)) +
  geom_point(position = position_jitter(seed = 1, width = 0.2), aes(color = genus_species)) +
  # facet_wrap(.~genus_species) +
  scale_color_manual(values = species_colors, guide = "none") +
  scale_x_discrete(limits = flc_levels, drop = FALSE) +
  theme_bw() +
  theme(strip.text = element_text(face = "italic", color = "black")) +
  theme(axis.text = element_text(size = 7.5, color = "black")) +
  xlab(expression(paste("Fluconazole MIC50, ", mu, "g/mL"))) +
  ylab("Mean SMG")

ggsave("Cglabrata_SMG_MIC_scatterplot.pdf", smg_mic_plot, width = 6, height = 5, units = "in")
