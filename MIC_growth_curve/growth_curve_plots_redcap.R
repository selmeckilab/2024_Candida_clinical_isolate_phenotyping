## ---------------------------
## Purpose: get growthcurve summaries and basic GC plots
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
## Notes:
## ---------------------------
options(scipen = 999)
## ---------------------------
## load packages
library(jsonlite)
library(paletteer)

# Get spreadsheet cleaning and gc metric functions
source("~/umn/clinical_isolate_redcap/MIC_growth_curve/gc_functions.R")

# Input vars
plate_reader_file <- "data/growth_curve/2023-10-11_MEC_GC30.xlsx"
plate_reader_sn <- ""
gc_date <- str_extract(plate_reader_file, "\\d+-\\d+-\\d+")
facet_colors <- c(paletteer_d("ggthemes::Tableau_20"), paletteer_d("ggsci::category20c_d3"))

api_url <- "https://redcap.ahc.umn.edu/redcap/api/"

## Read in plate, get sample metadata and check fit of logistic curve
plate_od <- clean_growthcurver(plate_reader_file) %>%
  filter(time <= 24.1) # make sure all plates are same time span

sample_data <- samples(plate_reader_file)

plate_data <- SummarizeGrowthByPlate(plate_od)

colnames(plate_data)[colnames(plate_data) == "sample"] <- "well"

# Tidy, add metadata and summarize values
plate_summary <- plate_data %>%
  inner_join(sample_data, by = "well") %>%
  group_by(primary_id) %>%
  filter(!toupper(primary_id) %in% toupper(c("water", "blank", "NA", "ypad"))) %>%
  filter(!is.na(primary_id)) %>%
  # filter(str_detect(well, "F|G|H")) %>% # to filter for specific rows
  summarize(
    mean_k = round(mean(k), digits = 3),
    sem_k = round(stderror(k), digits = 3),
    mean_n0 = round(mean(n0), digits = 3),
    sem_n0 = round(stderror(n0), digits = 3),
    mean_r = round(mean(r), digits = 3),
    sem_r = round(stderror(r), digits = 3),
    mean_tmid = round(mean(t_mid), digits = 3),
    sem_tmid = round(stderror(t_mid), digits = 3),
    mean_tgen = round(mean(t_gen), digits = 3),
    sem_tgen = round(stderror(t_gen), digits = 3),
    mean_auc_l = round(mean(auc_l), digits = 3),
    sem_auc_l = round(stderror(auc_l), digits = 3),
    mean_auc_e = round(mean(auc_e), digits = 3),
    sem_auc_e = round(stderror(auc_e), digits = 3)
  )

# Basic GC plot
input_data <- growth_curve(plate_reader_file)
full_gc_data <- link_metadata(input_data, sample_data, plate = 1)

blank_od <- full_gc_data %>%
  filter(toupper(primary_id) %in% c("YPAD", "BLANK")) %>%
  summarise(mean_blank = mean(OD))

gc_norm <- full_gc_data %>%
  filter(!is.na(primary_id)) %>%
  filter(!toupper(primary_id) %in% toupper(c("ypad", "water", "NA", "blank"))) %>%
  mutate(OD = OD - blank_od$mean_blank)

# Filter for wells and times when samples added at different points
# gc_norm <- gc_norm %>%
#    filter(str_detect(well, "F|G|H")) %>%
#    filter(time >27.18)

gc_mean <- gc_norm %>%
  filter(toupper(primary_id) != "BLANK") %>%
  group_by(primary_id, time) %>%
  mutate(mean_OD = mean(OD, na.rm = TRUE), sd_OD = sd(OD, na.rm = TRUE), se_OD = sd_OD / sqrt(length(OD))) %>%
  mutate(primary_id = as.factor(primary_id)) %>%
  mutate(primary_id = fct_relevel(primary_id, mixedsort))

# Plot by species
spec <- c("Candida glabrata")
ctrl <- "" # or empty str if not plotting
full_plot <- gc_mean %>%
  filter(species %in% spec | species == ctrl) %>%
  ggplot(aes(x = time, y = mean_OD, color = primary_id)) +
  geom_point() +
  geom_errorbar(aes(ymax = (mean_OD + sd_OD), ymin = (mean_OD - sd_OD)), alpha = 0.2, show.legend = FALSE) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0, 1.5)) +
  scale_x_continuous(
    breaks = c(0, 12, 24, 36, 48),
    labels = c("0", "12", "24", "36", "48")
  ) +
  scale_color_manual(values = facet_colors, name = "Isolate ID") +
  theme_grey(base_family = "sans") +
  xlab("Time (hours)") +
  ylab("Mean OD600") +
  labs(title = bquote(italic(.(spec)) ~ "growth in YPAD at 30C"))

ggsave(paste0("images/2023_growth_curves/", gc_date, "_C", str_split_i(spec, " ", 2), "_GC30.png"), full_plot, device = png, width = 8, height = 5.7, units = "in", dpi = 300)

# Facet plot of all isolates on plate
gc_mean %>%
  ggplot(aes(x = time, y = mean_OD, color = species)) +
  geom_point(show.legend = TRUE) +
  geom_errorbar(aes(ymax = (mean_OD + sd_OD), ymin = (mean_OD - sd_OD)), alpha = 0.2, show.legend = FALSE) +
  facet_wrap(~ as.factor(primary_id), nrow = 2) +
  scale_y_continuous(breaks = c(0, 0.5, 1.0, 1.5)) +
  scale_x_continuous(
    breaks = c(0, 12, 24, 36, 48),
    labels = c("0", "12", "24", "36", "48")
  ) +
  scale_color_manual(values = c("#0077bb", "#ee7733", "#8c8c8c", "#33bbee", "#cc3311", "#59A14FFF", "#CDB38B")) +
  xlab("Time (hours)") +
  ylab("Mean OD600") +
  labs(colour = "Species", title = "Growth in YPAD at 30C")

ggsave(paste0("images/2023_growth_curves/", gc_date, "_faceted_MEC_GC30.png"), facet_plot, width = 10, height = 4, units = "in", dpi = 300)

# Redcap upload
redcap_media <- case_when(
  toupper(sample_data$media[1]) == "YPAD" ~ 1,
  toupper(sample_data$media[1]) == "RPMI" ~ 3
)

redcap_temp <- case_when(
  sample_data$temp[1] == 30 ~ 0,
  sample_data$temp[1] == 37 ~ 1,
  sample_data$temp[1] == 35 ~ 2
)

redcap_reader <- case_when(
  plate_reader_sn == "1803065" ~ 0,
  plate_reader_sn == "23012423" ~ 1,
  plate_reader_sn == "23012512" ~ 2,
  plate_reader_sn == "151117B" ~ 3
)

for (i in 1:length(plate_summary$primary_id)) {
  record <- c(
    primary_id = plate_summary$primary_id[i],
    redcap_repeat_instrument = "growth_curve_data",
    redcap_repeat_instance = "new",
    gc_date = gc_date,
    gc_time = round(max(plate_od$time), digits = 2),
    gc_temp = redcap_temp,
    gc_media = redcap_media,
    drug_used = sample_data$drug[1],
    k = plate_summary$mean_k[i],
    n0 = plate_summary$mean_n0[i],
    r = plate_summary$mean_r[i],
    t_mid = plate_summary$mean_tmid[i],
    t_gen = plate_summary$mean_tgen[i],
    auc_l = plate_summary$mean_aucl[i],
    auc_e = plate_summary$mean_auce[i],
    gc_reader = redcap_reader,
    gc_file_name = basename(plate_reader_file),
    growth_curve_data_complete = 2
  )

  result_data <- toJSON(list(as.list(record)), auto_unbox = TRUE)

  formData <- list(
    "token" = Sys.getenv("redcap_api_key"),
    content = "record",
    action = "import",
    format = "json",
    type = "flat",
    overwriteBehavior = "normal",
    forceAutoNumber = "false",
    data = result_data,
    returnContent = "count",
    returnFormat = "json"
  )
  response <- httr::POST(api_url, body = formData, encode = "form")
  result <- httr::content(response)
  print(result)
}
