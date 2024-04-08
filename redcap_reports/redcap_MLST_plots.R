## ---------------------------
## Purpose: MLST bar plots
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
## Load packages
library(tidyverse)
library(paletteer)

# Redcap report IDs
samples <- '58043'
calb_mlst <- '58053'
cglab_mlst <- '58052'

patient_colors <- c(paletteer_d("ggsci::default_igv"),
                    "#838b8b",
                    "black",
                    "#cdcdb4",
                    "#155F83")

# Function to import report from redcap
import_report <- function(report_number) {
  url <- "https://redcap.ahc.umn.edu/redcap/api/"
  formData <- list("token"=Sys.getenv("redcap_api_key"),
                   content='report',
                   format='csv',
                   report_id=report_number,
                   csvDelimiter='',
                   rawOrLabel='label',
                   rawOrLabelHeaders='raw',
                   exportCheckboxLabel='true',
                   returnFormat='csv'
  )
  response <- httr::POST(url, body = formData, encode = "form")
  result <- httr::content(response, show_col_types = FALSE)
}

# Sample ID, species, series and cluster IDs
sample_info <- import_report(samples) %>%
    select(-c(starts_with('redcap_repeat'))) %>%
    filter(isolate_type == "clinical")

# C. albicans MLST types and plot
albicans_sts <- import_report(calb_mlst) %>% 
    select(primary_id, st, aat1a_exact_match, acc1_exact_match, adp1_exact_match,
           mpib_exact_match, sya1_exact_match, vps13_exact_match, zwf1b_exact_match) %>% 
    mutate(across(st:zwf1b_exact_match, as.character)) %>% 
    mutate(concat_alleles=paste(aat1a_exact_match, acc1_exact_match, 
                                adp1_exact_match, mpib_exact_match, 
                                sya1_exact_match, vps13_exact_match, 
                                zwf1b_exact_match, sep = "")) %>% 
    mutate(st = case_when(is.na(st) ~ concat_alleles, .default = st)) %>% 
    inner_join(sample_info %>% select(primary_id, patient_code), by=join_by(primary_id))

albicans_sts$patient_code <- as.character(albicans_sts$patient_code)

albicans_st_summary <- albicans_sts %>% 
    count(st) %>% 
    arrange(desc(n))

st_plot <- ggplot(albicans_sts %>% count(st, patient_code), aes(st, n, fill=patient_code))+
    geom_bar(stat = "identity") +
    scale_fill_manual(values=patient_colors, guide = "none") +
    scale_x_discrete(limits = albicans_st_summary$st) +
    theme_minimal() +
    ylab("Number of isolates") +
    xlab("Sequence type") +
    theme(axis.text.x = element_text(angle = 90))

ggsave("images/Calbicans/Calbicans_MLST_count.png", st_plot, 
       device = png, dpi=300, bg="white",
       width = 8, height = 5, units = "in")


# C. glabrata MLST types and plot
glabrata_sts <- import_report(cglab_mlst) %>% 
    select(primary_id, st, fks_exact_match, leu2_exact_match, nmt1_exact_match,
           trp1_exact_match, ugp1_exact_match, ura3_exact_match) %>% 
    mutate(across(st:ura3_exact_match, as.character)) %>%
    mutate(concat_alleles=paste(fks_exact_match, leu2_exact_match, 
                                nmt1_exact_match, trp1_exact_match, 
                                ugp1_exact_match, ura3_exact_match, sep = "")) %>% 
    mutate(st = case_when(is.na(st) ~ concat_alleles, .default = st)) %>% 
    inner_join(sample_info %>% select(primary_id, patient_code), by=join_by(primary_id))

glabrata_sts$patient_code <- as.character(glabrata_sts$patient_code)

glabrata_st_summary <- glabrata_sts %>% 
    count(st) %>% 
    arrange(desc(n))

st_plot <- ggplot(glabrata_sts %>% count(st, patient_code), aes(st, n, fill=patient_code))+
    geom_bar(stat = "identity") +
    scale_fill_manual(values=patient_colors, guide = "none") +
    scale_x_discrete(limits = glabrata_st_summary$st) +
    theme_minimal() +
    ylab("Number of isolates") +
    xlab("Sequence type") +
    theme(axis.text.x = element_text(angle = 90))

ggsave("images/Cglabrata/Cglabrata_MLST_count.png", st_plot, 
       device = png, dpi=300, bg="white",
       width = 8, height = 5, units = "in")