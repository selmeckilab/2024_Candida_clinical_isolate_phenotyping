## ---------------------------
## Script name: redcap_reports.R
##
## Purpose of script: Pull, filter, subset REDCap reports
##
## Author: Nancy Scott
##
## Date Created: 2023-11-28
##
## Email: scot0854@umn.edu
## ---------------------------
## load packages
library(tidyverse)
library(writexl)
library(paletteer)

# redcap report IDs
samples <- '58043'
mic_results <- '58044'
growth_curves <- '58045'
chef_data <- '58046'
genes <- '58048'
avail_seq_data <- '58050'
calb_mlst <- '58053'
cglab_mlst <- '58052'

token <- ''

patient_colors <- c(paletteer_d("ggsci::default_igv"),
                    "#838b8b",
                    "black",
                    "#cdcdb4",
                    "#155F83")

# Function to import report from redcap
import_report <- function(report_number) {
  url <- "https://redcap.ahc.umn.edu/api/"
  formData <- list("token"=token,
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
    filter(!primary_id %in% c("MEC103", "MEC113")) %>%
    filter(isolate_type == "clinical")

################################################################################
# CHEF gel results
chef_done <- import_report(chef_data) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gel_date, blot_prepared) %>%
    pivot_wider(names_from = "redcap_repeat_instance",
                values_from = c("gel_date", "blot_prepared"),
                names_vary = "slowest")

# No CHEF results yet
todo <- sample_info %>%
    filter(isolate_type == "clinical") %>%
    anti_join(chef_done)

# Export chef-todo list
#write_xlsx(todo,paste0(Sys.Date(),"_CHEF_todo.xlsx"))

################################################################################
# MIC and SMG results
mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA", mic_media == "RPMI") %>%
    select(primary_id, redcap_repeat_instance, drug, mic_date, mic50, eucast_breakpoint, smg) %>% 
    left_join((sample_info %>% select(primary_id, genus_species, series_id)), by=join_by(primary_id)) 

# MIC to-do list
mic_to_do <- sample_info %>% 
    filter(isolate_type == "clinical") %>%
    anti_join(mic_info)

################################################################################
# Growth curve results
gc <- import_report(growth_curves) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gc_date, drug_used, gc_temp, gc_time, k, r, t_gen, auc_l)

################################################################################
# MSI location if sequencing data exists
seq_info <- import_report(avail_seq_data) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, msi_path_r1, msi_path_r2, msi_long_read_path) %>%
    pivot_wider(names_from = "redcap_repeat_instance",
                values_from = c("msi_path_r1", "msi_path_r2", "msi_long_read_path"),
                names_vary = "slowest")

################################################################################
# SNP data
gene_vars <- import_report(genes) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gene, protein_change, alt_freq) %>%
    left_join((sample_info %>% select(primary_id, genus_species)))

################################################################################
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

####################################################################
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
