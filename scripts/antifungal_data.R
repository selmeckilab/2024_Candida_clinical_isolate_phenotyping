## ---------------------------
## Script name: antifungal_data.R
##
## Purpose of script: Analyze redcap and antifungal exposure data jointly
##
## Author: Nancy Scott
##
## Date Created: 2024-01-03
##
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

## load packages
library(tidyverse)
library(readxl)
library(writexl)
library(patchwork)

# local antifungal data
af_spreadsheet <- "data/metadata/MEC_antifungal_history.xlsx"
breakpoint_file <- "Candida_eucast_breakpoints.xlsx"

# redcap report IDs
samples <- '58043'
mic_results <- '58044'
growth_curves <- '58045'
genes <- '58048'

token <- ''

species_colors <- c("#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
                    "#882255","#BBBBBB", "#AA4499", "#DDCC77", "black")

# function to import report from redcap
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
################################################################################
# Read in reports.

# Sample ID, species, series and cluster IDs
sample_info <- import_report(samples) %>%
    select(-c(starts_with('redcap_repeat'))) %>%
    filter(!primary_id %in% c("MEC103", "MEC113")) %>%
    filter(isolate_type == "clinical")

# SNP data
gene_vars <- import_report(genes) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, gene, protein_change, alt_freq) %>%
    left_join((sample_info %>% select(primary_id, genus_species))) %>% 
    unite("gene_pos", gene:protein_change)

# MIC and SMG results
mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, drug, mic_media,mic_date, mic50, eucast_breakpoint, smg) %>% 
    left_join((sample_info %>% 
                   select(primary_id, genus_species, series_id)), 
              by=join_by(primary_id))

################################################################################
# Data wrangling and summary plots.

mic_info <- mic_info %>% 
    filter(mic_media=="RPMI", !(primary_id %in% c("AMS5122","AMS5123"))) %>%
    filter(!(primary_id %in% c("MEC216", "MEC217", "MEC218", "MIC219") & mic_date==as.Date("2023-08-03"))) %>% 
    inner_join(sample_info %>% select(primary_id, patient_code))

# for ordering species and colors for consistency
species_count <- sample_info %>%
    group_by(genus_species) %>%
    summarize(species_count=n()) %>%
    arrange(desc(species_count))

species_colors <- species_colors %>% 
    set_names(species_count$genus_species)

drugs <- as_labeller(c(fluconazole="Fluconazole", 
                       micafungin="Micafungin",
                       `amphotericin B` = "Amphotericin B"))
 
mic_info$mic50 <- factor(mic_info$mic50, levels=c("0.016", "0.032", "0.064", "0.125",
                                                  "0.256", "0.5", "1", ">1", "2", 
                                                  "4", "8", "16", "32", ">32"))

mic_info$genus_species <- factor(mic_info$genus_species, levels = species_count$genus_species)

# Plot MIC and SMG per species per drug.
flc <- ggplot(mic_info %>% filter(drug=="fluconazole"), aes(x=mic50)) + 
    geom_bar(aes(fill=genus_species), just = 1) +
    scale_fill_manual(values=species_colors, guide = "none") +
    facet_wrap(.~genus_species, nrow = 1)+
    theme_bw() +
    xlab(NULL)+
    #ylab(NULL)+
    #xlab("MIC50") +
    ylab("\nFluconazole\n") +
    theme(strip.text = element_text(face = "italic", size =8))+
    theme(panel.grid.major.x = element_blank()) + 
    theme(axis.text.x = element_text(angle=90, vjust = -0.3)) +
    geom_vline(data=filter(mic_info, genus_species %in% c("Candida albicans", 
                                                          "Candida lusitaniae", 
                                                          "Candida parapsilosis",
                                                          "Candida tropicalis", 
                                                          "Candida dubliniensis",
                                                          "Candida kefyr", 
                                                          "Candida nivariensis",
                                                          "Candida orthopsilosis",
                                                          "Candida utilis")), 
               aes(xintercept = "4"), linetype =2) +
    geom_vline(data = filter(mic_info, genus_species=="Candida glabrata"),
               aes(xintercept = "16"), linetype = 2)

mcf <- ggplot(mic_info %>% filter(drug=="micafungin"), aes(x=mic50)) + 
    geom_bar(aes(fill=genus_species), just = 1) +
    scale_fill_manual(values=species_colors, guide = "none") +
    facet_wrap(.~genus_species, nrow = 1)+
    theme_bw() +
    xlab(NULL)+
    ylab("Number of isolates\nMicafungin\n")+
    theme(strip.text = element_blank()) +
    theme(axis.text.x = element_text(angle=90, vjust = -0.4)) +
    theme(panel.grid.major.x = element_blank()) + 
    geom_vline(data=filter(mic_info, genus_species=="Candida albicans"), 
               aes(xintercept = "0.016"),linetype = 2) +
    geom_vline(data = filter(mic_info, genus_species=="Candida glabrata"),
               aes(xintercept = "0.032"), linetype = 2) +
    geom_vline(data = filter(mic_info, genus_species=="Candida parapsilosis"),
               aes(xintercept = "2"), linetype = 2) 

amb <- ggplot(mic_info %>% filter(drug=="amphotericin B"), aes(x=mic50)) + 
    geom_bar(aes(fill=genus_species), just = 1) +
    scale_fill_manual(values=species_colors, guide = "none") +
    facet_wrap(.~genus_species, nrow = 1)+
    theme_bw() +
    xlab("\nMIC value") +
    ylab("\nAmphotericin B\n") +
    #ylab(NULL) +
    theme(strip.text = element_blank())+
    theme(axis.text.x = element_text(angle=90, vjust=-0.5)) +
    theme(panel.grid.major.x = element_blank()) + 
    geom_vline(data=filter(mic_info, genus_species %in% c("Candida albicans", 
                                                          "Candida glabrata", 
                                                          "Candida parapsilosis", 
                                                          "Candida krusei", 
                                                          "Candida tropicalis", 
                                                          "Candida dubliniensis")), 
               aes(xintercept = "1"), linetype =2)

stack_plot <- flc/mcf/amb

ggsave(paste0("images/2023_MICs/",Sys.Date(),"_MEC_MIC_summary.png"),
       stack_plot,
       device=png, 
       bg="white", 
       dpi=300, 
       width=13.5, 
       height = 7.5, 
       units = "in")


smg <- ggplot(mic_info, aes(x=genus_species, y=smg, fill = genus_species)) + 
    geom_violin() +
    geom_point(data = filter(mic_info, genus_species=="Candida nivariensis")) +
    facet_grid(factor(drug, 
                      levels=c("fluconazole", "micafungin", "amphotericin B")) ~ ., 
               labeller = drugs,
               switch = "y") +
    theme_bw() +
    scale_fill_manual(values=species_colors, guide = "none") +
    theme(axis.text.x = element_text(angle = 30, 
                                     face="italic", 
                                     hjust = 1, 
                                     vjust = 1)) +
    xlab("Species") +
    ylab("SMG")

ggsave(paste0("images/2023_MICs/",Sys.Date(),"_MEC_SMG_summary.png"),
       smg,
       device=png, 
       bg="white", 
       dpi=300, 
       width = 11.5, 
       height=7, 
       units="in")

################################################################################
# Drug exposure data

# Process de-identified antifungal relative timeframe data
antifungals <- read_xlsx(af_spreadsheet)

remote_exposure <- antifungals %>% 
    filter(relative_start_days <= 0, relative_collection_day == 0)

no_results <- sample_info %>% filter(!(primary_id %in% antifungals$primary_id)) %>% 
    select(primary_id, genus_species, relative_days, patient_code, series_id, secondary_id)

################################################################################
# Resistance/variant data wrangling

mic_summary <- mic_info %>%
    group_by(drug) %>%
    count(genus_species)

resistant_isolates <- mic_info %>% 
    group_by(drug, genus_species) %>% 
    filter(eucast_breakpoint=="R")

res_summary <- resistant_isolates %>% 
    group_by(drug, genus_species) %>%
    count(patient_code)

res_frequency <- mic_summary %>% 
    inner_join(res_summary, by=join_by(drug,genus_species)) #%>% 
    #mutate(resistant = n.y/n.x)

patient_counts <- mic_info %>%
    group_by(drug, genus_species) %>%
    summarize(patients=length(unique(patient_code)))

patient_res <- resistant_isolates %>% 
    group_by(drug, genus_species) %>%
    summarize(patient_res=length(unique(patient_code)))

patient_freq <- patient_counts %>% 
    inner_join(patient_res, by=join_by(drug, genus_species)) %>% 
    mutate(patient_res_freq = patient_res/patients)
    
non_driver_vars <- gene_vars %>% 
    filter(!primary_id %in% resistant_isolates$primary_id)

possible_drivers <- gene_vars %>% 
    filter(!gene_pos %in% non_driver_vars$gene_pos)
