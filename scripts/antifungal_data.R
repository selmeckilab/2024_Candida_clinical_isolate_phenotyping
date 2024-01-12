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
af_spreadsheet <- "antifungal_relative_days_by_sample.xlsx"
breakpoint_file <- "Candida_eucast_breakpoints.xlsx"

# redcap report IDs
samples <- '58043'
mic_results <- '58044'
growth_curves <- '58045'
genes <- '58048'

token <- '' # no gh

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
    left_join((sample_info %>% select(primary_id, genus_species)))

# MIC and SMG results
mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, drug, mic_date, mic50, eucast_breakpoint, smg) %>% 
    left_join((sample_info %>% 
                   select(primary_id, genus_species, series_id)), 
              by=join_by(primary_id))

################################################################################
# Data wrangling and summary plots.

mic_info <- mic_info %>% 
    filter(!(is.na(eucast_breakpoint)), !(primary_id %in% c("AMS5122","AMS5123")))

# for ordering species and colors for consistency
species_count <- sample_info %>%
    group_by(genus_species) %>%
    summarize(species_count=n()) %>%
    arrange(desc(species_count))

species_colors <- species_colors %>% 
    set_names(species_count$genus_species)

drugs <- as_labeller(c(`amphotericin B` = "Amphotericin B", 
                       fluconazole="Fluconazole", 
                       micafungin="Micafungin"))

mic_info$mic50 <- factor(mic_info$mic50, levels=c("0.016", "0.032", "0.064", "0.125",
                                                  "0.256", "0.5", "1", ">1", "2", 
                                                  "4", "8", "16", "32", ">32"))

mic_info$genus_species <- factor(mic_info$genus_species, levels = species_count$genus_species)

# Plot MIC and SMG per species per drug.

amb <- ggplot(mic_info %>% filter(drug=="amphotericin B"), aes(x=mic50)) + 
    geom_bar(aes(fill=genus_species), just = 1) +
    scale_fill_manual(values=species_colors, guide = "none") +
    facet_wrap(.~genus_species, nrow = 1)+
    theme_bw() +
    xlab("\nMIC") +
    ylab("Amphotericin B\n") +
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

mcf <- ggplot(mic_info %>% filter(drug=="micafungin"), aes(x=mic50)) + 
    geom_bar(aes(fill=genus_species), just = 1) +
    scale_fill_manual(values=species_colors, guide = "none") +
    facet_wrap(.~genus_species, nrow = 1)+
    theme_bw() +
    xlab(NULL)+
    ylab("Micafungin\nNumber of isolates")+
    theme(strip.text = element_blank()) +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5)) +
    theme(panel.grid.major.x = element_blank()) + 
    geom_vline(data=filter(mic_info, genus_species=="Candida albicans"), 
               aes(xintercept = "0.016"),linetype = 2) +
    geom_vline(data = filter(mic_info, genus_species=="Candida glabrata"),
               aes(xintercept = "0.032"), linetype = 2) +
    geom_vline(data = filter(mic_info, genus_species=="Candida parapsilosis"),
               aes(xintercept = "2"), linetype = 2) 

flc <- ggplot(mic_info %>% filter(drug=="fluconazole"), aes(x=mic50)) + 
    geom_bar(aes(fill=genus_species), just = 1) +
    scale_fill_manual(values=species_colors, guide = "none") +
    facet_wrap(.~genus_species, nrow = 1)+
    theme_bw() +
    xlab(NULL)+
    #ylab(NULL)+
    #xlab("MIC50") +
    ylab("Fluconazole\n") +
    theme(strip.text = element_text(face = "italic"))+
    theme(panel.grid.major.x = element_blank()) + 
    theme(axis.text.x = element_text(angle=90, vjust = -0.5)) +
    geom_vline(data=filter(mic_info, genus_species %in% c("Candida albicans", 
                                                          "Candida lusitaniae", 
                                                          "Candida parapsilosis",
                                                          "Candida tropicalis", 
                                                          "Candida dubliniensis",
                                                          "Candida kefyr", 
                                                          "Candida nivariensis",
                                                          "Candida orthopsilosis")), 
               aes(xintercept = "4"), linetype =2) +
    geom_vline(data = filter(mic_info, genus_species=="Candida glabrata"),
               aes(xintercept = "16"), linetype = 2)

stack <- flc/mcf/amb

ggsave(paste0("images/2023_MICs/",Sys.Date(),"_MEC_MIC_summary.png"),
       stack,
       device=png, 
       bg="white", 
       dpi=300, 
       width=11.5, 
       height = 7, 
       units = "in")

smg <- ggplot(mic_info, aes(x=genus_species, y=smg)) + 
    geom_boxplot(aes(fill=genus_species), 
                 outlier.colour = NA, 
                 width = 0.4, 
                 alpha = 0.8) +
    geom_dotplot(binaxis = "y", 
                 binwidth = 0.015,
                 stackdir = "center", 
                 color = "grey36", 
                 fill="grey36") + 
    facet_grid(drug ~ ., labeller = drugs) +
    theme_bw() +
    scale_fill_manual(values=species_colors, guide = "none") +
    theme(axis.text.x = element_text(angle = 30, 
                                     face="italic", 
                                     hjust = 1, 
                                     vjust = 1)) +
    xlab("SMG by species") +
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
# Resistance/variant data wrangling

# Filter C. glabrata data for resistant isolates and resistant-only SNPs
mcf_res <- mic_info %>% filter(drug=="micafungin", eucast_breakpoint=="R")
mcf_sensitive <- mic_info %>% filter(drug=="micafungin", eucast_breakpoint=="S", genus_species=="Candida glabrata")

glabrata_vars <- gene_vars %>% filter(genus_species=="Candida glabrata")

glab_mcf <- mic_info %>% filter(genus_species=="Candida glabrata", drug=="micafungin")
mcf_res_vars <- glabrata_vars %>% filter(primary_id %in% mcf_res$primary_id) %>% select(gene, protein_change)
mcf_sens_vars <- glabrata_vars %>% filter(primary_id %in% mcf_sensitive$primary_id) %>% select(gene, protein_change)
possible_mcf <- mcf_res_vars %>% anti_join(mcf_sens_vars)

glab_flc <- mic_info %>% filter(genus_species=="Candida glabrata", drug=="fluconazole", !(is.na(eucast_breakpoint)))
flc_res <- glab_flc %>% filter(eucast_breakpoint=="R")
flc_int <- glab_flc %>% filter(eucast_breakpoint=="I")
flc_sens_vars <- glabrata_vars %>% filter(primary_id %in% flc_int$primary_id) %>% select(gene, protein_change)
flc_res_vars <- glabrata_vars %>% filter(primary_id %in% flc_res$primary_id) %>% select(gene, protein_change)

possible_mcf <- unique(mcf_res_vars %>% anti_join(mcf_sens_vars))
possible_mcf_joined <- gene_vars %>% inner_join(possible_mcf)

possible_flc <- flc_res_vars %>% anti_join(flc_sens_vars) 
possible_flc <- unique(flc_res_vars %>% anti_join(flc_sens_vars))
possible_flc_joined <- gene_vars %>% inner_join(possible_flc)

multi_res_vars <- glabrata_vars %>% filter(primary_id %in% mcf_res$primary_id, primary_id %in% flc_res$primary_id) %>% select(gene, protein_change)
possible_multi <- unique(multi_res_vars %>% anti_join(mcf_sens_vars)) %>% anti_join(flc_sens_vars)
possible_multi_joined <- gene_vars %>% inner_join(possible_multi)
################################################################################
# Drug exposure data

# Process de-identified antifungal relative timeframe data
antifungals <- read_xlsx(af_spreadsheet)
antifungals <- antifungals %>% 
    inner_join((sample_info %>% 
                    select(primary_id, genus_species, series_id, relative_days)), 
               by=join_by(sample==primary_id))

remote_exposure <- antifungals %>% 
    #group_by(series_id) %>% 
    filter(relative_start < -1, relative_days == 0)

no_results <- sample_info %>% filter(!(primary_id %in% antifungals$sample)) %>% 
    select(primary_id, genus_species, relative_days, series_id, secondary_id)
