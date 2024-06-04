## ---------------------------
## Purpose: Summarize and plot MEC isolate (EUCAST) MIC and SMG data by species 
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

## Load packages
library(tidyverse)
library(readxl)
library(patchwork)
library(ggbeeswarm)

# Redcap report IDs
samples <- '58043'
mic_results <- '58044'

species_colors <- c("#88CCEE", "#999933", "#CC6677", "#44AA99", "#117733", "#332288",
                    "#882255","#BBBBBB", "#AA4499", "#DDCC77", "black")

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
################################################################################
# Read in reports.

# Sample ID, species, series and cluster IDs
sample_info <- import_report(samples) %>%
    select(-c(starts_with('redcap_repeat'))) %>%
    filter(isolate_type == "clinical")

sample_info$genus_species <- str_replace(sample_info$genus_species, "Candida", "C.")

# MIC and SMG results
mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, 
           redcap_repeat_instance, 
           drug, 
           mic_media,
           mic_date, 
           mic50, 
           eucast_breakpoint, 
           mean_smg, 
           mean_no_drug_stationary_k,
           sd_no_drug_stationary_k,
           qc_ok) %>% 
    left_join((sample_info %>% 
                   select(primary_id, genus_species, series_id, patient_code)), 
              by=join_by(primary_id))

# Subset to RPMI data with valid control results
# Slice head removes repeated assays, assuming most recent is best (for now, filtering isn't working)
mic_info <- mic_info %>% 
    filter(mic_media %in% c("RPMI liquid", "RPMI agar, Etest"), 
           !(primary_id %in% c("AMS5122","AMS5123","AMS2401")), 
           qc_ok=="Yes") %>% 
    group_by(primary_id, drug) %>% 
    arrange(desc(mic_date)) %>% 
    slice_head()

# For ordering species, colors, drug labels and levels
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

################################################################################
# Plot MIC and SMG per species per drug.
flc <- ggplot(mic_info %>% filter(drug=="fluconazole"), aes(x=mic50)) + 
    geom_bar(aes(fill=genus_species), just = 1) +
    scale_fill_manual(values=species_colors, guide = "none") +
    facet_wrap(.~genus_species, nrow = 1)+
    theme_bw() +
    xlab(NULL)+
    #ylab(NULL)+
    #xlab("MIC50") +
    ylab("Fluconazole\n\n") +
    theme(strip.text = element_text(face = "italic", size =8))+
    theme(panel.grid.major.x = element_blank()) + 
    theme(axis.text.x = element_text(angle=90, vjust = -0.3)) +
    geom_vline(data=filter(mic_info, genus_species %in% c("C. albicans", 
                                                          "C. lusitaniae", 
                                                          "C. parapsilosis",
                                                          "C. tropicalis", 
                                                          "C. dubliniensis",
                                                          "C. kefyr", 
                                                          "C. nivariensis",
                                                          "C. orthopsilosis",
                                                          "C. utilis")), 
               aes(xintercept = "4"), linetype =2) +
    geom_vline(data = filter(mic_info, genus_species=="C. glabrata"),
               aes(xintercept = "16"), linetype = 2)

mcf <- ggplot(mic_info %>% filter(drug=="micafungin"), aes(x=mic50)) + 
    geom_bar(aes(fill=genus_species), just = 1) +
    scale_fill_manual(values=species_colors, guide = "none") +
    facet_wrap(.~genus_species, nrow = 1)+
    theme_bw() +
    xlab(NULL)+
    ylab("Micafungin\n\nCount of isolates")+
    theme(strip.text = element_blank()) +
    theme(axis.text.x = element_text(angle=90, vjust = -0.4)) +
    theme(panel.grid.major.x = element_blank()) + 
    geom_vline(data=filter(mic_info, genus_species=="C. albicans"), 
               aes(xintercept = "0.016"),linetype = 2) +
    geom_vline(data = filter(mic_info, genus_species=="C. glabrata"),
               aes(xintercept = "0.032"), linetype = 2) +
    geom_vline(data = filter(mic_info, genus_species=="C. parapsilosis"),
               aes(xintercept = "2"), linetype = 2) 

amb <- ggplot(mic_info %>% filter(drug=="amphotericin B"), aes(x=mic50)) + 
    geom_bar(aes(fill=genus_species), just = 1) +
    scale_fill_manual(values=species_colors, guide = "none") +
    facet_wrap(.~genus_species, nrow = 1)+
    theme_bw() +
    xlab("\nMIC, \u03BCg/mL") +
    ylab("Amphotericin B\n\n") +
    #ylab(NULL) +
    theme(strip.text = element_blank())+
    theme(axis.text.x = element_text(angle=90, vjust=-0.5)) +
    theme(panel.grid.major.x = element_blank()) + 
    geom_vline(data=filter(mic_info, genus_species %in% c("C. albicans", 
                                                          "C. glabrata", 
                                                          "C. parapsilosis", 
                                                          "C. krusei", 
                                                          "C. tropicalis", 
                                                          "C. dubliniensis")), 
               aes(xintercept = "1"), linetype =2)

full_mic_plot <- flc/mcf/amb

ggsave(paste0("images/2023_MICs/",Sys.Date(),"_MEC_MIC_summary.png"),
       full_mic_plot,
       device=png, 
       bg="white", 
       dpi=300, 
       width=13.5, 
       height = 7.5, 
       units = "in")


smg <- ggplot(mic_info, aes(x=genus_species, y=smg, color = genus_species)) + 
    geom_beeswarm(size=2, cex=0.8) +#geom_violin() +
    #geom_point(data = filter(mic_info, genus_species=="C. nivariensis")) +
    facet_grid(factor(drug, 
                      levels=c("fluconazole", "micafungin", "amphotericin B")) ~ ., 
               labeller = drugs,
               switch = "y") +
    theme_bw() +
    scale_color_manual(values=species_colors, guide = "none") +
    theme(axis.text.x = element_text(angle = 30, 
                                     face="italic", 
                                     hjust = 1, 
                                     vjust = 1)) +
    xlab(NULL) +
    ylab("SMG")

ggsave(paste0("images/2023_MICs/",Sys.Date(),"_MEC_SMG_summary.png"),
       smg,
       device=png, 
       bg="white", 
       dpi=300, 
       width = 11.5, 
       height=7, 
       units="in")
