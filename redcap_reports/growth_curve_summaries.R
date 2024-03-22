## ---------------------------
## Purpose: Summary plots of growth curve metrics
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999) 

## load packages
library(tidyverse)
library(readxl)
library(writexl)
library(patchwork)
library(beeswarm)

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
  url <- "https://redcap.ahc.umn.edu/redcap/api/"
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
    filter(isolate_type == "clinical")

sample_info$genus_species <- str_replace(sample_info$genus_species, "Candida", "C.")

mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, drug, mic_media,mic_date, mic50, eucast_breakpoint, smg, qc_ok) %>% 
    left_join((sample_info %>% 
                   select(primary_id, genus_species, series_id)), 
              by=join_by(primary_id))

mic_info <- mic_info %>% 
    filter(mic_media=="RPMI", !(primary_id %in% c("AMS5122","AMS5123")), qc_ok=="Yes")

mic_info$mic50 <- as.factor(mic_info$mic50)


mic_info$mic50 <- fct_relevel(mic_info$mic50, "0.016", "0.032", "0.064",
                              "0.125", "0.256", "0.5", "1",
                              ">1", "2", "4", "8", "16", "32", ">32")

mic_info <- mic_info %>% 
    group_by(primary_id, drug) %>% 
    filter(mic_date==max(mic_date))


mic_info <- mic_info %>%    
    select(primary_id, drug, mic_date, mic50, eucast_breakpoint, smg) %>% 
    pivot_wider(names_from = drug, values_from = c(mic50, smg, eucast_breakpoint, mic_date))
    

gc <- import_report(growth_curves) %>%
    filter(redcap_repeat_instrument != "NA", 
           !primary_id %in% c("AMS5122", "AMS5123")) %>%
    select(primary_id, redcap_repeat_instance, 
           gc_date, drug_used, 
           gc_temp, gc_time, 
           k, r, t_gen, auc_l, auc_e) 

gc <- gc %>% 
    group_by(primary_id) %>% 
    filter(drug_used=="None", gc_date==max(gc_date), gc_time < 24.1) %>% 
    left_join(sample_info %>% select(primary_id,genus_species,series_id), 
              by = join_by(primary_id))

gc <- gc %>% 
    left_join(mic_info)


################################################################################
# Data wrangling and summary plots.

# for ordering species and colors for consistency
species_count <- sample_info %>%
    group_by(genus_species) %>%
    summarize(species_count=n()) %>%
    arrange(desc(species_count))

species_colors <- species_colors %>% 
    set_names(species_count$genus_species)

gc$genus_species <- factor(gc$genus_species, levels = species_count$genus_species)

k_metric <- ggplot(gc, aes(x=genus_species, y=k, color = genus_species)) + 
    geom_beeswarm(cex = 0.7) +
    theme_bw() +
    scale_color_manual(values=species_colors, guide = "none") +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    #theme(axis.text.x = element_text(angle = 30, 
    #                                 face="italic", 
    #                                 hjust = 1, 
    #                                 vjust = 1)) +
    #xlab("Species") +
    ylab("Carrying capacity, k")

r_metric <- ggplot(gc, aes(x=genus_species, y=r, color = genus_species)) + 
    geom_beeswarm(cex = 0.7) +
    theme_bw() +
    scale_color_manual(values=species_colors, guide = "none") +
    theme(axis.text.x = element_text(angle = 30, 
                                     face="italic", 
                                     hjust = 1, 
                                     vjust = 1)) +
    xlab("Species") +
    ylab("Growth rate, r")

stack_plots <- k_metric/r_metric

auc_metric <- ggplot(gc, aes(x=genus_species, y=auc_e, color=genus_species)) + 
    geom_beeswarm(size = 0.9, cex = 0.7) +
    scale_color_manual(values=species_colors, guide = "none") +
    theme_bw() +
    scale_y_continuous(limits = c(0,25))+
    theme(axis.text.x = element_text(angle = 30, 
                                     face="italic", 
                                     hjust = 1, 
                                     vjust = 1)) +
    xlab("Species") +
    ylab("AUC, 24H")

ggsave(paste0("images/2023_growth_curves/",Sys.Date(),"_AUC_e.png"),
       auc_metric,
       device=png, 
       bg="white", 
       dpi=300, 
       width = 7, 
       height=5.5, 
       units="in")

auc_subset <- gc %>% 
    filter(genus_species %in% c("C. glabrata"), !is.na(eucast_breakpoint_micafungin)) %>% 
    ggplot(aes(x=eucast_breakpoint_micafungin, y=auc_e, color=genus_species)) + 
    geom_beeswarm() +
    scale_color_manual(values=species_colors, guide = "none") +
    theme_bw() +
    scale_y_continuous(limits = c(0,25))+
    #scale_x_discrete(labels =c("FLC Resistant", "FLC Sensitive, dose-dependent"), limits=c("R", "I"))+
    #scale_x_discrete(labels=c("FLC Resistant", "FLC Sensitive"))+
    scale_x_discrete(labels=(c("MCF Resistant", "MCF Sensitive"))) +
    theme(axis.text.x = element_text(#angle = 30, 
                                     #face="italic", 
                                     #hjust = 1, 
                                     #vjust = 1
        ),
          axis.title.x = element_text(face = "italic")) +
    xlab("C. glabrata") +
    ylab("AUC, 24H")

ggsave("images/2023_growth_curves/Cglabrata_auc_by_MCF_eucast.png", auc_subset, device = png, dpi=300, bg="white")

