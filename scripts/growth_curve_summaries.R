## ---------------------------
## Script name: growth_curve_summaries.R
##
## Purpose of script: Summary plots of growth curve metrics
##
## Author: Nancy Scott
##
## Date Created: 2024-01-15
##
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999) 

## load packages
library(tidyverse)
library(readxl)
library(writexl)
library(patchwork)

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
    filter(!primary_id %in% c("MEC103", "MEC113")) %>%
    filter(isolate_type == "clinical")

sample_info$genus_species <- str_replace(sample_info$genus_species, "Candida", "C.")

mic_info <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, drug, mic_media,mic_date, mic50, eucast_breakpoint, smg) %>% 
    left_join((sample_info %>% 
                   select(primary_id, genus_species, series_id)), 
              by=join_by(primary_id))

mic_info <- mic_info %>% 
    filter(mic_media=="RPMI", !(primary_id %in% c("AMS5122","AMS5123")))

mic_info$mic50 <- as.factor(mic_info$mic50)


mic_info$mic50 <- fct_relevel(mic_info$mic50, "0.016", "0.032", "0.064",
                              "0.125", "0.256", "0.5", "1",
                              ">1", "2", "4", "8", "16", "32", ">32")

mic_info <- mic_info %>% 
    group_by(primary_id) %>% 
    filter(mic_date==max(mic_date))


mic_info <- mic_info %>%    
    select(primary_id, drug, mic_date, mic50, eucast_breakpoint, smg) %>% 
    pivot_wider(names_from = drug, values_from = c(mic50, smg, eucast_breakpoint))
    

gc <- import_report(growth_curves) %>%
    filter(redcap_repeat_instrument != "NA", 
           !primary_id %in% c("AMS5122", "AMS5123")) %>%
    select(primary_id, redcap_repeat_instance, 
           gc_date, drug_used, 
           gc_temp, gc_time, 
           k, r, t_gen, auc_l) 

gc <- gc %>% 
    filter(drug_used=="None", gc_date > as.Date("2022-12-31")) %>% 
    left_join(sample_info %>% select(primary_id,genus_species,series_id), 
              by = join_by(primary_id))

gc <- gc %>% group_by(primary_id) %>% 
    filter(redcap_repeat_instance ==max(redcap_repeat_instance))

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

k_metric <- ggplot(gc, aes(x=genus_species, y=k)) + 
    geom_boxplot(aes(fill=genus_species), 
                 outlier.colour = NA, 
                 width = 0.4, 
                 alpha = 0.8) +
    geom_dotplot(binaxis = "y", 
                 binwidth = 0.015,
                 stackdir = "center", 
                 color = "grey36", 
                 fill="grey36") + 
    theme_bw() +
    scale_fill_manual(values=species_colors, guide = "none") +
    theme(axis.text.x = element_blank()) +
    theme(axis.title.x = element_blank()) +
    #theme(axis.text.x = element_text(angle = 30, 
    #                                 face="italic", 
    #                                 hjust = 1, 
    #                                 vjust = 1)) +
    #xlab("Species") +
    ylab("Carrying capacity, k")

r_metric <- ggplot(gc, aes(x=genus_species, y=r)) + 
    geom_boxplot(aes(fill=genus_species), 
                 outlier.colour = NA, 
                 width = 0.4, 
                 alpha = 0.8) +
    geom_dotplot(binaxis = "y", 
                 binwidth = 0.015,
                 stackdir = "center", 
                 color = "grey36", 
                 fill="grey36") + 
    theme_bw() +
    scale_fill_manual(values=species_colors, guide = "none") +
    theme(axis.text.x = element_text(angle = 30, 
                                     face="italic", 
                                     hjust = 1, 
                                     vjust = 1)) +
    xlab("Species") +
    ylab("Growth rate, r")

stack_plots <- k_metric/r_metric

tgen_metric <- ggplot(gc, aes(x=genus_species, y=t_gen, fill=genus_species)) + 
    geom_violin() +
    scale_fill_manual(values=species_colors, guide = "none") +
    geom_point(data = filter(gc, genus_species=="C. nivariensis")) +
    theme_bw() +
    scale_y_continuous(limits = c(0,3.5), breaks = c(0,1,2,3))+
    theme(axis.text.x = element_text(angle = 30, 
                                     face="italic", 
                                     hjust = 1, 
                                     vjust = 1)) +
    xlab("Species") +
    ylab("Doubling time (hr)")

ggsave(paste0("images/2023_growth_curves/",Sys.Date(),"_doubling_time.png"),
       tgen_metric,
       device=png, 
       bg="white", 
       dpi=300, 
       width = 7, 
       height=5, 
       units="in")

tgen_subset <- gc %>% 
    filter(genus_species %in% c("C. albicans"), !is.na(eucast_breakpoint_fluconazole)) %>% 
    ggplot(aes(x=eucast_breakpoint_fluconazole, y=t_gen, fill=genus_species)) + 
    geom_violin() +
    scale_fill_manual(values=species_colors, guide = "none") +
    geom_dotplot(binaxis = "y", 
                 binwidth = 0.015,
                 stackdir = "center", 
                 color = "grey36", 
                 fill="grey36") + 
    theme_bw() +
    scale_y_continuous(limits = c(0,3.5), breaks = c(0,1,2,3))+
    #scale_x_discrete(labels =c("FLC Resistant", "FLC Intermediate"), limits=c("R", "I"))+
    scale_x_discrete(labels=c("FLC Resistant", "FLC Sensitive"))+
    theme(axis.text.x = element_text(#angle = 30, 
                                     #face="italic", 
                                     #hjust = 1, 
                                     #vjust = 1
        ),
          axis.title.x = element_text(face = "italic")) +
    xlab("C. albicans") +
    ylab("Doubling time (hr)")

ggsave("images/2023_growth_curves/Calbicans_tgen_by_FLC_eucast.png", device = png, dpi=300, bg="white")

