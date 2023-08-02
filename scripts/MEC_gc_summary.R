## ---------------------------
## Script name: MEC_gc_summary.R
##
## Purpose of script: use growthcurver for gc metric calculation 
##
## Author: Nancy Scott
##
## Date Created: 2023-08-02
##
## Email: scot0854@umn.edu
## ---------------------------
## Notes: Adapted from C. lusitaniae growth curve scripts
##   
## ---------------------------
options(scipen = 999) # To view outputs in non-scientific notation
## ---------------------------
## load packages
library(tidyverse)
library(readxl)
library(writexl)
library(RCurl)
library(digest)
library(jsonlite)
library(ggplot2)
library(gtools)
library(ggpubr)
library(ggthemes)
library(patchwork)
################################################################################
## Get spreadsheet cleaning and gc metric functions
source("gc_functions.R")
###############################################################################
## Input vars
plate_reader_file <- "data/growth_curve/2023-07-19_MEC_GC30.xlsx"
gc_date <- "2023-07-19"

api_token <- ""
api_url <-  "https://redcap.ahc.umn.edu/api/"

################################################################################
## Read in plate, get sample metadata and check fit of logistic curve

plate_od <- clean_growthcurver(plate_reader_file) #%>%
    #filter(time < 23.48) # make sure all plates are same time span

samples <- samples(plate_reader_file)

plate_data <- SummarizeGrowthByPlate(plate_od)
colnames(plate_data)[colnames(plate_data)=="sample"] <- 'well'


###############################################################################
# tidy, add metadata and summarize values
plate_summary <- plate_data %>%
    inner_join(samples, by = "well") %>%
    group_by(primary_id) %>%
    filter(!primary_id %in% c("water", "blank", "NA", "ypad")) %>%
    summarize(mean_k = round(mean(k), digits = 3), 
              mean_n0 = round(mean(n0), digits = 3),
              mean_r = round(mean(r), digits = 3),
              mean_tmid = round(mean(t_mid), digits = 3),
              mean_tgen = round(mean(t_gen), digits = 3),
              mean_aucl = round(mean(auc_l), digits = 3), 
              mean_auce = round(mean(auc_e), digits = 3))

###############################################################################
# redcap upload work 

redcap_media <- case_when(toupper(samples$media[1])== "YPAD" ~ 1,
                         toupper(samples$media[1])== "RPMI" ~ 3)

redcap_temp <- case_when(samples$temp[1]== 30 ~ 0,
                         samples$temp[1] == 37 ~1)

redcap_reader <- 0 # serial number 1803065

for(i in 2:length(plate_summary$primary_id)){
    
    record <- c(
        primary_id = plate_summary$primary_id[i],
        redcap_repeat_instrument = "growth_curve_data",
        redcap_repeat_instance = "new",
        gc_date = gc_date,
        gc_time = round(max(plate_od$time), digits = 2),
        gc_temp = redcap_temp,
        gc_media = redcap_media,
        drug_used = samples$drug[1],
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
    
    result_data <- toJSON(list(as.list(record)), auto_unbox=TRUE)
    
    formData <- list("token"=api_token,
                     content='record',
                     action='import',
                     format='json',
                     type='flat',
                     overwriteBehavior='normal',
                     forceAutoNumber='false',
                     data=result_data,
                     returnContent='count',
                     returnFormat='json'
    )
   response <- httr::POST(api_url, body = formData, encode = "form")
   result <- httr::content(response)
   print(result)
}

###############################################################################
# basic GC plot

input_data <- growth_curve(plate_reader_file)
full_gc_data <- link_metadata(input_data, samples, plate = 1)

blank_od <- full_gc_data %>%
    filter(primary_id == "ypad") %>%
    summarise(mean_blank = mean(OD))

gc_norm <- full_gc_data %>%
    filter(!primary_id %in% c("ypad", "water", "NA")) %>%
    mutate(OD = OD - blank_od$mean_blank)

gc_mean <- gc_norm %>%
    group_by(primary_id, time) %>%
    mutate(mean_OD = mean(OD, na.rm = TRUE), sd_OD=sd(OD, na.rm=TRUE), se_OD=sd_OD/sqrt(length(OD))) %>%
    mutate(primary_id = as.factor(primary_id)) %>%
    mutate(primary_id = fct_relevel(primary_id, mixedsort))

full_plot <- gc_mean %>%
    ggplot(aes(x = time, y = mean_OD, color = primary_id)) +
    geom_point(show.legend = TRUE) +
    geom_errorbar(aes(ymax=(mean_OD + sd_OD), ymin=(mean_OD - sd_OD)), alpha=0.2, show.legend = FALSE) +
    scale_y_continuous(breaks = c(0,0.5,1.0,1.5)) +
    scale_x_continuous(breaks = c(0, 12, 24, 36, 48),
                       labels = c ("0", "12", "24", "36", "48")) +
    scale_color_tableau(palette = "Tableau 20", limits=levels(gc_mean$primary_id)) +
    xlab("Time (hours)") +
    ylab("Mean OD600") +
    labs(colour = "Isolate ID", title = "Growth in YPAD at 30C")

ggsave(paste0("images/",gc_date,"_combined_MEC_GC30.tiff"), full_plot, width = 9, height = 5.8, units = "in")

facet_plot <- gc_mean %>%
    ggplot(aes(x = time, y = mean_OD, color = primary_id)) +
    geom_point(show.legend = TRUE) +
    geom_errorbar(aes(ymax=(mean_OD + sd_OD), ymin=(mean_OD - sd_OD)), alpha=0.2, show.legend = FALSE) +
    facet_wrap(~as.factor(primary_id), nrow = 2)+
    scale_y_continuous(breaks = c(0,0.5,1.0,1.5)) +
    scale_x_continuous(breaks = c(0, 12, 24, 36, 48),
                       labels = c ("0", "12", "24", "36", "48")) +
    scale_color_tableau(palette = "Tableau 20", limits=levels(gc_mean$primary_id)) +
    xlab("Time (hours)") +
    ylab("Mean OD600") +
    labs(colour = "Isolate ID", title = "Growth in YPAD at 30C")

ggsave(paste0("images/",gc_date,"_faceted_MEC_GC30.tiff"), facet_plot, width = 11, height = 8.5, units = "in")
