## ---------------------------
## Script name: clin_isolate_MIC.R
##
## Purpose of script: Calculate MIC data for subset of MCF-evolved isolates
##
## Author: Nancy Scott
##
## Date Created: 2023-01-07
##
## Email: scot0854@umn.edu
## ---------------------------
## Notes:
## OD readings are read from one sheet and pivoted to "long" format, 
## and then metadata in long format is read from another sheet and concatenated
## ---------------------------
options(scipen = 999) # To view outputs in non-scientific notation
## ---------------------------
## load packages
library(readxl)
library(tidyverse)
library(gtools)
library(RCurl)
library(digest)
library(jsonlite)
library(patchwork)
## ---------------------------
## variables and functions
source("~/umn/mic_data/scripts/MIC_calc_functions.R")
source("~/umn/mic_data/scripts/MIC_heatmap.R")

api_token <- ""
input_drug <- "amb"
mic_date <- "2023-10-31"
mic_spreadsheet <-"data/MIC/2023-10-31_EW_MIC24_RPMI35.xlsx"
smg_spreadsheet <- "data/MIC/2023-11-01_EW_SMG48_RPMI35.xlsx"
#strains <- c("AMS5122",
#             "MEC161", "MEC151", "MEC150", "MEC149",
#             "MEC148", "MEC147", "MEC146", "MEC088",
#             "MEC118", "MEC117", "MEC116")

#exclude <- c("MEC219", "MEC218", "MEC217", "MEC216")
control_strains <- c("AMS5123", "AMS5122", "AMS2401")
api_url <-  "https://redcap.ahc.umn.edu/api/"

################################################################################
# MIC
# Load metadata
meta.frame <- read_excel(mic_spreadsheet, 
                         sheet = 2)
meta.frame$drug <- c(toupper(input_drug), rep(NA, times=11))
# Make sure at least one column is length 12
#stopifnot( "Neither strain or concentration are right number of cols" = 12 %in% c(length(meta.frame$strain[!is.na(meta.frame$strain)]), length(meta.frame$concentration[!is.na(meta.frame$concentration)])))

meta.frame$concentration <- case_when(meta.frame$drug[1] == "FLC" ~ c(0,0.5,1,2,4,8,16,32,NA,NA,NA,NA),
                                      meta.frame$drug[1] %in% c("MCF","AMB") ~ c(0,0.016,0.032,0.064,0.125,0.256,0.5,1,NA,NA,NA,NA))
#set vars for use below
meta_names = case_when(length(meta.frame$strain[!is.na(meta.frame$strain)]) == 12 ~ "strain",
                       length(meta.frame$concentration[!is.na(meta.frame$concentration)]) == 12 ~ "concentration")
meta_col = case_when(length(meta.frame$strain[!is.na(meta.frame$strain)]) != 12 ~ "strain",
                     length(meta.frame$concentration[!is.na(meta.frame$concentration)]) != 12 ~ "concentration")
    
cutoff <- case_when(meta.frame$drug[1] == "FLC" ~ 0.5,
                    meta.frame$drug[1] == "MCF" ~ 0.5,
                    meta.frame$drug[1] == "AMB" ~ 0.1)

max_concentration <- max(meta.frame$concentration, na.rm = TRUE)

x_axis_angle <- case_when(meta.frame$drug[1]== "FLC" ~ 0,
                          meta.frame$drug[1] == "MCF" ~ 90,
                          meta.frame$drug[1] == "AMB" ~ 90)

# want 3 dfs for mic
for(i in 1:3){
  d.frame <- read_excel(mic_spreadsheet, 
                       range = meta.frame[[tolower(input_drug)]][i])
  
  # column-wise vs row-wise strain layout:
  # whichever (strain or concentration) is length 12 is assumed to be colnames
  colnames(d.frame) <- case_when(length(meta.frame$strain[!is.na(meta.frame$strain)]) == 12 ~ c(meta.frame$strain),
                                 length(meta.frame$concentration[!is.na(meta.frame$concentration)]) == 12 ~ c(as.character(meta.frame$concentration)))
  
  # add (strain or concentration) not used as column names as a new column
  if(!all(colnames(d.frame) == meta.frame$strain)) { d.frame$strain = meta.frame$strain[!is.na(meta.frame$strain)]}
  if(!all(colnames(d.frame) == meta.frame$concentration)) { d.frame$concentration = as.character(meta.frame$concentration[!is.na(meta.frame$concentration)])}
  
  # add remaining metadata 
  d.frame <- d.frame %>%
    add_column(plate = i) %>%
    add_column(drug = meta.frame$drug[!is.na(meta.frame$drug)]) %>%
    add_column(media = meta.frame$media[!is.na(meta.frame$media)]) %>%
    add_column(temp = meta.frame$temp[!is.na(meta.frame$temp)])

  assign(paste("rep", i,sep = ""), d.frame)
}
# pull together and tidy data
drug_data <- rbind(rep1,rep2,rep3)
drug <- pivot_longer(drug_data, 
                    names_to = meta_names,
                    values_to = "OD600", 
                    cols = -c(all_of(meta_col), plate, drug, media, temp)) 
drug <- drug %>%
    mutate(concentration = case_when((!("blank" %in%  c(strain, concentration)) & drug == "FLC" & concentration == "32" & strain %in% control_strains) ~ "blank",
                                     (!("blank" %in%  c(strain, concentration)) & drug == "MCF" & concentration == "1" & strain %in% control_strains) ~ "blank",
                                     (!("blank" %in%  c(strain, concentration)) & drug == "AMB" & concentration == "1" & strain %in% control_strains) ~ "blank",
                                     .default = concentration))
drug$concentration <-factor(tolower(drug$concentration), levels = mixedsort(unique(drug$concentration))) 
strains <- c(unique(drug$strain[drug$strain %in% control_strains]), rev(unique(drug$strain[!(drug$strain %in% control_strains)])))

################################################################################
# SMG
# Load metadata
smg.meta <- read_excel(smg_spreadsheet, 
                         sheet = 2)

smg.meta$drug <- c(toupper(input_drug), rep(NA, times=11))
# Make sure at least one column is length 12
#stopifnot( "Neither strain or concentration are right number of cols" = 12 %in% c(length(smg.meta$strain[!is.na(smg.meta$strain)]), length(smg.meta$concentration[!is.na(smg.meta$concentration)])))

# make sure sheets are matched
stopifnot("MIC and SMG strains don't match" = smg.meta$strain %in% meta.frame$strain)
stopifnot("drugs don't match" = smg.meta$drug[1]==meta.frame$drug[1])

smg.meta$concentration <- case_when(smg.meta$drug[1]== "FLC" ~ c(0,0.5,1,2,4,8,16,32,NA,NA,NA,NA),
                                      smg.meta$drug[1] %in% c("MCF","AMB") ~ c(0,0.016,0.032,0.064,0.125,0.256,0.5,1,NA,NA,NA,NA))

#set var for names when pivoting tidy data below
smg_names = case_when(length(smg.meta$strain[!is.na(smg.meta$strain)]) == 12 ~ "strain",
                       length(smg.meta$concentration[!is.na(smg.meta$concentration)]) == 12 ~ "concentration")
smg_col = case_when(length(smg.meta$strain[!is.na(smg.meta$strain)]) != 12 ~ "strain",
                     length(smg.meta$concentration[!is.na(smg.meta$concentration)]) != 12 ~ "concentration")

for(i in 1:3){
    d.frame <- read_excel(smg_spreadsheet, 
                          range = smg.meta[[tolower(input_drug)]][i])
    
    # column-wise vs row-wise strain layout:
    # whichever (strain or concentration) is length 12 is assumed to be colnames
    colnames(d.frame) <- case_when(length(smg.meta$strain[!is.na(smg.meta$strain)]) == 12 ~ c(smg.meta$strain),
                                   length(smg.meta$concentration[!is.na(smg.meta$concentration)]) == 12 ~ c(as.character(smg.meta$concentration)))
    
    # add (strain or concentration) not used as column names as a new column
    if(!all(colnames(d.frame) == smg.meta$strain)) { d.frame$strain = smg.meta$strain[!is.na(smg.meta$strain)]}
    if(!all(colnames(d.frame) == smg.meta$concentration)) { d.frame$concentration = as.character(smg.meta$concentration[!is.na(smg.meta$concentration)])}
    
    # add remaining metadata 
    d.frame <- d.frame %>%
        add_column(plate = i) %>%
        add_column(drug = smg.meta$drug[!is.na(smg.meta$drug)]) %>%
        add_column(media = smg.meta$media[!is.na(smg.meta$media)]) %>%
        add_column(temp = smg.meta$temp[!is.na(smg.meta$temp)])
    
    assign(paste("rep", i,sep = ""), d.frame)
}
# pull together and tidy data
drug48_data <- rbind(rep1,rep2,rep3)
drug48 <- pivot_longer(drug48_data, 
                     names_to = smg_names,
                     values_to = "OD600", 
                     cols = -c(all_of(smg_col), plate, drug, media, temp)) 
drug48 <- drug48 %>%
    mutate(concentration = case_when((!("blank" %in%  c(strain, concentration)) & drug == "FLC" & concentration == "32" & strain %in% control_strains) ~ "blank",
                                     (!("blank" %in%  c(strain, concentration)) & drug == "MCF" & concentration == "1" & strain %in% control_strains) ~ "blank",
                                     (!("blank" %in%  c(strain, concentration)) & drug == "AMB" & concentration == "1" & strain %in% control_strains) ~ "blank",
                                     .default = concentration))
drug48$concentration <-factor(tolower(drug48$concentration), levels = mixedsort(unique(drug48$concentration))) 

################################################################################
# calculate OD and plot
drug_od <- calculate_od(drug)

drug_mic <- mic(drug_od, cutoff)
drug_smg_input <- smg_input(drug_od, cutoff)

drug_smg_od <- calculate_od(drug48)
drug_smg <- smg_subset(drug_smg_od, drug_smg_input)

drug_plotting_coords <- plotting_coords(drug_od, drug_mic, cutoff, strains)
drug_mic_plot <- mic_plot(drug_od, strains, drug_plotting_coords) +
        xlab(case_when(drug$drug[1] %in% c("AMB", "MCF") ~ paste("\n", drug$drug[1]," MIC"),
             .default = paste(drug$drug[1], " MIC"))) +
        theme(axis.text.x = element_text(angle = x_axis_angle),
          plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
          axis.title.y = element_text(hjust = 0.5,
                                      margin = margin(0,10,0,0))) 
    
drug_smg_plot <- smg_plot(left_join(drug_mic, drug_smg), strains) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(),
         plot.margin = unit(c(0,0.5,0.5,0), "cm"),
          axis.title.x = element_text(lineheight = 1.2)) +
    scale_x_continuous(limits = c(0,1), breaks = c(0, 0.5, 1), labels = c("0", "", "1")) +
    xlab(case_when(drug$drug[1] %in% c("AMB", "MCF")~ paste("\n\n","SMG"),
         .default = "SMG")) 
    

drug_full_plot <- (drug_mic_plot + drug_smg_plot) + plot_layout(guides = 'collect')
ggsave(paste0("images/2023_MICs/",mic_date,"_MEC_",drug$drug[1],"_MIC.png"), 
       drug_full_plot, 
       width = 6, 
       height = 8, 
       units = "in", 
       device = png, 
       bg = "white",
       dpi = 300)
################################################################################
# work on imports
drug_mic$concentration <- as.character(drug_mic$concentration)

drug_mic <- drug_mic %>%
    mutate(concentration = case_when((mean_norm_OD > cutoff & concentration == max_concentration) ~ paste0(">", max_concentration),
                                     .default = concentration))

redcap_drug_value <- case_when(drug$drug[1] == "FLC" ~ 0,
                               drug$drug[1] == "MCF" ~ 1,
                               drug$drug[1] == "AMB" ~ 2)

redcap_drug_man <- case_when(drug$drug[1] == "FLC" ~ "Alfa Aesar",
                             drug$drug[1] == "MCF" ~ "MedChemExpress",
                             drug$drug[1] == "AMB" ~ "Chem-Impex Intl Inc")

redcap_drug_lot <- case_when(drug$drug[1] == "FLC" ~ "P04F029",
                             drug$drug[1] == "MCF" ~ "22367",
                             drug$drug[1] == "AMB" ~ "001448-180229")

redcap_drug_solvent <- case_when(drug$drug[1] == "FLC" ~ "EtOH",
                                 drug$drug[1] == "MCF" ~ "DMSO",
                                 drug$drug[1] == "AMB" ~ "DMSO")

for(i in 1:length(drug_mic$strain)){
    primary_id=drug_mic$strain[i]
    
  record <- c(
    primary_id = primary_id,
    redcap_repeat_instrument = "mic_results",
    redcap_repeat_instance = "new",
    drug = redcap_drug_value,
    mic_date = mic_date,
    mic50 = as.character(drug_mic$concentration[drug_mic$strain==primary_id]),
    mic_file_name = basename(mic_spreadsheet),
    smg = as.character(drug_smg$mean_48[drug_smg$strain==primary_id]),
    smg_file_name = basename(smg_spreadsheet),
    mic_media = 1,
    mic_temp = 1,
    drug_source = redcap_drug_man,
    lot_number = redcap_drug_lot,
    solvent = redcap_drug_solvent,
    mic_results_complete = 2
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
