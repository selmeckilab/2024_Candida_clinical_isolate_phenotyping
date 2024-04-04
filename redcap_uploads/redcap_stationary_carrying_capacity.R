## ---------------------------
## Purpose: Calculate blank-corrected mean ODs for 24-hour static growth in RPMI
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

# API variables
api_token <- ""
api_url <-  "https://redcap.ahc.umn.edu/redcap/api/"

# Spreadsheet variables
input_drug <- c("flc","mcf", "amb")
replicates <- 3
od_tab <- 1  # Excel tab of OD values
metadata_tab <- 2 # Excel tab of metadata
mic_spreadsheet <-"data/MIC/2023-12-05_EW_MIC24_RPMI35.xlsx"

# Type either "strain" or "concentration" for col names, depending on plate layout
column_names <- "strain"

# Placeholder concentrations, only care about 0
concentration <- c(0, 0.016, 0.032, 0.064, 0.125, 0.256, 0.5, 1)

control_strains <- c("AMS5123", "AMS5122")

# Libraries
library(readxl)
library(jsonlite)

# Functions
read_meta <- function(spreadsheet,tab){
    read_excel(spreadsheet, 
               sheet = tab,
               range = "A1:F13")
}
source("scripts/MIC_heatmap.R")

# Load metadata

mic_date <-str_extract(mic_spreadsheet, "\\d+-\\d+-\\d+")

meta.frame <- read_meta(mic_spreadsheet, metadata_tab)

meta.frame$concentration <- c(concentration, rep(NA, times=12 - length(concentration)))

# Set vars for use below
meta_names = case_when(column_names=="strain" ~ "strain",
                       column_names=="concentration" ~ "concentration")
meta_col = case_when(column_names!="strain" ~ "strain",
                     column_names !="concentration" ~ "concentration")
    
# Loop over all plates and concatenate
for(j in 1:3){
  j_drug <- input_drug[j]
  for(i in 1:replicates){
    d.frame <- read_excel(mic_spreadsheet, 
                        sheet = od_tab,
                       range = meta.frame[[tolower(j_drug)]][i])
  
  # Colnames set from input variable 
    colnames(d.frame) <- as.character(pull(meta.frame,column_names))
  
  # Add (strain or concentration) not used as column names as a new column
    if(!all(meta.frame$strain[!is.na(meta.frame$strain)] %in% colnames(d.frame))) { d.frame$strain = meta.frame$strain[!is.na(meta.frame$strain)]}
    if(!all(meta.frame$concentration[!is.na(meta.frame$concentration)] %in% colnames(d.frame))) { d.frame$concentration = as.character(meta.frame$concentration[!is.na(meta.frame$concentration)])}
  
  # Add remaining metadata 
    d.frame <- d.frame %>%
      add_column(plate = i) %>%
      add_column(drug = toupper(j_drug)) %>%
      add_column(media = meta.frame$media[!is.na(meta.frame$media)]) %>%
      add_column(temp = meta.frame$temp[!is.na(meta.frame$temp)])

    assign(paste0("rep", j_drug, i), d.frame)
  }
}

drug_data <- bind_rows(mget(ls(pattern = "^rep\\w+\\d+$")))

# Tidy data
drug <- pivot_longer(drug_data, 
                    names_to = meta_names,
                    values_to = "OD530", 
                    cols = -c(all_of(meta_col), plate, drug, media, temp)) 

# BLANKS: this checks for "blank" labeled wells in metadata
# if none found, sets it to the highest concentration well of the control strain
drug <- drug %>%
    mutate(concentration = case_when((!("blank" %in%  c(strain, concentration)) & drug == "FLC" & concentration == "1" & strain %in% control_strains) ~ "blank",
                                     (!("blank" %in%  c(strain, concentration)) & drug == "MCF" & concentration == "1" & strain %in% control_strains) ~ "blank",
                                     (!("blank" %in%  c(strain, concentration)) & drug == "AMB" & concentration == "1" & strain %in% control_strains) ~ "blank",
                                     .default = concentration))
drug$concentration <-factor(tolower(drug$concentration), levels = mixedsort(unique(drug$concentration))) 

mic_boxplot(drug)

# Is there a lot of variation in RPMI compared to YPAD?
# Removing some of the higher vals before averaging
background <- drug %>% 
    filter(tolower(concentration)=="blank" | tolower(strain)=="blank", OD530 <0.17) %>% 
    summarize(mean_bg = mean(OD530))

background$mean_bg

# Get corrected values for upload
no_drug <- drug %>% 
    filter(concentration == 0, tolower(strain)!="blank") %>% 
    mutate(OD530 = OD530 - background$mean_bg) %>% 
    group_by(strain) %>% 
    summarise( mean_OD=round(mean(OD530), digits=3), sd_OD = round(sd(OD530), digits = 3))

# For each strain, create new record and send form 
rpmi_temp <- case_when(meta.frame$temp[1]==30 ~ 0,
                       meta.frame$temp[1]==35 ~ 1,
                       meta.frame$temp[1]==37 ~2)

for(i in 1:length(no_drug$strain)){
    primary_id=no_drug$strain[i]
    
    record <- c(
        primary_id = primary_id,
        redcap_repeat_instrument = "rpmi_growth_data",
        redcap_repeat_instance = "new",
        rpmi_date = mic_date,
        rpmi_temp = rpmi_temp,
        mean_static_od_24h = as.character(no_drug$mean_OD[no_drug$strain==primary_id]),
        sd_static_od_24h = as.character(no_drug$sd_OD[no_drug$strain==primary_id]),
        rpmi_growth_data_complete = 2
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
