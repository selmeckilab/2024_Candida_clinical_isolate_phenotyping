## ---------------------------
## Purpose: Carrying capacity (mean 24-hour OD, stationary growth in RPMI) 
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)

# Libraries
library(readxl)
library(jsonlite)

api_url <-  "https://redcap.ahc.umn.edu/redcap/api/"
mic_files <- '64707'

# Functions
import_report <- function(report_number) {
    url <- api_url
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

read_meta <- function(spreadsheet,tab){
    read_excel(spreadsheet, 
               sheet = tab,
               range = "A1:F13")
}

files_todo <- import_report(mic_files) %>% 
    filter(redcap_repeat_instrument != "NA",mic_media=="RPMI liquid",
           !mic_date %in% c(as.Date("2023-11-16"), as.Date("2024-02-20"), as.Date("2024-02-27"), as.Date("2024-03-07")) )

mic_files <- sort(unique(files_todo$mic_file_name))

smg_files <- sort(unique(files_todo$smg_file_name))

# Spreadsheet variables
input_drug <- c("flc","mcf", "amb")
replicates <- 3
od_tab <- 1  # Excel tab of OD values
metadata_tab <- 2 # Excel tab of metadata

# Iterate over list of files
for(k in 3:length(smg_files)){
  mic_spreadsheet <- paste0("data/MIC/", smg_files[k])

  # Type either "strain" or "concentration" for col names, depending on plate layout
  column_names <- "strain"

  # Placeholder concentrations, only care about 0
  concentration <- c(0, 0.016, 0.032, 0.064, 0.125, 0.256, 0.5, 1)

  control_strains <- c("AMS5123", "AMS5122")
 
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
              add_column(temp = meta.frame$temp[!is.na(meta.frame$temp)]) %>% 
              add_column(smg_file_name = smg_files[k])
          
          assign(paste0("rep", j_drug, i), d.frame)
      }
  }
  
  drug_data <- bind_rows(mget(ls(pattern = "^rep\\w+\\d+$")))
  
  # Tidy data
  drug <- pivot_longer(drug_data, 
                       names_to = meta_names,
                       values_to = "OD600", 
                       cols = -c(all_of(meta_col), plate, drug, media, temp, smg_file_name)) 
  
  # BLANKS: this checks for "blank" labeled wells in metadata
  # if none found, sets it to the highest concentration well of the control strain
  drug <- drug %>%
      mutate(concentration = case_when((!("blank" %in%  c(strain, concentration)) & drug == "FLC" & concentration == "1" & strain %in% control_strains) ~ "blank",
                                       (!("blank" %in%  c(strain, concentration)) & drug == "MCF" & concentration == "1" & strain %in% control_strains) ~ "blank",
                                       (!("blank" %in%  c(strain, concentration)) & drug == "AMB" & concentration == "1" & strain %in% control_strains) ~ "blank",
                                       .default = concentration))
  drug$concentration <-factor(tolower(drug$concentration), levels = mixedsort(unique(drug$concentration))) 
  
  # Is there a lot of variation in RPMI compared to YPAD?
  # Removing some of the higher vals before averaging
  background <- drug %>% 
      group_by(drug) %>% 
      filter(tolower(concentration)=="blank" | tolower(strain)=="blank", OD600 <0.17) %>% 
      summarize(mean_bg = mean(OD600))
  
  print(background$mean_bg)
  
  # Get corrected values for upload
  control_wells <- drug %>% 
      filter(concentration == 0, tolower(strain)!="blank") %>% 
      mutate(OD600 = OD600 - background$mean_bg) %>% 
      group_by(strain, drug, smg_file_name) %>% 
      summarise( mean_OD=round(mean(OD600), digits=3), sd_OD = round(sd(OD600), digits = 3))
  
  control_wells <- control_wells %>% 
      mutate(drug = replace(drug, drug=="AMB", "amphotericin B")) %>% 
      mutate(drug = replace(drug, drug=="FLC", "fluconazole")) %>% 
      mutate(drug= replace(drug, drug=="MCF", "micafungin")) %>% 
      inner_join(files_todo, by = c("strain" = "primary_id",
                                    "drug",
                                    "smg_file_name"))
  
  
  for(l in 1:length(control_wells$strain)){
      record <- c(
          primary_id = control_wells$strain[l],
          redcap_repeat_instrument = "mic_results",
          redcap_repeat_instance = as.character(control_wells$redcap_repeat_instance[l]),
          mean_stationary_k_48 = as.character(control_wells$mean_OD[l]),
          sd_stationary_k_48 = as.character(control_wells$sd_OD[l])
      )
      result_data <- toJSON(list(as.list(record)), auto_unbox=TRUE)
      formData <- list("token"=Sys.getenv("redcap_api_key"),
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
}
