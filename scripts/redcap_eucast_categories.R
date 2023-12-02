## ---------------------------
## Script name: redcap_eucast_calcs_api.R
##
## Purpose of script: Pull report CSVs and determine eucast categories when applicable
## Author: Nancy Scott
##
## Date Created: 2023-09-06
##
## Email: scot0854@umn.edu
## ---------------------------
## Notes:
## ---------------------------
## load packages
library(tidyverse)
library(lubridate)
#library(magrittr)
library(writexl)
#library(reshape2)
library(RCurl)
library(jsonlite)
# redcap report IDs
samples <- '58043'
mic_results <- '58044'

token <- "" # don't forget to delete before gh
api_url <-  "https://redcap.ahc.umn.edu/api/"

# add max concentration from MIC plates (below is standard MEC screening set-up)
# change if needed
max_flc <- 32
max_amb <- 1.00
max_mcf <- 1.00

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

# Sample report 
sample_info <- import_report(samples) %>%
    select(-c(starts_with('redcap_repeat')))

# make species df for use with breakpoint info
species <- sample_info %>%
    select(primary_id, genus_species)

non_eucast_species <- c("Candida lusitaniae", "Candida kefyr",
                        "Candida nivariensis", "Candida orthopsilosis", "Candida utilis")

# get MIC data and add species
eucast_mic <- import_report(mic_results) %>%
    filter(redcap_repeat_instrument != "NA") %>%
    select(primary_id, redcap_repeat_instance, drug, mic_date, mic50, mic_media, mic_temp, smg) %>%
    left_join(species, by = join_by("primary_id"))

eucast_mic$mic50 <- case_when(grepl(">", eucast_mic$mic50) ~ (as.numeric(sub("[^-.0-9]", "", eucast_mic$mic50)) + 0.001),
                              .default = (as.numeric(eucast_mic$mic50)))

eucast_mic <- eucast_mic %>%
    filter(mic_media == "RPMI" & mic_temp=="35C") %>%
    filter(mic_date %in% c(as.Date("2023-11-21"), as.Date("2023-11-22"))) %>%
    group_by(genus_species) %>%
    mutate(eucast_breakpoint = case_when(drug!= "fluconazole" & genus_species %in% non_eucast_species ~ 3,
                                         drug == "fluconazole" & genus_species == "Candida krusei" ~ 4,
                                         drug == "fluconazole" & genus_species == "Candida glabrata" & mic50 <= 16 ~ 1,
                                         drug == "fluconazole" & genus_species == "Candida glabrata" & mic50 > 16 ~ 2,
                                         drug == "fluconazole" & !(genus_species %in% c("Candida glabrata", "Candida krusei")) & mic50 <= 2 ~ 0,
                                         drug == "fluconazole" & !(genus_species %in% c("Candida glabrata", "Candida krusei"))  & between(mic50, 2.001, 4) ~ 1,
                                         drug == "fluconazole" & !(genus_species %in% c("Candida glabrata", "Candida krusei"))  & mic50 > 4 ~ 2,
                                         drug == "micafungin" & genus_species == "Candida albicans" & mic50 <= 0.016 ~ 0,
                                         drug == "micafungin" & genus_species == "Candida albicans" & mic50 > 0.016 ~ 2,
                                         drug == "micafungin" & genus_species == "Candida glabrata" & mic50 <= 0.032 ~ 0,
                                         drug == "micafungin" & genus_species == "Candida glabrata" & mic50 > 0.032 ~ 2,
                                         drug == "micafungin" & genus_species == "Candida parapsilosis"  & mic50 <= max_mcf & max_mcf < 2 ~ 0,
                                         drug == "micafungin" & genus_species == "Candida parapsilosis" & mic50 > 2 ~ 2,
                                         drug == "micafungin" & genus_species %in% c("Candida krusei", "Candida tropicalis") ~ 3,
                                         drug == "amphotericin B" & genus_species %in% c("Candida albicans",
                                                                                         "Candida dubliniensis",
                                                                                         "Candida glabrata",
                                                                                         "Candida krusei",
                                                                                         "Candida parapsilosis",
                                                                                         "Candida tropicalis") & mic50 <= 1 ~ 0,
                                         drug == "amphotericin B" & genus_species %in% c("Candida albicans",
                                                                                         "Candida dubliniensis",
                                                                                         "Candida glabrata",
                                                                                         "Candida krusei",
                                                                                         "Candida parapsilosis",
                                                                                         "Candida tropicalis") & mic50 > 1 ~ 2))


# Upload EUCAST categories to REDCap

for(i in 1:length(eucast_mic$primary_id)){
    
    record <- c(
        primary_id = eucast_mic$primary_id[i],
        redcap_repeat_instrument = "mic_results",
        redcap_repeat_instance = eucast_mic$redcap_repeat_instance[i],
        eucast_breakpoint = eucast_mic$eucast_breakpoint[i]
    )
    
    result_data <- toJSON(list(as.list(record)), auto_unbox=TRUE)
    
    formData <- list("token"=token,
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
