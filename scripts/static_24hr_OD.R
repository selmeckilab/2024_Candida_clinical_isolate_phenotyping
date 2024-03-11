## ---------------------------
## Purpose: Calculate blank-corrected mean ODs for 24-hour static growth in RPMI
## Author: Nancy Scott
## Email: scot0854@umn.edu
## ---------------------------
options(scipen = 999)
## ---------------------------
read_meta <- function(spreadsheet,tab){
    read_excel(spreadsheet, 
               sheet = tab,
               range = "A1:F13")
}
source("scripts/MIC_heatmap.R")

# Enter drug and spreadsheet data
input_drug <- c("flc","mcf", "amb")
replicates <- 3
od_tab <- 1  # Excel tab of OD values
metadata_tab <- 2 # Excel tab of metadata
mic_spreadsheet <-"data/MIC/2023-12-19_EW_MIC24_RPMI35.xlsx"

# Type either "strain" or "concentration" for col names, depending on plate layout
column_names <- "strain"

concentration <- c(0, 0.016, 0.032, 0.064, 0.125, 0.256, 0.5, 1)

control_strains <- c("AMS5123", "AMS5122")

mic_date <-str_extract(mic_spreadsheet, "\\d+-\\d+-\\d+")
################################################################################
# Load metadata
meta.frame <- read_meta(mic_spreadsheet, metadata_tab)

meta.frame$concentration <- c(concentration, rep(NA, times=12 - length(concentration)))

# Set vars for use below
meta_names = case_when(column_names=="strain" ~ "strain",
                       column_names=="concentration" ~ "concentration")
meta_col = case_when(column_names!="strain" ~ "strain",
                     column_names !="concentration" ~ "concentration")
    
################################################################################
# Loop over each plate
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
# Pull together and tidy data
drug_data <- bind_rows(mget(ls(pattern = "^rep\\w+\\d+$")))
drug <- pivot_longer(drug_data, 
                    names_to = meta_names,
                    values_to = "OD600", 
                    cols = -c(all_of(meta_col), plate, drug, media, temp)) 

# BLANKS: this checks for "blank" labeled wells in metadata
# if none found, sets it to the highest concentration well of the control strain
drug <- drug %>%
    mutate(concentration = case_when((!("blank" %in%  c(strain, concentration)) & drug == "FLC" & concentration == "1" & strain %in% control_strains) ~ "blank",
                                     (!("blank" %in%  c(strain, concentration)) & drug == "MCF" & concentration == "1" & strain %in% control_strains) ~ "blank",
                                     (!("blank" %in%  c(strain, concentration)) & drug == "AMB" & concentration == "1" & strain %in% control_strains) ~ "blank",
                                     .default = concentration))
drug$concentration <-factor(tolower(drug$concentration), levels = mixedsort(unique(drug$concentration))) 

# consider min bg instead of mean
# re-do for 2024 vals or just skip?
background <- drug %>% 
    filter(concentration=="blank", OD600 <0.17) %>% 
    summarize(mean_bg = mean(OD600))

no_drug <- drug %>% 
    filter(concentration == 0) %>% 
    mutate(OD600 = OD600 - background$mean_bg) %>% 
    group_by(strain) %>% 
    summarise( mean_OD=mean(OD600), sd_od = sd(OD600))

mic_boxplot(drug)

write_csv(no_drug, paste0(mic_date, "_no_drug_norm_OD.csv"))
