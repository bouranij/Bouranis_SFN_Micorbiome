# Environment -------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(lubridate)
library(cowplot)
library(PKNCA)
setwd('~/Documents/Manuscripts/BSS_Nitriles/Bouranis_SFN_Microbiome/')
here::i_am('./PK_Analysis_Script.R')
library(here)

#Set Up Theme:
theme_sfn <- function(x){
  cowplot::theme_cowplot() %+replace%
  theme(
    update_geom_defaults('point', list(shape = 16, size = 4)),
    update_geom_defaults('path', list(size = 2))) 
}

#ggplot shortcuts because I'm lazy
drop_legend <- theme(legend.position = 'none')
drop_gap <- scale_y_continuous(expand = c(0,0))

#Set Cohort Colors:
chrt_col <- c('#E6122F', '#FAD20A', '#ED6700', '#6140B5', '#0A64FF', '#04B44A', '#49A9BA', '#C954A6')
chrt_1 <- c('#740615', '#9E0519', '#E6122F', '#B83E4F')
chrt_2 <- c('#FFEB85', '#FAD20A','#B59F31', '#806D12')
chrt_3 <- c('#ED6700', '#EB7E2B', '#B55810', '#F6A465')
chrt_4 <- c('#4F3494', '#54496E', '#705F9C', '#6140B5', '#8263CF', '#BAA2F6')
chrt_5 <- c('#0A64FF', '#003CA3', '#81A9EE')
chrt_6 <- c('#0C6F33', '#188C46', '#04B44A', '#1BDA67', '#86E9AE', '#96DFA3')
chrt_7 <- c('#15444D', '#0D7082', '#52939E', '#7DBECA', '#49A9BA', '#52D4EB', '#91E4F2', '#C8F2F9')
chrt_8 <- c('#9F1977', '#AB498E', '#C954A6')

BSS_col <- c(chrt_1, chrt_2, chrt_3, chrt_4, chrt_5, chrt_6, chrt_7, chrt_8)
BSS_col <- c(chrt_3, chrt_4, chrt_5, chrt_6, chrt_7, chrt_8)

pie(rep(1, 30), col = BSS_col)

#Set Broccoli and Alfalfa Colors:
BA_col <- c('#3AE151', '#DDE712')

#Set Label and Unlabeled Colors:
LUL_col <- c('#3A43E0', '#E81368')

#Useful Functions:
se <- function(x){
  sd(x, na.rm = T)/sqrt(length(x))
}

#Helpful Functions
log_helper <- function(x, min.val) {
  log2((x + sqrt(x ^ 2 + min.val ^ 2)) / 2)
}

#Pareto Scaling:
PS_helper <- function(x) {
  (x - mean(x)) / sqrt(sd(x, na.rm = T))
}	

#Auto Scaling:
AS_helper <- function(x) {
  (x - mean(x)) / sd(x, na.rm = T)
} 

auto_scale <- function(mtb){
  mtb_scaled <- apply(mtb, 2, AS_helper) 
  return(mtb_scaled)
}

#Transformation Functions:
#Log Scaling:
log_transform <- function(mtb){
  mtb_nz <- mtb[ ,which(apply(mtb, 2, sum) != 0)]
  min.val <- min(abs(mtb_nz[mtb_nz!=0]))/10
  mtb_log_trans <- apply(mtb_nz, 2, log_helper, min.val)
  return(mtb_log_trans)
}


parse_diet_logs <- function(path){
  #Grab the subject ID
  ID <- suppressMessages(read_tsv(path, col_names = F)[[1,1]]) %>%
   gsub('Spreadsheet: ', '', .) %>%
   gsub('[ -]', '', .)
  #Read in the log
  testlog <- suppressMessages(read_tsv(path, skip = 2, na = '--'))
  #Pull out just the summary info
  summary_info <- testlog[str_detect(testlog$`Item Name`, 'Day'), ]
  summary_clean <- summary_info[,-c(2,3)] %>%
   mutate(rdsplit = str_split(string = `Item Name`, pattern = '\\(')) %>%
   mutate(relative_day = map_chr(rdsplit, function(x) gsub(pattern = ' ', replacement = '', x = x[1], ))) %>%
   dplyr::select(-`Item Name`, -rdsplit, -`Wgt (g)`) %>%
   pivot_longer(ends_with(')'), names_to = 'nutrient') %>%
   cbind(ID, .)
  return(summary_clean)
}

read_in_logs <- function(folder){
  map_df(folder, parse_diet_logs)
}

# Data Load In ------------------------------------------------------------

## SFN Data + Cleaning -----------------------------------------------------

#Metadata:
metadata_urine <- read_csv(here('Data/Metadata/Meta_urine.csv'))
metadata_sprout <- read_csv(here('Data/Metadata/sprout_meta.csv')) %>%
  modify_at(c('cohort'), as.factor)
metadata_treatment <- read_csv(here('Data/Metadata/meta_treatment.csv'))
metadata_fecal <- read_csv(here('Data/Metadata/meta_fecal_broconly.csv'))

#Raw Data:
urine_raw <- read_csv(here('Data/SFN_Data/Urine_R_Format.csv'))
fecal_raw <- read_csv(here('Data/SFN_Data/Fecal_R_format_update.csv'))
plasma_raw <- read_csv(here('Data/SFN_Data/Plasma_R_Format.csv'))

#Set Subject ID Color Names
subid <- metadata_treatment %>%
  filter(treatment %in% c('BU', 'BL')) %>%
  filter(cohort %in% 3:8) %>%
  arrange(cohort) %>%
  pull(subject_id)
names(BSS_col) <- subid
pal_bss <- scale_color_manual(values = BSS_col)

#Set Cohort Color Names
names(chrt_col) <- 1:8
pal_chrt <- scale_color_manual(values = chrt_col)

#Set Broccoli Alfalfa Names:
names(BA_col) <- c('Broccoli', 'Alfalfa')
pal_ba <- scale_color_manual(values = BA_col)

#Set Label/Unlabel Names:
names(LUL_col) <- c('BL', 'BU')
pal_lb <- scale_color_manual(values = LUL_col)


#Urine:
#Clean the metadata
mdata_urine_clean <- metadata_urine %>%
  #Weight (grams) used as a  proxy for mL, convert to L
  mutate(across(c(4:9), ~.x/1000)) %>%
  #Pivot the data longer
  pivot_longer(cols = starts_with('vol'), names_to = 'time', values_to = 'volume') %>%
  #Rename the variables so they can treated as factors
  mutate(time = gsub('vol_', '', time)) %>%
  mutate(time = gsub('h','', time)) 

#Clean the urine data
urine_clean <- urine_raw %>%
  #Rename the time points so they match the metadata
  mutate(time = gsub('h','', time)) %>%
  #Join the data with the metadata
  left_join(., mdata_urine_clean) %>%
  group_by(subject_id, time) %>%
  #Multiply all the uM by the (proxy) volumes to convert to uMol recovered
  mutate(across(starts_with('SFN'), ~.x*volume)) %>%
  #Convert to factors for analysis + graphing
  modify_at(c('subject_id', 'time', 'cohort'), as.factor)
#Relevel the time variables to make it work
urine_clean$time %<>% fct_relevel(c('0', '3', '6', '24', '48', '72'))

#Drop cohorts 1 and 2
urine_clean_sub <- urine_clean %>%
  filter(!cohort %in% c(1,2))

#Plasma:
plasma_clean <- plasma_raw %>%
  mutate(time = gsub('h','', time)) %>%
  #Join the data with the metadata
  left_join(., metadata_treatment) %>%
  #Convert to factors for analysis + graphing
  modify_at(c('subject_id', 'time', 'cohort'), as.factor) %>%
  modify_at(4:10, as.numeric) 
plasma_clean$time %<>% fct_relevel(c('0', '3', '6', '24', '48', '72'))

#Drop cohorts 1 and 2
plasma_clean_sub <- plasma_clean %>%
  filter(!cohort %in% c(1,2))

#Fecal:
fecal_clean <- fecal_raw %>%
  rename_at(c('Treatment', 'Time'), tolower) %>%
  mutate(time = gsub('h','', time)) %>%
  #Join the data with the metadata
  left_join(., metadata_treatment) %>%
  #Convert to factors for analysis + graphing
  modify_at(c('subject_id', 'time', 'cohort'), as.factor) %>%
  modify_at(5:10, as.numeric) 
fecal_clean$time %<>% fct_relevel(c('0','24', '48', '72'))

#Drop cohorts 1 and 2
fecal_clean_sub <- fecal_clean %>%
  filter(!cohort %in% c(1,2))

#Make the alfalfa data only
urine_clean_alf <- filter(urine_clean_sub, treatment %in% c('AL', 'AU'))
plasma_clean_alf <- filter(plasma_clean_sub, treatment %in% c('AL', 'AU'))
fecal_clean_alf <- filter(fecal_clean_sub, treatment %in% c('AL', 'AU'))

#Filter our data to remove alfalfa folks
urine_clean_sub %<>% filter(treatment %in% c('BL', 'BU'))
plasma_clean_sub %<>% filter(treatment %in% c('BL', 'BU'))
fecal_clean_sub %<>% filter(treatment %in% c('BL', 'BU'))

aurine <- urine_clean_alf %>%
  pivot_longer(cols = c(starts_with('SFN'),-SFN_Tot), names_to = 'metabolite', values_to = 'uMol') %>%
  group_by(time, metabolite) %>%
  summarise(mean = mean(uMol, na.rm = T)) %>%
  pivot_wider(names_from = 'time', values_from = 'mean')

aplasma <- plasma_clean_alf %>%
  pivot_longer(cols = c(starts_with('SFN'),-SFN_Tot), names_to = 'metabolite', values_to = 'uMol') %>%
  group_by(time, metabolite) %>%
  summarise(mean = mean(uMol, na.rm = T)) %>%
  pivot_wider(names_from = 'time', values_from = 'mean')

afecal <- fecal_clean_alf %>%
  pivot_longer(cols = c(starts_with('SFN'),-SFN_Tot), names_to = 'metabolite', values_to = 'uMol') %>%
  group_by(time, metabolite) %>%
  summarise(mean = mean(uMol, narm = T)) %>%
  pivot_wider(names_from = 'time', values_from = 'mean')

#Make data tidy for plotting and analysis
urine_tidy <- urine_clean_sub %>%
  pivot_longer(cols = c(starts_with('SFN'),-SFN_Tot), names_to = 'metabolite', values_to = 'uMol')

#Make data tidy for plotting and analysis
plasma_tidy <- plasma_clean_sub %>%
  pivot_longer(cols = c(starts_with('SFN'),-SFN_Tot), names_to = 'metabolite', values_to = 'uMol')

#Make data tidy for plotting and analysis
fecal_tidy <- fecal_clean_sub %>%
  pivot_longer(cols = c(starts_with('SFN'),-SFN_Tot), names_to = 'metabolite', values_to = 'uMol')


# PK Analysis -------------------------------------------------------------
library(PKNCA)

## Plasma Data --------------------------------------------------------------
pkdata <- plasma_clean_sub %>%
  mutate(DTC_Tot = SFN + SFN_Cys + SFN_NAC + SFN_GSH + SFN_CG) %>%
  #dplyr::select(subject_id, time, DTC_Tot, SFN_NIT, SFN_Tot, treatment) %>%
  magrittr::inset('dose', value = 100) %>%
  modify_if(is.factor, as.character) %>%
  modify_at('time', as.integer)

dose_data <- pkdata %>%
  filter(time == 0) %>%
  dplyr::select(subject_id, time, dose, treatment)

conc_obj_DTC <- PKNCAconc(pkdata, DTC_Tot~time|subject_id)
conc_obj_NIT <- PKNCAconc(pkdata, SFN_NIT~time|subject_id)
conc_obj_SFN <- PKNCAconc(pkdata, SFN~time|subject_id)
conc_obj_SFN_Cys <- PKNCAconc(pkdata, SFN_Cys~time|subject_id)
conc_obj_SFN_CG <- PKNCAconc(pkdata, SFN_CG~time|subject_id)
conc_obj_SFN_GSH <- PKNCAconc(pkdata, SFN_GSH~time|subject_id)
conc_obj_SFN_NAC <- PKNCAconc(pkdata, SFN_NAC~time|subject_id)



conc_obj_DTC_label <- PKNCAconc(pkdata, DTC_Tot~time|subject_id/treatment)
conc_obj_NIT_label <- PKNCAconc(pkdata, SFN_NIT~time|subject_id/treatment)
conc_obj_SFN_label <- PKNCAconc(pkdata, SFN_Tot~time|subject_id/treatment)

dose_obj <- PKNCAdose(dose_data, dose~time|subject_id)
dose_obj_label <- PKNCAdose(dose_data, dose~time|subject_id)

intervals_manual <- data.frame(start = 0, end = c(72, Inf),
                               cmax = c(TRUE, TRUE),
                               tmax = c(TRUE, TRUE),
                               auclast = c(TRUE, TRUE),
                               aucinf.obs = c(TRUE, TRUE),
                               half.life = c(TRUE, TRUE),
                               aucinf.pred = c(TRUE, TRUE),
                               aucall = c(TRUE,TRUE))

auto_DTC <- PKNCAdata(conc_obj_DTC, dose_obj, intervals = intervals_manual)
auto_NIT <- PKNCAdata(conc_obj_NIT, dose_obj, intervals = intervals_manual)
auto_SFN <- PKNCAdata(conc_obj_SFN, dose_obj, intervals = intervals_manual)
auto_SFN_Cys <- PKNCAdata(conc_obj_SFN_Cys, dose_obj, intervals = intervals_manual)
auto_SFN_CG <- PKNCAdata(conc_obj_SFN_CG, dose_obj, intervals = intervals_manual)
auto_SFN_GSH <- PKNCAdata(conc_obj_SFN_GSH, dose_obj, intervals = intervals_manual)
auto_SFN_NAC <- PKNCAdata(conc_obj_SFN_NAC, dose_obj, intervals = intervals_manual)

auto_DTC_label <- PKNCAdata(conc_obj_DTC_label, dose_obj, intervals = intervals_manual)
auto_NIT_label <- PKNCAdata(conc_obj_NIT_label, dose_obj, intervals = intervals_manual)
auto_SFN_label <- PKNCAdata(conc_obj_SFN_label, dose_obj, intervals = intervals_manual)

pknca_units_table(concu = 'µM', doseu = 'µMol SFN Equivalents', timeu = 'hr')
PKNCA.options(min.hl.points = 2, allow.tmax.in.half.life = T)


DTC_final <- pk.nca(auto_DTC)
NIT_final <- pk.nca(auto_NIT)
SFN_final <- pk.nca(auto_SFN)
SFN_Cys_final <- pk.nca(auto_SFN_Cys)
SFN_CG_final <- pk.nca(auto_SFN_CG)
SFN_NAC_final <- pk.nca(auto_SFN_NAC)
SFN_GSH_final <- pk.nca(auto_SFN_GSH)

DTC_final_label <- pk.nca(auto_DTC_label)
NIT_final_label <- pk.nca(auto_NIT_label)
SFN_final_label <- pk.nca(auto_SFN_label)

mean_na <- function(x){
  mean(x, na.rm = T)
}

sd_na <- function(x){
  sd(x, na.rm = T)
}


PKNCA.set.summary(name = c('auclast', 'cmax', 'tmax', 'half.life', 'aucinf.obs', 'aucinf.pred', 'aucall'),
                  point = mean_na,
                  spread = sd_na,
                  description = 'Means and SDs')


summary(DTC_final)
summary(NIT_final)
summary(SFN_final)
summary(SFN_GSH_final)
summary(SFN_CG_final)
summary(SFN_Cys_final)
summary(SFN_NAC_final)

summary(DTC_final_label)
summary(NIT_final_label)
summary(SFN_final_label)

#Calculate total excreted mmols in urine
urine_clean_sub %>%
  mutate(DTC_Tot = SFN + SFN_Cys + SFN_NAC + SFN_GSH + SFN_CG) %>%
  pivot_longer(cols = c('SFN', 'SFN_Cys', 'SFN_NAC', 'SFN_CG', 'SFN_GSH', 'DTC_Tot', 'SFN_NIT'), 
               names_to = 'metab', values_to = 'value') %>%
  group_by(subject_id, metab) %>%
  summarise(total = sum(value)) %>%
  ungroup() %>%
  group_by(metab) %>%
  summarise(mean = mean_na(total),
            sd = sd_na(total)) %>%
  mutate(across(2:3, ~.x/100, .names = "{.col}_frac"))







