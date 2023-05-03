# Environment -------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(lubridate)
library(plotly)
library(cowplot)
library(phyloseq)
setwd('~/Documents/Manuscripts/BSS_Nitriles/Bouranis_SFN_Microbiome/')
here::i_am('./Figure_Script.R')
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
metadata_health <- read_csv(here('Data/Metadata/meta_health.csv'))
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

## Microbiome Data ---------------------------------------------------------
ps_raw <- readRDS('./Data/Plot_Data/ps_raw.RDS')

#Taxa Filtering:
#ps_raw has 4060 ASVs
#Remove taxa that do not appear at least 3 times in at least 20% of all samples
ps_counts <- ps_raw %>% filter_taxa(function(x) sum(x > 3) > (0.2*length(x)), TRUE) #317
#Convert from counts to relative abundance
ps_relab <- ps_counts %>% transform_sample_counts(function(x) x / sum(x) ) #317
#Filter out low abundance (>1e-5) taxa
ps <- ps_relab %>% filter_taxa(function(x) mean(x) > 1e-5, TRUE) #317

ps #317

## Cell Culture Data -------------------------------------------------------

se <- function(x){
  sd(x)/sqrt(length(x))
}

hdac <- read_csv(here('./Data/Plot_Data/HDAC_Acitivty.csv')) %>%
  modify_at('treatment', ~gsub(pattern = '-[0-9]', replacement = '', x = .x)) %>%
  modify_at('treatment', factor, levels = c('Vehicle ctl (0.27% DMSO) ', 'SFN 15uM ', 'SFN-nitrile 15uM ', 'SFN-nitrile 45uM ', 'SFN-nitrile 135uM ')) 

hdac_summary <- hdac %>%
  group_by(treatment) %>%
  summarise(mean = mean(FC), 
            ser = se(FC)) 

qpcr_24 <- read_csv(here('./Data/Plot_Data/24h_qpcr.csv')) %>%
  modify_at('treatment', ~gsub(pattern = '-[0-9]', replacement = '', x = .x)) %>%
  modify_at('treatment', factor, levels = c('Vehicle ctl (0.27% DMSO) ', 'SFN 15uM ', 'SFN-nitrile 15uM ', 'SFN-nitrile 45uM ', 'SFN-nitrile 135uM ')) %>%
  inset('timepoint', value = '24h')

qpcr_48 <- read_csv(here('./Data/Plot_Data/48h_qpcr.csv')) %>%
  modify_at('treatment', ~gsub(pattern = '-[0-9]', replacement = '', x = .x)) %>%
  modify_at('treatment', factor, levels = c('Vehicle ctl (0.27% DMSO) ', 'SFN 15uM ', 'SFN-nitrile 15uM ', 'SFN-nitrile 45uM ', 'SFN-nitrile 135uM ')) %>%
  inset('timepoint', value = '48h')

qpcr_summary <- rbind(qpcr_24, qpcr_48) %>%
  pivot_longer(c('HO1', 'NQO1'), names_to = 'gene') %>%
  group_by(treatment, gene, timepoint) %>%
  summarise(mean = mean(value), 
            ser = se(value)) 

nqdata <- qpcr_summary %>%
  filter(gene == 'NQO1')

hodata <- qpcr_summary %>%
  filter(gene == 'HO1')

## sPLS Model --------------------------------------------------------------

ps_spls <- ps_counts %>% 
  filter_taxa(function(x) mean(x / sum(x)) > 1e-5, TRUE) %>%
  subset_samples(cohort %in% 3:8) %>%
  subset_samples(time %in% 0) %>%
  subset_samples(veg == 'broc')

clr_spls <- ps_spls %>%
  otu_table() %>%
  data.frame() %>%
  vegan::decostand('clr', pseudocount = 1)

#Metadata
meta_full <- data.frame(sample_data(ps_spls))

#Microbiome
micro_meta <- clr_spls %>%
  rownames_to_column('sample') %>%
  left_join(meta_full) %>%
  filter(time == 0) %>%
  filter(veg == 'broc')

#Urine SFN Data
urine_data <- urine_clean %>%
  dplyr::select(-sample) %>%
  filter(cohort %in% 3:8) %>%
  mutate(DTC_Tot = SFN + SFN_Cys + SFN_NAC + SFN_CG + SFN_GSH) %>%
  group_by(subject_id) %>%
  summarise(Tot_DTC = sum(DTC_Tot, na.rm = T),
            Tot_NIT = sum(SFN_NIT, na.rm = T), 
            Tot_SFN = sum(SFN_Tot, na.rm = T))

#Dietary Data
logs <- list.files('./Data/Metadata/Diet_logs/', full.names = T)
all_logs <- read_in_logs(list.files(here('Data/Metadata/Diet_logs/'), full.names = T)) 

#Diet Data
diet_summary <- all_logs %>%
  filter(relative_day %in% paste0('Day', 1:7)) %>%
  group_by(ID, nutrient) %>%
  summarise(mean_prestudy = mean(value, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = 'nutrient', values_from = 'mean_prestudy') %>%
  rename_with(~gsub('\\)', '', gsub('\\(', '', gsub(' ', '_', .x)))) %>%
  rename_with(~paste0('diet_', .x)) %>%
  rename('subject_id' = 'diet_ID') %>%
  mutate(across(starts_with('diet_'), ~ .x/diet_Cals_kcal)) %>%
  dplyr::select(-diet_Cals_kcal) 
  
#Remove redundant information
diet_small <- diet_summary %>%
  dplyr::select(-diet_Fib16_g, -diet_FatCals_kcal, -diet_SolFib16_g, -diet_SatCals_kcal, -diet_Carb_g, -diet_Fat_g, -diet_Folate_mcg)

#And the small version
full_block <- left_join(micro_meta, diet_small) %>%
  left_join(urine_data)

#Create the microbiome block
micro_block <- full_block %>%
  dplyr::select(starts_with('ASV')) %>%
  auto_scale()

#Create the diet blocks
diet_block <- full_block %>%
  dplyr::select(starts_with('diet_'), -diet_condensed) %>%
  log_transform() %>%
  auto_scale()

SFN_block <- full_block %>%
  dplyr::select(starts_with('Tot')) %>%
  log_transform() %>%
  auto_scale()

#Now the raw blocks for Spearman's rho:
micro_raw <- full_block %>%
  dplyr::select(starts_with('ASV')) 

diet_raw <- full_block %>%
  dplyr::select(starts_with('diet_'), -diet_condensed) 

SFN_raw <- full_block %>%
  dplyr::select(starts_with('Tot')) 

#Combine into a final list
blocks <- list(micro = micro_block, diet = diet_block, SFN = SFN_block)

model_final <- readRDS(here('Data/Plot_Data/final_sPLS.RDS'))

# Figure 1 ----------------------------------------------------------------

#Next we want to know the distribution of each metabolite in each tissue:
metab_col <- c('#EF476F', '#F8961E', '#06D6A0', '#118AB2', '#7400B8', '#073B4C')
names(metab_col) <- rev(unique(plasma_tidy$metabolite))
metabcol2 <- metab_col[c(1,4,2,5,3,6)]

plasma_dist <- plasma_tidy %>%
  #filter(time != 0) %>%
  group_by(metabolite, time) %>%
  summarise(mean = mean(uMol, na.rm = TRUE),
            ser = se(uMol)) %>%
  ungroup() %>%
  modify_at('metabolite', factor, levels = rev(c('SFN_NIT', 'SFN_NAC', 'SFN_GSH', 'SFN_Cys', 'SFN_CG', 'SFN')))

pdish <- ggplot(plasma_dist, aes(x = time, y = mean, fill = metabolite)) +
  geom_col() +
  theme_sfn() +
  scale_fill_manual(values = metabcol2, name = 'Metabolite') + 
  xlab('Time (Hours)') +
  ylab('Mean Total SFN Metabolites (µM)') +
  drop_legend +
  drop_gap +
  theme(axis.text.x = element_text(size = 35), 
        axis.text.y = element_text(size = 35),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size = 38),
        plot.title = element_text(size = 35, face = 'plain'),
        legend.title = element_text(size = 38, face = 'plain'),
        legend.text = element_text(size = 35)) +
  guides(fill = guide_legend(override.aes = list(size = 30)))

urine_dist <- urine_tidy %>%
  #filter(time != 0) %>%
  drop_na() %>%
  group_by(metabolite, time) %>%
  summarise(mean = mean(uMol, na.rm = TRUE),
            ser = se(uMol)) %>%
  ungroup() %>%
  modify_at('metabolite', factor, levels = rev(c('SFN_NIT', 'SFN_NAC', 'SFN_GSH', 'SFN_Cys', 'SFN_CG', 'SFN')))

udish <- ggplot(urine_dist, aes(x = time, y = mean, fill = metabolite)) +
  geom_col() +
  theme_sfn() +
  scale_fill_manual(values = metabcol2, 
                    name = 'Metabolite') +
  xlab('Time (Hours)') +
  ylab('Mean Total SFN Metabolites (µmol)') +
  drop_legend +
  drop_gap +
  theme(axis.text.x = element_text(size = 35),
        axis.text.y = element_text(size = 35),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size = 38),
        plot.title = element_text(size = 35, face = 'plain'),
        legend.title = element_text(size = 38, face = 'plain'),
        legend.text = element_text(size = 35)) +
  guides(shape = guide_legend(override.aes = list(size = 30)))

fecal_dist <- fecal_tidy %>%
  #filter(time != 0) %>%
  drop_na() %>%
  group_by(metabolite, time) %>%
  summarise(mean = mean(uMol, na.rm = T),
            ser = se(uMol)) %>%
  ungroup() %>%
  modify_at('metabolite', factor, levels = rev(c('SFN_NIT', 'SFN_NAC', 'SFN_GSH', 'SFN_Cys', 'SFN_CG', 'SFN')))


fdish <- ggplot(fecal_dist, aes(x = time, y = mean, fill = metabolite)) +
  geom_col() +
  theme_sfn() +
  scale_fill_manual(values = metabcol2, name = 'Metabolite',
                    guide = guide_legend(reverse = T),
                    labels = rev(c('SFN-NIT', 'SFN-NAC', 'SFN-GSH', 'SFN-Cys', 'SFN-CG', 'SFN'))) +
  xlab('Time (Hours)') +
  ylab('Mean Total SFN Metabolites (nmol/g)') +
  theme(legend.position = c(0.55,0.8)) +
  drop_gap  +
  theme(axis.text.x = element_text(size = 35),
        axis.text.y = element_text(size = 35),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size = 38),
        plot.title = element_text(size = 35, face = 'plain'),
        legend.title = element_text(size = 38, face = 'plain'),
        legend.text = element_text(size = 35)) +
  guides(fill = guide_legend(override.aes = list(size = 20)))

plot_grid(pdish, udish, fdish, nrow = 1, labels = c('A', 'B', 'C'), label_size = 40)

# Figure 2 ----------------------------------------------------------------

#Plasma - SFNNIT:
fixgap <- function(x,f = 5){
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,max(x, na.rm = TRUE)+sd(x, na.rm = TRUE)/f))
}

pntc2 <- ggplot(plasma_clean_sub, aes(x = time, y = SFN_NIT, group = subject_id, color = subject_id)) +
  geom_point(alpha = 1, size = 8) +
  #geom_path(alpha = 1, size = 1) +
  geom_path(aes(group = 1),stat = 'summary', fun = mean, color = 'black', alpha = 1, size = 5) +
  theme_sfn() +
  scale_color_manual(values = BSS_col, name = 'Subject') +
  xlab('Time (Hours)') +
  ylab('SFN-NIT (µM)') +
  fixgap(plasma_clean_sub$SFN_NIT) +
  drop_legend +
  theme(axis.text.x = element_text(size = 35),
        axis.text.y = element_text(size = 35),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size = 38),
        plot.title = element_text(size = 35, face = 'plain'),
        legend.title = element_text(size = 38, face = 'plain'),
        legend.text = element_text(size = 35)) +
  guides(fill = guide_legend(override.aes = list(size = 20))) #+
  #annotate('text', x = 5, y = 2, label = 'Plasma SFN-NIT', size = 20)

pctc2 <- ggplot(plasma_clean_sub, aes(x = time, y = SFN_Cys, group = subject_id, color = subject_id)) +
  geom_point(alpha = 1, size = 8) +
  #geom_path() +
  geom_path(aes(group = 1),stat = 'summary', fun = mean, color = 'black', alpha = 1, size = 5) +
  theme_sfn() +
  scale_color_manual(values = BSS_col, name = 'Subject') +
  xlab('Time (Hours)') +
  ylab('SFN-Cys (µM)') +
  fixgap(plasma_clean_sub$SFN_Cys) +
  drop_legend +
  theme(axis.text.x = element_text(size = 35),
        axis.text.y = element_text(size = 35),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size = 38),
        plot.title = element_text(size = 35, face = 'plain'),
        legend.title = element_text(size = 38, face = 'plain'),
        legend.text = element_text(size = 35)) +
  guides(fill = guide_legend(override.aes = list(size = 20))) #+

pnactc2 <- ggplot(plasma_clean_sub, aes(x = time, y = SFN_NAC, group = subject_id, color = subject_id)) +
  geom_point(alpha = 1, size = 8) +
  geom_path(aes(group = 1),stat = 'summary', fun = mean, color = 'black', alpha = 1, size = 5) +
  theme_sfn() +
  scale_color_manual(values = BSS_col, name = 'Subject') +
  xlab('Time (Hours)') +
  ylab('SFN-NAC (µM)') +
  fixgap(plasma_clean_sub$SFN_NAC)  +
  drop_legend +
  theme(axis.text.x = element_text(size = 35),
        axis.text.y = element_text(size = 35),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size = 38),
        plot.title = element_text(size = 35, face = 'plain'),
        legend.title = element_text(size = 38, face = 'plain'),
        legend.text = element_text(size = 35)) +
  guides(fill = guide_legend(override.aes = list(size = 20))) #+

#Urine - SFNNIT:
untc2 <- ggplot(urine_clean_sub, aes(x = time, y = SFN_NIT, group = subject_id, color = subject_id)) +
  geom_point(alpha = 1, size = 8) +
  geom_path(aes(group = 1),stat = 'summary', fun = mean, color = 'black', alpha = 1, size = 5) +
  #geom_path() +
  theme_sfn() +
  scale_color_manual(values = BSS_col, name = 'Subject') +
  xlab('Time (Hours)') +
  ylab('SFN-NIT (µmol)') +
  fixgap(urine_clean_sub$SFN_NIT) +
  drop_legend +
  theme(axis.text.x = element_text(size = 35),
        axis.text.y = element_text(size = 35),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size = 38),
        plot.title = element_text(size = 35, face = 'plain'),
        legend.title = element_text(size = 38, face = 'plain'),
        legend.text = element_text(size = 35)) +
  guides(fill = guide_legend(override.aes = list(size = 20))) #+


unactc2 <- ggplot(urine_clean_sub, aes(x = time, y = SFN_NAC, group = subject_id, color = subject_id)) +
  geom_point(alpha = 1, size = 8) +
  geom_path(aes(group = 1),stat = 'summary', fun = mean, color = 'black', alpha = 1, size = 5) +
  theme_sfn() +
  scale_color_manual(values = BSS_col, name = 'Subject') +
  xlab('Time (Hours)') +
  ylab('SFN-NAC (µmol)') +
  fixgap(urine_clean_sub$SFN_NAC) +
  #drop_legend +
  theme(axis.text.x = element_text(size = 35),
        axis.text.y = element_text(size = 35),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size = 38),
        plot.title = element_text(size = 35, face = 'plain'),
        legend.title = element_text(size = 38, face = 'plain'),
        legend.text = element_text(size = 35)) +
  guides(color = guide_legend(override.aes = list(size = 10),
                              ncol = 1)) #+

uleg <- get_legend(unactc2)
unactc2 <- unactc2 +
  drop_legend


fecal_dif2 <- metadata_fecal %>%
  modify_at(3:11, mdy_hm) %>%
  mutate(across(c(`0h_produced`, `24h_produced`, `48h_produced`, `72h_produced`), 
                ~as.numeric(difftime(.x, time_consumed, units = 'hours')))) %>%
  dplyr::select(where(is.character), where(is.numeric)) %>%
  pivot_longer(cols = ends_with('produced'), names_to = 'time', values_to = 'time_passed') %>%
  mutate(time = gsub('h_produced', '', time))

fecal_dif2$time %<>% factor(., levels = c(0,24,48,72))

fecal_time2 <- left_join(fecal_clean_sub, fecal_dif2, by = c('subject_id' = 'ID', 'time'))


fs2 <- ggplot(fecal_time2, aes(x = time_passed, y = SFN_Tot, color = subject_id)) +
  #geom_point(alpha = 0.7) +
  geom_point(alpha = 1, size = 1) +
  geom_path(alpha = 1, size = 1) +
  #geom_path(aes(group = 1),stat = 'summary', fun = mean, color = 'black', alpha = 1, size = 3) +
  theme_sfn() +
  scale_color_manual(values = BSS_col) +
  #geom_vline(xintercept = -24, linetype = 'dashed') +
  #geom_vline(xintercept = 0, linetype = 'dashed') +
  #geom_vline(xintercept = 24, linetype = 'dashed') +
  #geom_vline(xintercept = 48, linetype = 'dashed') +
  #geom_vline(xintercept = 72, linetype = 'dashed') +
  scale_x_continuous(breaks = c(-24,0,24,48,72), limits = c(-26, 75)) +
  #scale_x_continuous(breaks = c(-24,0,24,48,72)) +
  xlab('Time (Hours)') +
  ylab('Total SFN Metabolites (nmol/g)') +
  fixgap(fecal_time2$SFN_Tot) +
  drop_legend +
  theme(axis.text.x = element_text(size = 35),
        axis.text.y = element_text(size = 35),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size = 38),
        plot.title = element_text(size = 35, face = 'plain'),
        legend.title = element_text(size = 38, face = 'plain'),
        legend.text = element_text(size = 35)) +
  guides(color = guide_legend(override.aes = list(size = 20))) #+

ugh <- fecal_clean_sub %>%
  dplyr::select(subject_id, treatment, cohort, time, SFN_Tot) %>%
  pivot_wider(names_from = 'time', values_from = 'SFN_Tot') %>%
  magrittr::inset('-24', value = NA) %>%
  pivot_longer(cols = c('0','24','48','72','-24'), names_to = 'time', values_to = 'SFN_Tot') %>%
  modify_at('time', as.numeric)

#Fecal - SFN:
fntc <- ggplot(ugh, aes(x = time, y = SFN_Tot, group = subject_id, color = subject_id)) +
  #geom_bar(aes(group = time),stat = 'summary', fun = mean, fill = 'grey', color = 'black', alpha = 0.5) +
  #geom_point(alpha = 0.7) +
  geom_path(aes(group = 1),stat = 'summary', fun = mean, color = 'black', alpha = 1, size = 5) +
  theme_void() +
  scale_x_continuous(breaks = c(-24,0,24,48,72), limits = c(-26, 75)) +
  #scale_x_discrete(breaks = c(-24,0,24,48,72)) +
  scale_color_manual(values = BSS_col, name = 'Subject') +
  xlab('Time (Hours)') +
  ylab('Total Metabolites (nmol/g)') +
  fixgap(fecal_clean_sub$SFN_Tot) +
  drop_legend

fs22 <- ggplot(fecal_time2, aes(x = time_passed, y = SFN_Tot, color = subject_id)) +
  geom_point(alpha = 1, size = 8) +
  theme_sfn() +
  scale_color_manual(values = BSS_col) +
  #geom_vline(xintercept = -24, linetype = 'dashed') +
  #geom_vline(xintercept = 0, linetype = 'dashed') +
  #geom_vline(xintercept = 24, linetype = 'dashed') +
  #geom_vline(xintercept = 48, linetype = 'dashed') +
  #geom_vline(xintercept = 72, linetype = 'dashed') +
  scale_x_continuous(breaks = c(-24,0,24,48,72), limits = c(-26, 75)) +
  #scale_x_continuous(breaks = c(-24,0,24,48,72)) +
  xlab('Time (Hours)') +
  ylab('Total Metabolites (nmol/g)') +
  fixgap(fecal_time2$SFN_Tot) +
  drop_legend +
  theme(axis.text.x = element_text(size = 35),
        axis.text.y = element_text(size = 35),
        axis.title.x = element_text(size = 38),
        axis.title.y = element_text(size = 38),
        plot.title = element_text(size = 35, face = 'plain'),
        legend.title = element_text(size = 38, face = 'plain'),
        legend.text = element_text(size = 35)) +
  guides(color = guide_legend(override.aes = list(size = 20))) #+

fecal_aligned <- align_plots(fntc, fs22, align="hv", axis="tblr")
falgn <- ggdraw(fecal_aligned[[2]]) + draw_plot(fecal_aligned[[1]])


pplots2 <- plot_grid(pntc2, pctc2, pnactc2, ncol = 3, labels = c('A', 'B', 'C'), label_size =  40)
pbox <- ggplot(plasma_clean_sub, aes(x = time, y = SFN_NAC, group = subject_id, color = subject_id)) +
  theme_void() +
  annotate('rect', xmin = 0, xmax = 10, ymin = 0, ymax = 5, fill = 'grey87') +
  annotate('text', x = 5, y = 2.5, label = 'Plasma', size = 22) 
pfin <- plot_grid(pbox, pplots2, ncol = 1, align = 'hv', rel_heights = c(1,10))


fplots22 <- plot_grid(untc2, unactc2, falgn, ncol = 3, labels = c('D', 'E', 'F'))
uplots <- plot_grid(untc2, unactc2, ncol = 2, labels = c('D', 'E'), label_size = 40)
ubox <- ggplot(urine_clean_sub, aes(x = time, y = SFN_NAC, group = subject_id, color = subject_id)) +
  theme_void() +
  annotate('rect', xmin = 0, xmax = 10, ymin = 0, ymax = 5, fill = 'grey87') +
  annotate('text', x = 5, y = 2.5, label = 'Urine', size = 22) 
ufin <- plot_grid(ubox, uplots, ncol = 1, align = 'hv', rel_heights = c(1,10))


fbox <- ggplot(fecal_time2, aes(x = time_passed, y = SFN_Tot, color = subject_id)) +
  theme_void() +
  annotate('rect', xmin = 0, xmax = 10, ymin = 0, ymax = 5, fill = 'grey87') +
  annotate('text', x = 5, y = 2.5, label = 'Stool', size = 22) 
ffin <- plot_grid(fbox, falgn, ncol = 1, align = 'hv', rel_heights = c(1,10), labels = c('', 'F'), label_size = 40)
bfin <- plot_grid(ufin, ffin, nrow = 1, rel_widths = c(2, 1))
newplots <- plot_grid(pfin, bfin, ncol = 1)
plot_grid(newplots, uleg, ncol = 2, rel_widths = c(10,1))


# Figure 3 ----------------------------------------------------------------

hplot <- ggplot(hdac_summary, aes(x = treatment, y = mean)) +
  geom_errorbar(aes(ymin = mean - ser, ymax = mean + ser), size = 3.5, width = 0.5) +
  geom_col(position = 'dodge', color = 'black', fill = 'grey40') +
  cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.3)) +
  xlab('') +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, 'Set1'),
                    name = 'Time',
                    labels = c('24h', '48h'))  +
  ylab('Relative HDAC \n Activity') +
  scale_x_discrete(labels = c('Control (0.27% DMSO)', 'SFN 15µM', 'SFN-NIT 15µM', 'SFN-NIT 45µM', 'SFN-NIT 135µM')) +
  theme(legend.position = c(0.7,0.7),
        axis.text.y = element_text(size = 35),
        axis.text.x = element_text(size = 35),
        axis.title.y = element_text(size = 40),
        legend.text = element_text(size = 34), 
        legend.title = element_text(size = 38), 
        legend.key.size = unit(22, 'mm')) +
  annotate('text', x = 2.0, y = 0.95, label = '*', size = 42)
hplot



nplot <- ggplot(nqdata, aes(x = treatment, y = mean, group = timepoint, fill = timepoint)) +
  geom_errorbar(aes(ymin = mean - ser, ymax = mean + ser),size = 3.5, position = position_dodge(width = 0.9), width = 0.5) +
  geom_col(position = 'dodge', color = 'black') +
  cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 4.95)) +
  xlab('') +
  scale_fill_manual(values = c('grey85', 'grey40'),
                    name = 'Time',
                    labels = c('24h', '48h'))  +
  ylab('Relative NQO1 \n Expression') +
  scale_x_discrete(labels = c('Control (0.27% DMSO)', 'SFN 15µM', 'SFN-NIT 15µM', 'SFN-NIT 45µM', 'SFN-NIT 135µM')) +
  theme(legend.position = c(0.9, 0.7),
        axis.text.y = element_text(size = 35),
        axis.text.x = element_text(size = 35),
        axis.title.y = element_text(size = 40),
        legend.text = element_text(size = 34), 
        legend.title = element_text(size = 38), 
        legend.key.size = unit(22, 'mm')) +
  annotate('text', x = 1.78, y = 2.8, label = '*', size = 32) +
  annotate('text', x = 2.23, y = 4.5, label = '*', size = 32) 
nplot


hoplot <- ggplot(hodata, aes(x = treatment, y = mean, fill = timepoint)) +
  geom_errorbar(aes(ymin = mean - ser, ymax = mean + ser),size = 3.5, position = position_dodge(width = 0.9), width = 0.5) +
  geom_col(position = 'dodge', color = 'black') +
  cowplot::theme_cowplot() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 4.5)) +
  xlab('') +
  scale_fill_manual(values = c('grey85', 'grey40'),
                    name = 'Time',
                    labels = c('24h', '48h'))  +
  ylab('Relative HO-1 \n Expression') +
  scale_x_discrete(labels = c('Control (0.27% DMSO)', 'SFN 15µM', 'SFN-NIT 15µM', 'SFN-NIT 45µM', 'SFN-NIT 135µM')) +
  theme(legend.position = c(0.9,0.7),
        axis.text.y = element_text(size = 35),
        axis.text.x = element_text(size = 35),
        axis.title.y = element_text(size = 40),
        legend.text = element_text(size = 34), 
        legend.title = element_text(size = 38), 
        legend.key.size = unit(22, 'mm')) +
  annotate('text', x = 1.78, y = 4.16, label = '*', size = 32) + 
  annotate('text', x = 2.23, y = 3.1, label = '*', size = 32) 
hoplot


cowplot::plot_grid(hplot, nplot, hoplot, nrow = 3, labels = c('A', 'B', 'C'), label_size = 40)


# Figure 4 ----------------------------------------------------------------

alltax <- data.frame(tax_table(ps_spls)) %>%
  rownames_to_column('ASV') %>%
  mutate(spasv = paste0(Genus, '_', Species))

final_ASV <- selectVar(model_final)$micro$value %>%
  filter(abs(value.var) >= 0.1) %>%
  rownames()

final_diet <- selectVar(model_final)$diet$value %>%
  filter(abs(value.var) >= 0.1) %>%
  rownames()

ASVtest <- micro_raw %>%
  dplyr::select(all_of(final_ASV)) 

diettest <- diet_raw %>%
  dplyr::select(all_of(final_diet)) 

#ASV to SFN
data.frame(cor(ASVtest, SFN_raw, method = 'sp')) %>%
  rownames_to_column('ASV') %>% 
  left_join(alltax) %>%
  column_to_rownames('spasv') %>%
  dplyr::select(where(is.numeric)) %>%
  pheatmap::pheatmap(display_numbers = F, labels_col = c('Total DTCs', 'SFN-NIT', 'All Metabolites'), angle_col = 270,
                     fontsize_row = 35, fontsize_col = 40, legend = F)

ah2 <- data.frame(cor(ASVtest, SFN_raw, method = 'sp')) %>%
  rownames_to_column('ASV') %>% 
  left_join(alltax) %>%
  column_to_rownames('spasv') %>%
  dplyr::select(where(is.numeric)) %>%
  pheatmap::pheatmap(display_numbers = F, labels_col = c('Total DTCs', 'SFN-NIT', 'All Metabolites'), angle_col = 270,
                     fontsize_row = 35, fontsize_col = 40, legend = T)

ggdraw(ah2$gtable$grob[[6]])


# Figure S1 ---------------------------------------------------------------

library(vegan)
dps <- ps_counts

#Make the rarefaction curves
d <- dps %>%  
  otu_table()

S <- specnumber(d) # observed number of species
(raremax <- min(rowSums(d)))
Srare <- rarefy(d, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
out <- rarecurve(d, step = 20, sample = raremax, col = "blue", cex = 0.6)

rare <- lapply(out, function(x){b <- as.data.frame(x) 
b <- data.frame(clase = b[,1], raw.read = rownames(b)) 
b$raw.read <- as.numeric(gsub("N", "",  b$raw.read)) 
return(b)})

#convert to data frame:
rare <- map_dfr(rare, function(x){
  z <- data.frame(x) 
  return(z)
}, .id = "Sample") 

#Pull metadata
md <- data.frame(sample_data(dps)) %>%
  dplyr::select(treatment, subject_id, sample)

rare %<>% modify_at('Sample', as.integer)

#Pull out the data to plot
uh <- rare %>% 
  dplyr::group_by(Sample) %>% 
  nest() %>%
  cbind(md) %>%
  unnest(cols = 'data') %>%
  modify_at('treatment', as.character) %>%
  mutate(veg = ifelse(treatment %in% c('BL', 'BU'), 'Broccoli', 'Alfalfa'))

#Plot that beautiful data
ggplot(uh, aes(x=raw.read, y=clase, color = veg, group = Sample))+
  geom_path(size = 3) + 
  labs(x = 'Sample Size', y = 'Species') +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c('#D90368', '#04A777'),
                    labels = c('Alfalfa', 'Broccoli'),
                    name = 'Vegetable \n Type') +
  theme(legend.position = c(0.8,0.3),
        axis.text.y = element_text(size = 35),
        axis.text.x = element_text(size = 35),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        legend.text = element_text(size = 34), 
        legend.title = element_text(size = 38), 
        legend.key.size = unit(22, 'mm')) 

# Figure S2 ---------------------------------------------------------------

#Broccoli vs Alfalfa - 24hr
ps_ab <- ps_counts %>%
  #subset_samples(treatment %in% c('BL', 'BU', 'AL', 'AU')) %>%
  subset_samples(time %in% c(24)) %>%
  subset_samples(cohort %in% 3:8)
ordBC_ab <- ordinate(ps_ab, method = 'PCoA', distance = 'jaccard')
vegplot <- plot_ordination(ps_ab, ordBC_ab, color = 'veg') +
  geom_point(size = 6) +
  #ggtitle('Broccoli vs Alfalfa - 24hr') +
  theme_cowplot() +
  scale_color_manual(values = c('#D90368', '#04A777'),
                    labels = c('Alfalfa', 'Broccoli'),
                    name = 'Vegetable \n Type') +
  theme(axis.text.y = element_text(size = 35),
        axis.text.x = element_text(size = 35),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        legend.text = element_text(size = 34), 
        legend.title = element_text(size = 38), 
        legend.key.size = unit(22, 'mm'))

ps_ball <- ps_counts %>%
  subset_samples(veg %in% c('broc')) %>%
  subset_samples(time %in% c(0,24, 48, 72)) %>%
  subset_samples(cohort %in% 3:8)
ordBC_ball <- ordinate(ps_ball, method = 'PCoA', distance = 'jaccard')
timeplot <- plot_ordination(ps_ball, ordBC_ball, color = 'time') +
  geom_point(size = 6) +
  #ggtitle('Broccoli - 0hr vs 24hr') +
  theme_cowplot() +
  scale_color_manual(values = c('#2176FF', '#F79824', '#51CB20', '#730071'),
                    labels = c('0 Hour', '24 Hour', '48 Hour', '72 Hour'),
                    name = 'Time Post \n Consumption') +
  theme(axis.text.y = element_text(size = 35),
        axis.text.x = element_text(size = 35),
        axis.title.y = element_text(size = 40),
        axis.title.x = element_text(size = 40),
        legend.text = element_text(size = 34), 
        legend.title = element_text(size = 38), 
        legend.key.size = unit(22, 'mm'))

cowplot::plot_grid(vegplot, timeplot, ncol = 2, labels = c('A', 'B'), label_size = 40)



# Figure S3 ---------------------------------------------------------------


data.frame(cor(ASVtest, diettest, method = 'sp')) %>%
  rownames_to_column('ASV') %>% 
  left_join(alltax) %>%
  column_to_rownames('spasv') %>%
  dplyr::select(where(is.numeric)) %>%
  pheatmap::pheatmap(display_numbers = F, fontsize_col = 40, fontsize_row = 35, legend = F)

meh <- data.frame(cor(ASVtest, diettest, method = 'sp')) %>%
  rownames_to_column('ASV') %>% 
  left_join(alltax) %>%
  column_to_rownames('spasv') %>%
  dplyr::select(where(is.numeric)) %>%
  pheatmap::pheatmap(display_numbers = F, fontsize_col = 40, fontsize_row = 35, legend = T)

ggdraw(meh$gtable$grob[[6]])

# Figure S4 ---------------------------------------------------------------

cloTest <- micro_raw %>%
  dplyr::select(ASV125, ASV220, ASV276) 

data.frame(cor(cloTest, SFN_raw, method = 'sp')) %>%
  rownames_to_column('ASV') %>% 
  left_join(alltax) %>%
  column_to_rownames('spasv') %>%
  dplyr::select(where(is.numeric)) %>%
  pheatmap::pheatmap(display_numbers = F, labels_col = c('Total DTCs', 'SFN-NIT', 'All Metabolites'), angle_col = 270,
                     fontsize_row = 35, fontsize_col = 40, legend = F)

ah3 <- data.frame(cor(cloTest, SFN_raw, method = 'sp')) %>%
  rownames_to_column('ASV') %>% 
  left_join(alltax) %>%
  column_to_rownames('spasv') %>%
  dplyr::select(where(is.numeric)) %>%
  pheatmap::pheatmap()

ggdraw(ah3$gtable$grob[[6]])


# Figure S5 ---------------------------------------------------------------


ASV_SFN <- cbind(ASVtest, SFN_raw)

makeCorPlot <- function(df, x, y, xlab = NULL){
  xOut <- df[,x]
  yOut <- df[,y]
  fr <- formula(paste0(quo_name(y), '~', quo_name(x)))
  temp_lm <- coef(lm(fr, df))
  temp_cor <- round(cor(xOut, yOut, method = 'sp'), 2)
  plot <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], color = .data[[y]])) +
    geom_point(size = 3) +
    geom_abline(intercept = temp_lm[1], slope = temp_lm[2]) +
    scale_color_gradient(high = 'green', low = 'blue', 
                         name = switch(y, 'Tot_SFN' = 'All Metabolites \n (µMol)', 'Tot_DTC' = 'Total DTCs \n (µMol)', 'Tot_NIT' = 'SFN-NIT \n (µMol)')) +
    cowplot::theme_cowplot() +
    ylab(switch(y, 'Tot_SFN' = 'All Metabolites (µMol)', 'Tot_DTC' = 'Total DTCs (µMol)', 'Tot_NIT' = 'SFN-NIT (µMol)')) +
    xlab(ifelse(is.null(xlab), x, xlab))
  return(plot)
}


## Roseburia ASV112 --------------------------------------------------------
roseplot <- makeCorPlot(ASV_SFN, 'ASV112', 'Tot_SFN', xlab = 'Unannotated Roseburia (ASV112)')

## Dorea ASV28 -------------------------------------------------------------
doreaplot <- makeCorPlot(ASV_SFN, 'ASV28', 'Tot_SFN', xlab = 'Dorea longicatena')

## Alistipes ASV313 --------------------------------------------------------
aliplot <- makeCorPlot(ASV_SFN, 'ASV313', 'Tot_DTC', xlab = 'Unannotated Alistipes (ASV313)')

## Blautia ASV145 ----------------------------------------------------------
blautiaplot <- makeCorPlot(ASV_SFN, 'ASV145', 'Tot_SFN', xlab = 'Unannotated Blautia (ASV145)')

## Bifidobacterium ASV247 --------------------------------------------------
bifidoplot <- makeCorPlot(ASV_SFN, 'ASV247', 'Tot_NIT', xlab = 'Unannotated Bifidobacterium (ASV247)')

## Ruminococcus ASV527 -----------------------------------------------------
ruminoplot <- makeCorPlot(ASV_SFN, 'ASV527', 'Tot_NIT', xlab = 'Unannotated Ruminococcus \n Torques Group (ASV527)')

cowplot::plot_grid(roseplot, doreaplot, blautiaplot, bifidoplot, ruminoplot, aliplot, labels = c(LETTERS[1:6]))



