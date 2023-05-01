#Analysis of microbiome data for BSS-targeted SFN Analysis
#Analysis completed by Yanni Bouranis, dada2 preprocessing completed by Ed Davis
# Environment Setup -------------------------------------------------------
library(tidyverse)
library(phyloseq)
library(magrittr)
library(vegan)
library(cowplot)
library(here)

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

set.seed(120)
setwd('~/Documents/Manuscripts/BSS_Nitriles/Bouranis_SFN_Microbiome/')
here::i_am('./Microbiome_Analysis_Script.R')


# Data Load ---------------------------------------------------------------
#Prep data for loading into phyloseq
asvtab <- readRDS(here('./Data/Microbiome/seqtab_nochim.rds'))
taxtab <- readRDS(here('./Data/Microbiome/tax.rds'))
metadata <- as.data.frame(read_csv(here('./Data/Microbiome/metadata_microbiome.csv')))
metadata$time %<>% factor(levels = c(0, 3, 24, 48, 72))
metadata$cohort %<>% factor()
#Load in information on how much sprouts each group ate
brocinfo <- read_csv(here('./Data/Metadata/sprout_meta.csv')) %>%
  modify_at('cohort', as.factor)
#Load in SFN-info
load(here('./Data/Microbiome/SFN_Targeted.RData'))
rownames(metadata) <- metadata$sample
rownames(asvtab) <- metadata$sample
ps_raw <- phyloseq(otu_table(asvtab, taxa_are_rows = FALSE),
               sample_data(metadata),
               tax_table(taxtab))
#save(ps_raw, urine_clean, plasma_clean, fecal_clean, metadata, brocinfo, file = 'CompiledData.RData')

#Load in the relevant data:
#load(here('human_analysis/Targeted_SFN/Final_Scripts/Data/CompiledData.RData'))

# Microbiome Cleaning -----------------------------------------------------

#Give arbitrary names to the taxa as opposed to keeping as just DNA-sequences which identify them
taxa_names(ps_raw) <- paste0("ASV", seq(ntaxa(ps_raw)))

#Fix the metadata
mdata_veg <- data.frame(sample_data(ps_raw)) %>%
  mutate(veg = ifelse(treatment %in% c('BL', 'BU'), 'broc', 'alf')) 
sample_data(ps_raw) <- mdata_veg

#Fill in missing genus names:
renames <- rownames(tax_table(ps_raw)[is.na(tax_table(ps_raw)[, 'Genus'])])
taxdf <- tax_table(ps_raw)[renames,]
renamed_genus <- unname(sapply(taxa_names(taxdf), function(x) paste0('f_', taxdf[x, 'Family'], '_', x)))
tax_table(ps_raw)[renames, 'Genus'] <- renamed_genus

#Fill in missing species names
renames_sp <- rownames(tax_table(ps_raw)[is.na(tax_table(ps_raw)[, 'Species'])])
taxdf_sp <- tax_table(ps_raw)[renames_sp,]
renamed_species <- unname(sapply(taxa_names(taxdf_sp), function(x) paste0('g_', taxdf_sp[x, 'Genus'], '_', x)))
tax_table(ps_raw)[renames_sp, 'Species'] <- renamed_species

#ps_raw has 4060 ASVs
#Remove taxa that do not appear at least 3 times in at least 20% of all samples
ps_counts <- ps_raw %>% filter_taxa(function(x) sum(x > 3) > (0.2*length(x)), TRUE) #317
#Convert from counts to relative abundance
ps_relab <- ps_counts %>% transform_sample_counts(function(x) x / sum(x) ) #317
#Filter out low abundance (>1e-5) taxa
ps <- ps_relab %>% filter_taxa(function(x) mean(x) > 1e-5, TRUE) #317

ps

# Alpha Diversity Analysis ------------------------------------------------

adiv <- ps_raw %>%
  rarefy_even_depth(rngseed = 12) %>%
  subset_samples(cohort %in% 3:8) %>%
  subset_samples(time %in% c(0,24,48,72)) %>%
  estimate_richness(measures = c('Observed', 'Shannon', 'Simpson')) %>%
  rownames_to_column('sample') %>%
  left_join(.,as.data.frame(sample_data(ps_raw))) %>%
  left_join(brocinfo) %>%
  pivot_longer(cols = c('Observed', 'Shannon', 'Simpson'), names_to = 'measure') %>%
  #mutate(veg = ifelse(treatment %in% c('BL', 'BU'), 'broc', 'alf')) %>%
  drop_na() %>%
  group_by(measure) %>%
  nest() %>%
  #Run all three models 
  mutate(model = purrr::map(data, function(x) lmerTest::lmer(value ~ time*veg + grams_fed + (1|subject_id), data = x))) %>%
  #Contrast 1: For broccoli, 0 vs 24hr
  mutate(broc_0v24 = map_dbl(model, function(x) lmerTest::contest(x, c(0,1,0,0,0,0,1,0,0))$`Pr(>F)`)) %>%
  mutate(broc_0v48 = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,1,0,0,0,0,1,0))$`Pr(>F)`)) %>%
  mutate(broc_0v72 = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,0,1,0,0,0,0,1))$`Pr(>F)`)) %>%
  mutate(broc_24v48 = map_dbl(model, function(x) lmerTest::contest(x, c(0,1,-1,0,0,0,1,-1,0))$`Pr(>F)`)) %>%
  mutate(broc_24v72 = map_dbl(model, function(x) lmerTest::contest(x, c(0,1,0,-1,0,0,1,0,-1))$`Pr(>F)`)) %>%
  mutate(broc_48v72 = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,1,-1,0,0,0,1,-1))$`Pr(>F)`)) %>%
  mutate(alf_0v24 = map_dbl(model, function(x) lmerTest::contest(x, c(0,1,0,0,0,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(alf_0v48 = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,1,0,0,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(alf_0v72 = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,0,1,0,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(alf_24v48 = map_dbl(model, function(x) lmerTest::contest(x, c(0,1,-1,0,0,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(alf_24v72 = map_dbl(model, function(x) lmerTest::contest(x, c(0,1,0,-1,0,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(alf_48v72 = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,1,-1,0,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(bva_0h = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,0,0,1,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(bva_24h = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,0,0,1,0,1,0,0))$`Pr(>F)`)) %>%
  mutate(bva_48h = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,0,0,1,0,0,1,0))$`Pr(>F)`)) %>%
  mutate(bva_72h = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,0,0,1,0,0,0,1))$`Pr(>F)`))  

adiv %>%
  column_to_rownames('measure') %>%
  dplyr::select(-data, -model) %>%
  apply(., 2, function(x) sum(x < 0.05)) 

adiv_adjusted <- adiv %>%
  ungroup() %>%
  pivot_longer(where(is.numeric), names_to = 'test', values_to = 'pval') %>%
  mutate(padj = p.adjust(pval, method = 'BH')) %>%
  dplyr::select(-data, -model, -pval) %>%
  pivot_wider(names_from = 'test', values_from = 'padj')

adiv_adjusted %>%
  column_to_rownames('measure') %>%
  apply(., 2, function(x) sum(x < 0.05)) 

# Beta Diversity Analysis -------------------------------------------------

#Cohort Effect
ps_co <- ps %>%
  subset_samples(time %in% c(0)) %>%
  subset_samples(cohort %in% 3:8)
ordBC_co <- ordinate(ps_co, method = 'PCoA', distance = 'bray')
plot_ordination(ps_co, ordBC_co, color = 'cohort', shape = 'treatment') +
  geom_point(size = 2) +
  ggtitle('Broccoli - 0h, L vs UL')
distmat_co <- phyloseq::distance(ps_co, method = 'bray')
mdata_co <- data.frame(sample_data(ps_co))%>%
  left_join(brocinfo)
#Run betadisper to verify the distrubtion of our groups is equal, an underlying assumption of PERMANOVA
bdisp_co <- betadisper(distmat_co, mdata_co$cohort)
#Evaluate using a permutation test
permutest(bdisp_co) #p = 0.95
adonis2(distmat_co~cohort, data = mdata_co) #p = .297
#No changes
adonis2(distmat_co~grams_fed, data = mdata_co) #p = .664

#Broccoli vs Alfalfa - 0hr
ps_ab0h <- ps_counts %>%
  #subset_samples(treatment %in% c('BL', 'BU', 'AL', 'AU')) %>%
  subset_samples(time %in% c(0)) %>%
  subset_samples(cohort %in% 3:8)
ordBC_ab0h <- ordinate(ps_ab0h, method = 'PCoA', distance = 'jaccard')
plot_ordination(ps_ab0h, ordBC_ab0h, color = 'veg') +
  geom_point(size = 2) +
  #ggtitle('Broccoli vs Alfalfa - 24hr') +
  theme_cowplot() +
  scale_color_manual(values = c('#D90368', '#04A777'),
                    labels = c('Alfalfa', 'Broccoli'),
                    name = 'Vegetable Type')
distmat_ab0h <- phyloseq::distance(ps_ab0h, method = 'jaccard')
mdata <- data.frame(sample_data(ps_ab0h)) 
#Run betadisper to verify the distrubtion of our groups is equal, an underlying assumption of PERMANOVA
bdisp_ab0h <- betadisper(distmat_ab0h, mdata$veg)
#Evaluate using a permutation test
permutest(bdisp_ab0h) #p = 0.211
adonis2(distmat_ab0h~veg, data = mdata) #p = 0.967 
#No changes

#Broccoli vs Alfalfa - 24hr
ps_ab <- ps_counts %>%
  #subset_samples(treatment %in% c('BL', 'BU', 'AL', 'AU')) %>%
  subset_samples(time %in% c(24)) %>%
  subset_samples(cohort %in% 3:8)
ordBC_ab <- ordinate(ps_ab, method = 'PCoA', distance = 'jaccard')
plot_ordination(ps_ab, ordBC_ab, color = 'veg') +
  geom_point(size = 2) +
  #ggtitle('Broccoli vs Alfalfa - 24hr') +
  theme_cowplot() +
  scale_color_manual(values = c('#D90368', '#04A777'),
                    labels = c('Alfalfa', 'Broccoli'),
                    name = 'Vegetable Type')
distmat_ab <- phyloseq::distance(ps_ab, method = 'jaccard')
mdata <- data.frame(sample_data(ps_ab)) 
#Run betadisper to verify the distrubtion of our groups is equal, an underlying assumption of PERMANOVA
bdisp_ab <- betadisper(distmat_ab, mdata$veg)
#Evaluate using a permutation test
permutest(bdisp_ab) #p = 0.211
adonis2(distmat_ab~veg, data = mdata) #p = 0.967 
#No changes

#Broccoli vs Alfalfa - 48hr
ps_ab48h <- ps_counts %>%
  #subset_samples(treatment %in% c('BL', 'BU', 'AL', 'AU')) %>%
  subset_samples(time %in% c(48)) %>%
  subset_samples(cohort %in% 3:8)
ordBC_ab48h <- ordinate(ps_ab48h, method = 'PCoA', distance = 'jaccard')
plot_ordination(ps_ab48h, ordBC_ab48h, color = 'veg') +
  geom_point(size = 2) +
  #ggtitle('Broccoli vs Alfalfa - 24hr') +
  theme_cowplot() +
  scale_color_manual(values = c('#D90368', '#04A777'),
                    labels = c('Alfalfa', 'Broccoli'),
                    name = 'Vegetable Type')
distmat_ab48h <- phyloseq::distance(ps_ab48h, method = 'jaccard')
mdata <- data.frame(sample_data(ps_ab48h)) 
#Run betadisper to verify the distrubtion of our groups is equal, an underlying assumption of PERMANOVA
bdisp_ab48h <- betadisper(distmat_ab48h, mdata$veg)
#Evaluate using a permutation test
permutest(bdisp_ab48h) #p = 0.211
adonis2(distmat_ab48h~veg, data = mdata) #p = 0.967 
#No changes

#Broccoli vs Alfalfa - 72hr
ps_ab72h <- ps_counts %>%
  #subset_samples(treatment %in% c('BL', 'BU', 'AL', 'AU')) %>%
  subset_samples(time %in% c(72)) %>%
  subset_samples(cohort %in% 3:8)
ordBC_ab72h <- ordinate(ps_ab72h, method = 'PCoA', distance = 'jaccard')
plot_ordination(ps_ab72h, ordBC_ab72h, color = 'veg') +
  geom_point(size = 2) +
  #ggtitle('Broccoli vs Alfalfa - 24hr') +
  theme_cowplot() +
  scale_color_manual(values = c('#D90368', '#04A777'),
                    labels = c('Alfalfa', 'Broccoli'),
                    name = 'Vegetable Type')
distmat_ab72h <- phyloseq::distance(ps_ab72h, method = 'jaccard')
mdata <- data.frame(sample_data(ps_ab72h)) 
#Run betadisper to verify the distrubtion of our groups is equal, an underlying assumption of PERMANOVA
bdisp_ab72h <- betadisper(distmat_ab72h, mdata$veg)
#Evaluate using a permutation test
permutest(bdisp_ab72h) #p = 0.211
adonis2(distmat_ab72h~veg, data = mdata) #p = 0.967 
#No changes

#Broccoli: 0 vs 24hr
ps_b <- ps_counts %>%
  subset_samples(veg %in% c('broc')) %>%
  subset_samples(time %in% c(0,24)) %>%
  subset_samples(cohort %in% 3:8)
ordBC_b <- ordinate(ps_b, method = 'PCoA', distance = 'jaccard')
plot_ordination(ps_b, ordBC_b, color = 'time') +
  geom_point(size = 2) +
  #ggtitle('Broccoli - 0hr vs 24hr') +
  theme_cowplot() +
  scale_color_manual(values = c('#2176FF', '#F79824'),
                    labels = c('0 Hour', '24 Hour'),
                    name = 'Time Post \n Consumption')
distmat_b <- phyloseq::distance(ps_b, method = 'jaccard')
mdata <- data.frame(sample_data(ps_b)) 
#Run betadisper to verify the distrubtion of our groups is equal, an underlying assumption of PERMANOVA
bdisp_b <- betadisper(distmat_b, mdata$time)
#Evaluate using a permutation test
permutest(bdisp_b) #p = 0.539
adonis2(distmat_b~time, data = mdata) #p = 0.967 
#No changes

#Alfalfa: 0 vs 24hr
ps_a <- ps_counts %>%
  subset_samples(veg %in% c('alf')) %>%
  subset_samples(time %in% c(0,24)) %>%
  subset_samples(cohort %in% 3:8)
ordBC_a <- ordinate(ps_a, method = 'PCoA', distance = 'jaccard')
plot_ordination(ps_a, ordBC_a, color = 'time') +
  geom_point(size = 2) +
  #ggtitle('Broccoli - 0hr vs 24hr') +
  theme_cowplot() +
  scale_color_manual(values = c('#2176FF', '#F79824'),
                    labels = c('0 Hour', '24 Hour'),
                    name = 'Time Post \n Consumption')
distmat_a <- phyloseq::distance(ps_a, method = 'jaccard')
mdata <- data.frame(sample_data(ps_a)) 
#Run betadisper to verify the distrubtion of our groups is equal, an underlying assumption of PERMANOVA
bdisp_a <- betadisper(distmat_a, mdata$time)
#Evaluate using a permutation test
permutest(bdisp_a) #p = 0.539
adonis2(distmat_a~time, data = mdata) #p = 0.967 
#No changes

#Broccoli: All Timepoints 
ps_ball <- ps_counts %>%
  subset_samples(veg %in% c('broc')) %>%
  subset_samples(time %in% c(0,24, 48, 72)) %>%
  subset_samples(cohort %in% 3:8)
ordBC_ball <- ordinate(ps_ball, method = 'PCoA', distance = 'jaccard')
plot_ordination(ps_ball, ordBC_ball, color = 'time') +
  geom_point(size = 2) +
  #ggtitle('Broccoli - 0hr vs 24hr') +
  theme_cowplot() +
  scale_color_manual(values = c('#2176FF', '#F79824', '#51CB20', '#730071'),
                    labels = c('0 Hour', '24 Hour', '48 Hour', '72 Hour'),
                    name = 'Time Post \n Consumption')
distmat_ball <- phyloseq::distance(ps_ball, method = 'jaccard')
mdata <- data.frame(sample_data(ps_ball)) 
#Run betadisper to verify the distrubtion of our groups is equal, an underlying assumption of PERMANOVA
bdisp_ball <- betadisper(distmat_ball, mdata$time)
#Evaluate using a permutation test
permutest(bdisp_ball) #p = 0.539
adonis2(distmat_ball~time, data = mdata) #p = 0.967 
#No changes

#Alfalfa: All Timepoints 
ps_alfall <- ps_counts %>%
  subset_samples(veg %in% c('alf')) %>%
  subset_samples(time %in% c(0,24, 48, 72)) %>%
  subset_samples(cohort %in% 3:8)
ordBC_alfall <- ordinate(ps_alfall, method = 'PCoA', distance = 'jaccard')
plot_ordination(ps_alfall, ordBC_alfall, color = 'time') +
  geom_point(size = 2) +
  #ggtitle('Broccoli - 0hr vs 24hr') +
  theme_cowplot() +
  scale_color_manual(values = c('#2176FF', '#F79824', '#51CB20', '#730071'),
                    labels = c('0 Hour', '24 Hour', '48 Hour', '72 Hour'),
                    name = 'Time Post \n Consumption')
distmat_alfall <- phyloseq::distance(ps_alfall, method = 'jaccard')
mdata <- data.frame(sample_data(ps_alfall)) 
#Run betadisper to verify the distrubtion of our groups is equal, an underlying assumption of PERMANOVA
bdisp_alfall <- betadisper(distmat_alfall, mdata$time)
#Evaluate using a permutation test
permutest(bdisp_alfall) #p = 0.539
adonis2(distmat_alfall~time, data = mdata) #p = 0.967 
#No changes


# Differential Abundance Analysis ------------------------------------------

#Differential Abundance Testing
ps_c <- ps_counts %>% 
  filter_taxa(function(x) mean(x / sum(x)) > 1e-5, TRUE) %>%
  subset_samples(cohort %in% 3:8) %>%
  subset_samples(time %in% c(0,24,48,72)) %>%
  rarefy_even_depth(rngseed = 12) 
 
#Form the matrix of ASV counts
microdata <- ps_c %>%
  otu_table() %>%
  as.data.frame() 

#CLR transform the data
micro_clr <- microdata %>%
  vegan::decostand('clr', pseudocount = 1)

#Make nice and tidy
micro_tidy <- micro_clr %>%
  rownames_to_column('sample') %>%
  pivot_longer(cols = starts_with('ASV'), names_to = 'ASV', values_to = 'clr') %>%
  left_join(., as.data.frame(sample_data(ps_c))) %>%
  left_join(brocinfo)  

micro_lm <- micro_tidy %>%
  group_by(ASV) %>%
  nest() %>%
  mutate(model = purrr::map(data, function(x) lmerTest::lmer(clr ~ time*veg + grams_fed + (1|subject_id), data = x))) %>%
  #Contrast 1: For broccoli, 0 vs 24hr
  mutate(broc_0v24 = map_dbl(model, function(x) lmerTest::contest(x, c(0,1,0,0,0,0,1,0,0))$`Pr(>F)`)) %>%
  mutate(broc_0v48 = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,1,0,0,0,0,1,0))$`Pr(>F)`)) %>%
  mutate(broc_0v72 = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,0,1,0,0,0,0,1))$`Pr(>F)`)) %>%
  mutate(broc_24v48 = map_dbl(model, function(x) lmerTest::contest(x, c(0,1,-1,0,0,0,1,-1,0))$`Pr(>F)`)) %>%
  mutate(broc_24v72 = map_dbl(model, function(x) lmerTest::contest(x, c(0,1,0,-1,0,0,1,0,-1))$`Pr(>F)`)) %>%
  mutate(broc_48v72 = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,1,-1,0,0,0,1,-1))$`Pr(>F)`)) %>%
  mutate(alf_0v24 = map_dbl(model, function(x) lmerTest::contest(x, c(0,1,0,0,0,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(alf_0v48 = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,1,0,0,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(alf_0v72 = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,0,1,0,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(alf_24v48 = map_dbl(model, function(x) lmerTest::contest(x, c(0,1,-1,0,0,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(alf_24v72 = map_dbl(model, function(x) lmerTest::contest(x, c(0,1,0,-1,0,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(alf_48v72 = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,1,-1,0,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(bva_0h = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,0,0,1,0,0,0,0))$`Pr(>F)`)) %>%
  mutate(bva_24h = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,0,0,1,0,1,0,0))$`Pr(>F)`)) %>%
  mutate(bva_48h = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,0,0,1,0,0,1,0))$`Pr(>F)`)) %>%
  mutate(bva_72h = map_dbl(model, function(x) lmerTest::contest(x, c(0,0,0,0,1,0,0,0,1))$`Pr(>F)`))  

micro_lm %>%
  column_to_rownames('ASV') %>%
  dplyr::select(where(is.numeric)) %>%
  apply(., 2, function(x) sum(x < 0.05)) 

lm_adjusted <- micro_lm %>%
  ungroup() %>%
  pivot_longer(where(is.numeric), names_to = 'test', values_to = 'pval') %>%
  mutate(padj = p.adjust(pval, method = 'BH')) %>%
  dplyr::select(-data, -model, -pval) %>%
  pivot_wider(names_from = 'test', values_from = 'padj')

lm_adjusted %>%
  column_to_rownames('ASV') %>%
  apply(., 2, function(x) sum(x < 0.05)) 


# Multi-Omic Data Prep ----------------------------------------------------

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
logs <- list.files(here('./Data/Metadata/Diet_logs/', full.names = T))
all_logs <- read_in_logs(list.files(here('./Data/Metadata/Diet_logs/'), full.names = T)) 

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

# sPLS Model --------------------------------------------------------------

## Model Tuning ------------------------------------------------------------
library(mixOmics)


#Model Design Tuning:

#Exploratory Pairwise Analysis:
keepX <- c(25,25)
keepY <- c(25,25)

#Diet and Microbiome
pls_dm <- spls(blocks$diet, blocks$micro, keepX = keepX, keepY = keepY, ncomp = 2, mode = 'canonical')
plotVar(pls_dm, cutoff = 0.5, title = 'Diet vs Microbiome', legend = c('Diet', 'Microbiome'), 
        var.names = F, pch = c(16,18), col = c('darkorchid', 'brown')) #Both positive and negative correlations
cor(pls_dm$variates$X, pls_dm$variates$Y) #0.83, 0.84

#Microbiome and SFN
pls_ms <- spls(blocks$micro, blocks$SFN, keepX = keepX,  ncomp = 2, mode = 'canonical')
plotVar(pls_ms, cutoff = 0.5, title = 'Microbiome vs SFN', legend = c('Microbiome', 'GLSHP'), 
        var.names = F, pch = c(18,17), col = c('brown', 'lightgreen')) #Not super strong
cor(pls_ms$variates$X, pls_ms$variates$Y) #0.83, 0.85

#Diet and SFN:
pls_ds <- spls(blocks$diet, blocks$SFN, keepX = keepX, ncomp = 2)
plotVar(pls_ds, cutoff = 0.5, title = 'Diet vs SFN', legend = c('Diet', 'GLSHP'), 
        var.names = F, pch = c(16,17), col = c('darkorchid', 'lightgreen')) #Positive Correlation but not strong
cor(pls_ds$variates$X, pls_ds$variates$Y) #0.68, 0.61

#Connect Diet and SFN to Microbiome but not each other
design <- matrix(c(0,1,1,#M
                   1,0,1,#D
                   1,1,0), ncol = 3, nrow = 3)


#Tune the number of components for Diet vs Microbiome:
dm_perf <- perf(pls_dm, validation = 'loo', auc = T)
dm_perf_mf <- perf(pls_dm, validation = 'Mfold', auc = T, folds = 5, nrepeat = 10)
plot(dm_perf, criterion = 'Q2.total') #1 Component
plot(dm_perf_mf, criterion = 'Q2.total') #1 Component

#Tune the number of components for Microbiome vs SFN:
ms_perf <- perf(pls_ms, validation = 'loo', auc = T)
ms_perf_mf <- perf(pls_ms, validation = 'Mfold', auc = T, folds = 5, nrepeat = 10)
plot(ms_perf, criterion = 'Q2.total') #1 Component
plot(ms_perf_mf, criterion = 'Q2.total') #1 Component

keepX_m <- c(seq(10,100, 5),
             seq(10,100, 5))
keepX_d <- c(seq(10,50, 5),
             seq(10,50, 5))
#Diet vs Microbiome Tuning
tune_dm <- tune.spls(blocks$diet, blocks$micro, test.keepX = keepX_d, test.keepY = keepX_m,  ncomp = 2, mode = 'canonical', 
                    folds =5, nrepeat = 10, measure = 'cor')
plot(tune_dm)
tune_dm$choice.keepX #Diet: 25,30
tune_dm$choice.keepY #Micro: 75,10

#Microbiome vs SFN tuning
tune_ms <- tune.spls(blocks$micro, blocks$SFN, test.keepX = keepX_m, ncomp = 2, mode = 'canonical', 
                    folds =5, nrepeat = 10, measure = 'cor')
plot(tune_ms)
tune_ms$choice.keepX #Micro: 10,10


diet_test <- c(seq(5,50, 5))
micro_test <- c(seq(20,100, 5))
finalstats <- matrix(ncol = 5)
colnames(finalstats) <- c('comp1', 'comp2', 'blocks', 'param_d', 'param_m')
for(dt in diet_test){
  for(mt in micro_test){
    keepX_temp <- list(micro = rep(mt, 2),
                       diet = rep(dt, 2),
                       SFN = rep(3,2))
    tempmod <- block.spls(X = blocks, indY =  3, ncomp = 2, keepX = keepX_temp, mode = 'canonical', design = design)
  outstats <- cbind(rbind(diag(cor(tempmod$variates$micro, tempmod$variates$diet)),
      diag(cor(tempmod$variates$micro, tempmod$variates$SFN))),
      data.frame(blocks = c('MvD', 'MvS'),
                 param_d = rep(dt, 2),
                 param_m = rep(mt, 2)))
  #print(outstats)
  finalstats <- rbind(finalstats, outstats)
  }
}

finalclean <- finalstats %>%
  drop_na() %>%
  pivot_longer(cols = c('comp1', 'comp2'), names_to = 'comp', values_to = 'cor')
outplot <- ggplot(finalclean, aes(x = param_d, y = param_m, size = cor, color = cor)) +
  geom_point() +
  facet_wrap(~blocks+comp) +
  viridis::scale_color_viridis(option = 'C')
plotly::ggplotly(outplot) #MvD wants high but MvS wants low - Split the difference

# Final Model -------------------------------------------------------------
library(mixOmics)

keepX_final <- list(micro = c(55,5),
                    diet = c(25,5),
                    SFN = c(3,3))
design_final <- matrix(c(0,1,1,
                         1,0,0,
                         1,0,0), ncol = 3, nrow = 3)
model_final <- block.spls(X = blocks, indY =  3, ncomp = 2, keepX = keepX_final, mode = 'canonical', design = design_final)
plotVar(model_final, cutoff = 0.5, var.names = F, pch = c(18, 16, 17), col = c('brown', 'darkorchid', 'lightgreen'),
        legend = T)
plotLoadings(model_final, comp = 1, ndisplay = 20, size.name = 1.3)
selectVar(model_final, comp = 1)
cor(model_final$variates$micro, model_final$variates$diet) #0.81, 0.8
cor(model_final$variates$micro, model_final$variates$SFN) #0.76, 0.25

final_ASV <- selectVar(model_final)$micro$value %>%
  filter(abs(value.var) >= 0.1) %>%
  rownames()

final_diet <- selectVar(model_final)$diet$value %>%
  filter(abs(value.var) >= 0.1) %>%
  rownames()

final_loadings_ASVs <- data.frame(tax_table(ps)) %>%
  rownames_to_column('ASV') %>%
  right_join(selectVar(model_final)$micro$value %>% rownames_to_column('ASV')) %>%
  filter(abs(value.var) >= 0.1) 
  
final_loadings_diet <- selectVar(model_final)$diet$value %>%
  filter(abs(value.var) >= 0.1) %>%
  rownames_to_column('Nutrient')

dim(final_loadings_ASVs)
dim(final_loadings_diet)


# Correlation Analysis ----------------------------------------------------
alltax <- data.frame(tax_table(ps_spls)) %>%
  rownames_to_column('ASV') %>%
  mutate(spasv = paste0(Genus, '_', Species, '_', ASV))

ASVtest <- micro_raw %>%
  dplyr::select(all_of(final_ASV)) 

diettest <- diet_raw %>%
  dplyr::select(all_of(final_diet)) 

data.frame(cor(diettest, method = 'sp')) %>%
  pheatmap::pheatmap(display_numbers = T)


#ASV to SFN
data.frame(cor(ASVtest, SFN_raw, method = 'sp')) %>%
  rownames_to_column('ASV') %>% 
  left_join(alltax) %>%
  column_to_rownames('spasv') %>%
  dplyr::select(where(is.numeric)) %>%
  pheatmap::pheatmap(display_numbers = T)

data.frame(cor(ASVtest, SFN_raw, method = 'sp')) %>%
  rownames_to_column('ASV') %>% 
  left_join(alltax) %>%
  column_to_rownames('spasv') %>%
  dplyr::select(where(is.numeric)) %>%
  pheatmap::pheatmap(display_numbers = T)

#ASV to Diet
data.frame(cor(ASVtest, diettest, method = 'sp')) %>%
  rownames_to_column('ASV') %>% 
  left_join(alltax) %>%
  column_to_rownames('spasv') %>%
  dplyr::select(where(is.numeric)) %>%
  pheatmap::pheatmap(display_numbers = T)


# Regression Analysis -----------------------------------------------------

gfinfo <- micro_meta %>%
  left_join(brocinfo) %>%
  dplyr::select(grams_fed)

stepData_gf <- data.frame(cbind(micro_block, SFN_block)) %>%
  dplyr::select(all_of(final_ASV), Tot_SFN) %>%
  cbind(gfinfo)

fullMod_gf <- lm(Tot_SFN ~ ., stepData_gf)
nullMod_gf <- lm(Tot_SFN ~ grams_fed, stepData_gf)
stepFit_gf <- MASS::stepAIC(fullMod_gf, direction = 'backward', trace = T, scope = list(upper = fullMod_gf, lower = nullMod_gf))

finalData <- data.frame(coef(stepFit_gf)) %>%
  rownames_to_column('ASV') %>%
  left_join(alltax)

# Clostridia Only Analysis ------------------------------------------------

cloTest <- micro_raw %>%
  dplyr::select(ASV125, ASV220, ASV276) 

data.frame(cor(cloTest, SFN_raw, method = 'sp')) %>%
  pheatmap::pheatmap(display_numbers = T)

