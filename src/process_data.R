options(stringsAsFactors = F)
if(interactive()) setwd('~/d/sci/src/nhp_deepbrain_aso/')
library(tidyverse)
library(janitor)
library(openxlsx)

name_mapping = read_tsv('data/inputs/name_mapping.tsv', col_types=cols())

region_meta = tibble(region=c('thalamus','caudate','putamen'),
                     order = 1:3)

if(exists('all_ana')) rm(all_ana)
for (each_region in region_meta$region) {
  ana_original = read_tsv(paste0('data/inputs/ana_',each_region,'.tsv'), 
                          col_types=cols()) %>%
    select(-assignments, -prnp_umis, -rnaseh1_umis) %>%
    rename(umap1 = x_coordinate,
           umap2 = y_coordinate)
  amt = read_csv(paste0('data/inputs/amt_',each_region,'.csv'), col_types=cols()) %>%
    clean_names() %>% 
    rename(assigned = 2) %>%
    inner_join(name_mapping, by='assigned') %>%
    select(-assigned)
  ana = inner_join(ana_original, amt, by='barcode') %>%
    mutate(region = each_region)
  if(exists('all_ana')) {
    all_ana = rbind(all_ana, ana)
  } else {
    all_ana = ana
  }
}
all_ana %>% 
  relocate(region) %>%
  mutate(animal = as.integer(substr(sample_id,1,4))) -> all_ana
all_ana %>%
  group_by(cell_type) %>%
  summarize(.groups='keep', 
            n = n(),
            n_regions = length(unique(region)),
            which_regions = toString(unique(region))) %>%
  ungroup() -> count_smry

count_smry

write_tsv(count_smry, 'data/modeling_outputs/count_smry.tsv')
write_tsv(all_ana, 'data/analytical/ana.tsv')

pk = read.xlsx('received/Malat1 NHP PK PD data file.xlsx', cols = 1:5, startRow=2) %>%
  as_tibble() %>%
  rename(region = 1) %>%
  pivot_longer(cols=2:5) %>%
  rename(animal = name, pk = value)
pd = read.xlsx('received/Malat1 NHP qPCR values.xlsx') %>%
  rename(region = 1) %>%
  pivot_longer(cols=2:7) %>%
  rename(animal = name) %>%
  mutate(qpcr = value / 100) %>%
  select(-value)
write_tsv(pk, 'data/jafar-nejad-powers-2021/pk.tsv')
write_tsv(pd, 'data/jafar-nejad-powers-2021/pd.tsv')
