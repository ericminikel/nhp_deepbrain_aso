if(interactive()) setwd('~/d/sci/src/nhp_deepbrain_aso')
library(pzfx)
library(tidyverse)
library(janitor)
library(drc)
library(RColorBrewer)
select = dplyr::select
summarize = dplyr::summarize

source('../helper.R')

tabes = pzfx_tables('received/mouse PK PD Figure1.pzfx')
tabes = tabes[tabes != "Mouse Striatum RNA-PK"] # remove one corrupted table from the PZFX
if (exists('mo')) rm(mo)
for (tabe in tabes) {
  tissue = tolower(gsub('Mouse ','',gsub(' RNA-PK','',tabe)))
  mo_raw = read_pzfx('received/mouse PK PD Figure1.pzfx', table=tabe)
  mo_raw %>%
    rename(dose=ROWTITLE) %>%
    pivot_longer(cols = -dose) %>%
    mutate(tissue = tissue) %>%
    mutate(measurement = gsub('_.+','',name)) %>%
    mutate(animal = gsub('.+_','',name)) %>%
    select(-name) %>%
    pivot_wider(values_from=value, names_from=measurement) %>%
    mutate(RNA = RNA / 100) -> mo_this

  if(exists('mo')) {
    mo = rbind(mo, mo_this)
  } else {
    mo = mo_this
  }
}

write_tsv(mo, 'data/jafar-nejad-powers-2021/mouse_pkpd.tsv', na='')






tabes = pzfx_tables('received/Rat PK PD Figure3.pzfx')
if (exists('ra')) rm(ra)
for (tabe in tabes) {
  tissue = tolower(gsub('Rat ','',gsub(' RNA-PK','',tabe)))
  ra_raw = read_pzfx('received/Rat PK PD Figure3.pzfx', table=tabe)
  ra_raw %>%
    rename(dose=ROWTITLE) %>%
    pivot_longer(cols = -dose) %>%
    mutate(tissue = tissue) %>%
    mutate(measurement = gsub('_.+','',name)) %>%
    mutate(animal = gsub('.+_','',name)) %>%
    select(-name) %>%
    pivot_wider(values_from=value, names_from=measurement) %>%
    mutate(RNA = RNA / 100) -> ra_this
    
    if(exists('ra')) {
      ra = rbind(ra, ra_this)
    } else {
      ra = ra_this
    }
}
write_tsv(ra, 'data/jafar-nejad-powers-2021/rat_pkpd.tsv', na='')


