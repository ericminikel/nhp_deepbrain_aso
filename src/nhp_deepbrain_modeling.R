options(stringsAsFactors=F)
if(interactive()) setwd('~/d/sci/src/nhp_deepbrain_aso/')
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(MASS)); select = dplyr::select; summarize = dplyr::summarize

# statistical parameters
ci_size = .95
alpha = 1 - ci_size # 0.05
qnormval = qnorm(1-alpha/2) # ~1.96

treatments = read_tsv('data/analytical/treatments.tsv', col_types=cols())
animals = read_tsv('data/analytical/animals.tsv', col_types=cols())
regions = read_tsv('data/analytical/region_meta.tsv', col_types=cols())
ana = read_tsv('data/analytical/ana.tsv', col_types=cols())
types = read_tsv('data/analytical/cell_types.tsv', col_types=cols()) %>%
  mutate(y = max(x) - x + 1)

if (exists('afi')) rm(afi)
if (exists('aco')) rm(aco)

for (this_region in regions$region) {
  
  ana %>%
    filter(region==this_region) %>%
    inner_join(types, by='cell_type') %>%
    inner_join(animals, by='animal') %>%
    inner_join(treatments, by='treatment') -> this_ana

  cell_type_smry <- this_ana %>%
    group_by(cell_type, category, x, y) %>%
    summarize(.groups='keep',
              n_cells=n(),
              sum_total_umis = sum(total_umis),
              ctl_cell_ratio = mean(tx =='control'),
              ctl_read_ratio = sum(total_umis[tx =='control'])/sum(total_umis)) %>%
    ungroup() %>%
    rename(n_umis = sum_total_umis) %>%
    arrange(x)
  
  categories <- cell_type_smry %>%
    group_by(category) %>%
    summarize(.groups='keep',
              midx=mean(x), minx=min(x), maxx=max(x),
              midy=mean(y), miny=min(y), maxy=max(y)) 
  
  by_amt <- this_ana %>%
    group_by(h5_num, sample_id, treatment, shortname, cell_type, category, x, y) %>%
    summarize(.groups='keep',
              n_cells = n(),
              total_umis = sum(total_umis),
              malat1_umis = sum(malat1_umis)) %>%
    ungroup() %>%
    filter(n_cells >= 10) %>%
    mutate(subtype_fct = as.factor(cell_type)) %>%
    mutate(malat1_rpm = 1e6 * malat1_umis / total_umis) %>%
    mutate(total_umis_pseudocount = ifelse(total_umis==0, 1, total_umis))
  
  this_model = paste0(this_region,'-Malat1ASO')
    
  fit = glm.nb(malat1_umis ~ subtype_fct + subtype_fct:shortname + offset(log(total_umis_pseudocount)), data=by_amt)
    
  tx_colors = treatments %>% distinct(color, tx, treatment, shortname)
    
  # note that the join solely on subtype here creates a brittleness where subtype
  # has to be unique, so you can't have celltype/subtype combinations such as
  # "neuron/other" and "glia/other", you'll get a many-to-many join here
  
  by_amt %>%
    group_by(cell_type, shortname) %>%
    summarize(.groups='keep', 
              txc_cells = sum(n_cells)) %>%
    ungroup() %>%
    group_by(cell_type) %>%
    mutate(total_cells = sum(txc_cells)) %>%
    ungroup() %>%
    select(shortname, cell_type, n_cells=txc_cells, total_cells) -> cell_counts
  
  as_tibble(summary(fit)$coefficients, rownames='coefname') %>%
    select(coefname, estimate=Estimate, se=`Std. Error`, z=`z value`, p=`Pr(>|z|)`) %>%
    mutate(cell_type = as.factor(gsub('\\(Intercept\\)',levels(by_amt$subtype_fct)[1],gsub('subtype_fct','',gsub(':.*','',coefname))))) %>%
    mutate(group = ifelse(grepl('shortname',coefname),gsub('.*shortname','',coefname),'control')) %>%
    mutate(log_normed = ifelse(group!='control', estimate, 0)) %>%
    mutate(normed = exp(log_normed),
           l95 = exp(log_normed - qnormval*se),
           u95 = exp(log_normed + qnormval*se)) %>%
    inner_join(tx_colors %>% select(-tx), by=c('group'='shortname')) %>%
    inner_join(types[,c('category','cell_type','x','y')], by=c('cell_type')) %>%
    inner_join(cell_counts, by = c("group" = "shortname", "cell_type"))%>%
    rename(shortname=group)-> ba_coefs
  
  by_amt$fit_resid = fit$residuals
  
  by_amt %>%
    select(-x, -y) %>%
    inner_join(select(ba_coefs, -n_cells, -total_cells), by=c('cell_type','shortname','category','treatment')) %>%
    mutate(point_estimate = exp(log_normed + fit_resid)) %>%
    select(-coefname, -estimate, -se, -z, -p, -log_normed, -normed, -l95, -u95) -> ba_fits
  
  ba_coefs$model = this_model
  ba_fits$model = this_model
  
  if(exists('aco')) {
    aco = rbind(aco, ba_coefs)
  } else {
    aco = ba_coefs
  }
  
  if (exists('afi')) {
    afi = rbind(afi, ba_fits)
  } else {
    afi = ba_fits
  }
    
}
  write_tsv(aco, 'data/modeling_outputs/aco.tsv')
  write_tsv(afi, 'data/modeling_outputs/afi.tsv')
  