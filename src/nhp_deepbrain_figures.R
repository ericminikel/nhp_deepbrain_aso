# STARTUP #### 

overall_start_time = Sys.time()
tell_user = function(...) { cat(file=stderr(), paste0(...)); flush.console() }

# DEPENDENCIES ####

tell_user('Loading required packages...')


options(stringsAsFactors=F)
if(interactive()) setwd('~/d/sci/src/nhp_deepbrain_aso')
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(openxlsx))
suppressMessages(library(gridExtra))
suppressMessages(library(ggrepel))
suppressMessages(library(Seurat))
suppressMessages(library(drc))
suppressMessages(library(RColorBrewer))
suppressMessages(library(MASS)); select = dplyr::select; summarize = dplyr::summarize

# OUTPUT STREAMS #### 

tell_user('done.\nCreating output streams...')

text_stats_path = 'display_items/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T
write_stats = function(...) {
  write(paste(list(...),collapse='',sep=''),text_stats_path,append=T)
  write('\n',text_stats_path,append=T)
}

supplement_path = 'display_items/supplement.xlsx'
supplement = createWorkbook()
# options("openxlsx.numFmt" = "0.00") # this looks better for residuals but terrible for p values and weeks post-dose
supplement_directory = tibble(name=character(0), title=character(0))
write_supp_table = function(tbl, title='') {
  # write Excel sheet for supplement
  table_number = length(names(supplement)) + 1
  table_name = paste0('s',formatC(table_number,'d',digits=0,width=2,flag='0'))
  addWorksheet(supplement,table_name)
  bold_style = createStyle(textDecoration = "Bold")
  writeData(supplement,table_name,tbl,headerStyle=bold_style,withFilter=T)
  freezePane(supplement,table_name,firstRow=T)
  saveWorkbook(supplement,supplement_path,overwrite = TRUE)
  
  # also write tab-sep version for GitHub repo
  write_tsv(tbl,paste0('display_items/table-',table_name,'.tsv'), na='')
  
  # and save the title in the directory tibble for later
  assign('supplement_directory',
         supplement_directory %>% add_row(name=table_name, title=title),
         envir = .GlobalEnv)
}


# FUNCTIONS & CONSTANTS #### 

tell_user('done.\nDefining functions & constants...')

ctl_color = '#979797'
mal_color = '#FF00FF'

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x <= 0, '', ifelse(signed, '+', '')),formatC(100*x,format='f',digits=digits),'%'))

upper = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) + sds*sd(x)/sqrt(sum(!is.na(x)))
}
lower = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) - sds*sd(x)/sqrt(sum(!is.na(x)))
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

clipdist = function(x, minx, maxx) {
  return (pmin(maxx,pmax(minx,x)))
}

rbind_files = function(path, grepstring) {
  all_files = list.files(path, full.names=T)
  these_files = all_files[grepl(grepstring,all_files)]
  for (this_file in these_files) {
    this_tbl = read_delim(this_file, col_types=cols()) %>% clean_names()
    this_tbl$file = gsub('.*\\/','',gsub('\\.[tc]sv','',this_file))
    if (exists('rbound_table')) {
      rbound_table = rbind(rbound_table, this_tbl)
    } else {
      rbound_table = this_tbl
    }
  }
  return (rbound_table)
}

format_p = function(p) {
  formatted_p = paste0(' = ',formatC(signif(p,2), format='g', digits=2))
  formatted_p[signif(p,2) < 0.1] = paste0(' = ',formatC(signif(p[signif(p,2) < 0.1],2), format='f', digits=3))
  formatted_p[signif(p,2) < 0.01] = paste0(' = ',formatC(signif(p[signif(p,2) < 0.1],2), format='e', digits=1))
  formatted_p[signif(p,2) == 1] = ' = 1.0'
  formatted_p[signif(p,2) < 1e-15] = ' < 1.0e-15'
  return (formatted_p)
}

# DATA #### 

tell_user('done.\nReading in data...')

treatments = read_tsv('data/analytical/treatments.tsv', col_types=cols())
animals = read_tsv('data/analytical/animals.tsv', col_types=cols())
regions = read_tsv('data/analytical/region_meta.tsv', col_types=cols())
ana = read_tsv('data/analytical/ana.tsv', col_types=cols())
types = read_tsv('data/analytical/cell_types.tsv', col_types=cols()) %>%
  mutate(y = max(x) - x + 1)

metrics = rbind_files('data/metrics_summary/','.csv') %>% 
  separate_wider_delim(file, '_', names=c('animal','region','drop1','drop2')) %>%
  select(-drop1, -drop2) %>%
  relocate(animal, region)

write_supp_table(metrics, 'Single nucleus sequencing quality control metrics.')

# FIGURE 1 #### 

tell_user('done.\nCreating Figure 1...')

resx=600
png(paste0('display_items/figure-1.png'),width=6.5*resx,height=4.5*resx,res=resx)

layout_matrix = matrix(c(1,2,3,4), nrow=1, byrow=T)
layout(layout_matrix, widths=c(1,.3,.6,.5))


panel = 1

pk_raw = read_tsv('data/jafar-nejad-powers-2021/nhp_pk.tsv', col_types=cols())
region_names = read_tsv('data/jafar-nejad-powers-2021/region_name_mapping.tsv', col_types=cols())
pk_raw %>%
  inner_join(region_names, by=c('region'='pk_name')) %>%
  inner_join(animals, by='animal') %>%
  group_by(meta_region, animal, pch) %>%
  summarize(.groups='keep', pk = mean(pk)) %>%
  ungroup() -> pk

#setdiff(pk_raw$region, region_names$pk_name)
#setdiff(region_names$pk_name, pk_raw$region)


pk %>%
  group_by(meta_region) %>%
  summarize(.groups='keep',
            n=n(),
            mean = mean(pk),
            l95 = lower(pk),
            u95 = upper(pk)) %>%
  ungroup() %>%
  arrange(desc(mean)) %>%
  mutate(x = rank(-mean), y=rank(mean)) %>%
  mutate(font = case_when(meta_region %in% c('caudate','putamen','thalamus') ~ 2,
                          TRUE ~ 1)) -> pk_smry





pd_raw = read_tsv('data/jafar-nejad-powers-2021/nhp_pd.tsv', col_types=cols())
region_names = read_tsv('data/jafar-nejad-powers-2021/region_name_mapping.tsv', col_types=cols())
pd_raw %>%
  inner_join(animals, by='animal') %>%
  filter(treatment=='Malat1 ASO') %>%
  inner_join(region_names, by=c('region'='qpcr_name')) %>%
  group_by(meta_region, animal, pch) %>%
  summarize(.groups='keep', qpcr=mean(qpcr)) %>%
  ungroup() -> pd

pd %>%
  group_by(meta_region) %>%
  summarize(.groups='keep',
            mean = mean(qpcr),
            l95 = lower(qpcr),
            u95 = upper(qpcr)) %>%
  ungroup() %>%
  arrange(desc(mean)) %>%
  mutate(x = rank(mean), y=rank(-mean)) %>%
  mutate(font = case_when(meta_region %in% c('caudate','putamen','thalamus') ~ 2,
                          TRUE ~ 1)) -> pd_smry


pd_smry %>%
  inner_join(pk_smry, by='meta_region', suffix=c('_pd','_pk')) %>%
  mutate(y = rank(y_pd)) -> pkpd

pd$y = pkpd$y[match(pd$meta_region, pkpd$meta_region)]
pk$y = pkpd$y[match(pk$meta_region, pkpd$meta_region)]

par(mar=c(3,4,3,1))
xlims = c(0,50)
xbigs = 0:5*10
xats = 0:50
ybigs = 0:2/2
yats = 0:10/10
ylims = c(0,1.05)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xbigs, tck=-0.04, labels=NA)
axis(side=1, at=xbigs, lwd=0, labels= xbigs, line=-0.35)
axis(side=1, at=xats, tck=-0.02, labels=NA)
mtext(side=1, line=1.8, text='tissue ASO concentration µg/g', cex=0.9)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, lwd=0, labels=percent(ybigs), las=2, line=-0.25)
axis(side=2, at=yats, tck=-0.025, labels=NA)
mtext(side=2, line=2.4, text=expression(paste(italic('Malat1'),' (% residual)')), cex=0.9)
points(pkpd$mean_pk, pkpd$mean_pd, pch=20)
par(xpd=T)
text(pkpd$mean_pk, pkpd$mean_pd, labels=gsub('temporal cortex','\ntemporal\ncortex',pkpd$meta_region), font=pkpd$font_pd, cex=0.6,
      pos = case_when(pkpd$mean_pk < 12.5 & pkpd$mean_pd < 0.25 ~ 2, 
                      pkpd$mean_pd < 0.09 ~ 1,
                      TRUE ~ 4))
par(xpd=F)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.4)
panel = panel + 1

par(mar=c(3,0,3,0))
ylims = range(pkpd$y) + c(-0.5, 0.5)
xlims = 0:1
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
mtext(side=4, line=-0.1, at=pkpd$y, text=pkpd$meta_region, las=2, adj=1, cex=0.6,
      font = case_when(pkpd$meta_region %in% regions$region ~ 2, TRUE ~ 1))


par(mar=c(3,0,3,1))
ylims = range(pkpd$y) + c(-0.5, 0.5)
xlims = c(0, 1.35)
xats = 0:20/10
xbigs = 0:4/2
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
abline(v=0, lwd=2)
abline(v=1, lwd=0.75, lty=3)
axis(side=1, at=xats, tck=-0.02, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, labels=percent(xbigs), lwd=0, line=-0.5, cex.axis=0.8)
mtext(side=1, line=2.0, text=expression(paste(italic('Malat1'),' (% residual)')), cex=0.9)
barwidth = 0.4
rect(xleft=rep(0, nrow(pkpd)), xright=pkpd$mean_pd, ybottom=pkpd$y-barwidth, ytop=pkpd$y+barwidth, col=alpha(mal_color,ci_alpha), border=NA)
set.seed(1)
points(x=pd$qpcr, y=pd$y, pch=pd$pch, col=mal_color, bg='#FFFFFF')
for (anml in animals$animal) {
  subs = pd %>% subset(animal == anml) %>% arrange(y)
  points(x=subs$qpcr, y=subs$y, col=mal_color, type='l', lwd=0.5)
}
leg = animals %>% filter(treatment=='Malat1 ASO')
legend('topright', as.character(leg$animal), title='animal', pch=leg$pch, col=mal_color, bty='n', cex=0.8)
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.4)
panel = panel + 1


par(mar=c(3,0,3,1))
ylims = range(pkpd$y) + c(-0.5, 0.5)
xlims = c(0.07, max(pk$pk)) * 1.2
xats = rep(1:9, 3) * 10^rep(-1:1,each=9)
xbigs = 10^(-1:1)
plot(NA, NA, xlim=xlims, ylim=ylims, log='x', axes=F, ann=F, xaxs='i', yaxs='i')
abline(v=0.1)
axis(side=1, at=xats, tck=-0.02, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, lwd=0, line=-0.5)
mtext(side=1, at=0.1, line=1.6, text='LLQ',cex=0.8)
mtext(side=1, line=1.8, text='PK (µg/g)', cex=0.9)
barwidth = 0.4
rect(xleft=rep(0.1, nrow(pkpd)), xright=pkpd$mean_pk, ybottom=pkpd$y-barwidth, ytop=pkpd$y+barwidth, col=alpha(mal_color,ci_alpha), border=NA)
set.seed(1)
points(x=pk$pk, y=pk$y, pch=pk$pch, col=mal_color, bg='#FFFFFF')
for (anml in animals$animal) {
  subs = pk %>% subset(animal == anml) %>% arrange(y)
  points(x=subs$pk, y=subs$y, col=mal_color, type='l', lwd=0.5)
}
mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.4)
panel = panel + 1

silence = dev.off()





# FIGURE S1 ####

tell_user('done.\rCreating Figure S1...')
resx=600
png('display_items/figure-s1.png',width=6.5*resx,height=8*resx,res=resx)
par(mfrow=c(2,1), mar=c(3,3,3,12))


panel = 1

mo = read_tsv('data/jafar-nejad-powers-2021/mouse_pkpd.tsv', col_types=cols())
mo %>%
  distinct(tissue) -> mo_tmeta
mo_tmeta$color = brewer.pal(n = nrow(mo_tmeta), name = 'Dark2')

mo$color = mo_tmeta$color[match(mo$tissue, mo_tmeta$tissue)]


xlims = c(0.001, 100)
xats = rep(1:9, 6) * rep(10^(-3:2), each=9)
xbigs = 10^(-3:2)
xbiglabs = c('1 ng/g', '10 ng/g', '100 ng/g', '1 µg/g', '10 µg/g', '100 µg/g')
ylims = c(0, 1.5)
yats = 0:6/4
ybigs = 0:3/2
ybiglabs = percent(0:3/2)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, log='x', xaxs='i', yaxs='i')
axis(side=1, at=xats, tck=-0.02, labels=NA)
axis(side=1, at=xbigs, tck=-0.04, labels=NA)
axis(side=1, at=xbigs, labels = xbiglabs, lwd=0, line=-0.5, cex.axis=0.6)
mtext(side=1, line =1.6, text='drug tissue concentration (µg/g)')
axis(side=2, at=yats, tck=-0.02, labels=NA)
axis(side=2, at=ybigs, tck=-0.04, labels=NA)
axis(side=2, at=ybigs, labels = ybiglabs, las=2, lwd=0, line=-0.5, cex.axis=0.6)
mtext(side=2, line =2.2, text='residual Malat1 (%)')
abline(h=c(0.5,1.0),  lty=2, lwd=0.5)
points(mo$PK, mo$RNA, col=mo$color, pch=20)
mo_tmeta$ic50 = as.numeric(NA)
for (i in 1:nrow(mo_tmeta)) {
  this_tissue = mo_tmeta$tissue[i]
  model = drm(RNA ~ PK, data = subset(mo, tissue==this_tissue), fct=LL.4(fixed=c(b=NA,c=0,d=1,e=NA)))
  x = xats 
  y = predict(model, newdata=data.frame(PK=x))
  points(x, y, type='l', lwd=1, col=mo_tmeta$color[i])
  mo_tmeta$ic50[i] = coefficients(model)['e:(Intercept)']
}

mo_tmeta %>%
  arrange(desc(ic50)) -> mo_tmeta
mo_tmeta$disp = as.character(NA)
for (i in 1:nrow(mo_tmeta)) {
  mo_tmeta$disp[i] = paste0(mo_tmeta$tissue[i], paste0(rep(' ',max(nchar(mo_tmeta$tissue))-nchar(mo_tmeta$tissue[i])+1),collapse=''), formatC(mo_tmeta$ic50[i], format='f', digits=2), ' µg/g')
}

par(xpd=T)
par(family='mono')
legend(x=max(xlims),y=max(ylims), mo_tmeta$disp, col=mo_tmeta$color, pch=20, lwd=1, bty='n', cex=0.8)
par(family='sans')
par(xpd=F)

mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.4)
mtext('mouse', side=3, line = 1)
panel = panel + 1

write_supp_table(mo, 'Mouse PK/PD data from Jafar-Nejad & Powers 2021 Supplementary Table 3.')
write_supp_table(mo_tmeta, 'Mouse IC50 values by tissue.')


ra = read_tsv('data/jafar-nejad-powers-2021/rat_pkpd.tsv', col_types=cols())

ra %>%
  distinct(tissue) %>%
  filter(tissue != 'striatum') -> ra_tmeta # remove striatum, too many missing values, model does not converge
ra_tmeta$color = c(brewer.pal(n = 6, name = 'Dark2'), brewer.pal(n = nrow(ra_tmeta)-6, name = 'Accent'))

ra$color = ra_tmeta$color[match(ra$tissue, ra_tmeta$tissue)]

plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, log='x', xaxs='i', yaxs='i')
axis(side=1, at=xats, tck=-0.02, labels=NA)
axis(side=1, at=xbigs, tck=-0.04, labels=NA)
axis(side=1, at=xbigs, labels = xbiglabs, lwd=0, line=-0.5, cex.axis=0.6)
mtext(side=1, line =1.6, text='drug tissue concentration (µg/g)')
axis(side=2, at=yats, tck=-0.02, labels=NA)
axis(side=2, at=ybigs, tck=-0.04, labels=NA)
axis(side=2, at=ybigs, labels = ybiglabs, las=2, lwd=0, line=-0.5, cex.axis=0.6)
mtext(side=2, line =2.2, text='residual Malat1 (%)')
abline(h=c(0.5,1.0),  lty=2, lwd=0.5)
points(ra$PK, ra$RNA, col=ra$color, pch=20)
ra_tmeta$ic50 = as.numeric(NA)
for (i in 1:nrow(ra_tmeta)) {
  this_tissue = ra_tmeta$tissue[i]
  model = drm(RNA ~ PK, data = subset(ra, tissue==this_tissue), fct=LL.4(fixed=c(b=NA,c=0,d=1,e=NA)))
  x = xats 
  y = predict(model, newdata=data.frame(PK=x))
  points(x, y, type='l', lwd=1, col=ra_tmeta$color[i])
  ra_tmeta$ic50[i] = coefficients(model)['e:(Intercept)']
}

ra_tmeta %>%
  arrange(desc(ic50)) -> ra_tmeta
ra_tmeta$disp = as.character(NA)
for (i in 1:nrow(ra_tmeta)) {
  ra_tmeta$disp[i] = paste0(ra_tmeta$tissue[i], paste0(rep(' ',max(nchar(ra_tmeta$tissue))-nchar(ra_tmeta$tissue[i])+1),collapse=''), formatC(ra_tmeta$ic50[i], format='f', digits=2), ' µg/g')
}



par(xpd=T)
par(family='mono')
legend(x=max(xlims),y=max(ylims), ra_tmeta$disp, col=ra_tmeta$color, pch=20, lwd=1, bty='n', cex=0.8)
par(family='sans')
par(xpd=F)

mtext(LETTERS[panel], side=3, cex=1.5, adj = 0.0, line = 0.4)
mtext('rat', side=3, line = 1)
panel = panel + 1

write_supp_table(ra, 'Rat PK/PD data from Jafar-Nejad & Powers 2021 Supplementary Table 4.')
write_supp_table(ra_tmeta, 'Rat IC50 values by tissue.')


silence = dev.off()








# FIGURE S2 #### 

tell_user('done.\nCreating Figure S2...')


png(paste0('display_items/figure-s2.png'),width=6.5*resx,height=4.5*resx,res=resx)

pd %>%
  group_by(meta_region) %>%
  mutate(pd_diff = qpcr - mean(qpcr)) %>%
  ungroup() -> pd_offset
pk %>%
  group_by(meta_region) %>%
  mutate(pk_ratio = pk / mean(pk)) %>%
  ungroup() -> pk_offset
pd_offset %>%
  inner_join(pk_offset, by=c('animal', 'meta_region', 'pch')) %>%
  select(animal, meta_region, pd_diff, pk_ratio, pch) -> pd_pk_offset

xlims = c(0.02, 4)
ylims = c(-.5, 0.8)
ybigs = c(-.5, 0, .5, 1)
yats = -5:10/10
xats = rep(1:9, 4) * 10^rep(-2:1,each=9)
xbigs = 10^(-2:1)
plot(NA, NA, xlim=xlims, ylim=ylims, log='x', ann=F, axes=F, xaxs='i', yaxs='i')
axis(side=1, at=xats, tck=-0.02, labels=NA)
axis(side=1, at=xbigs, tck=-0.05, labels=NA)
axis(side=1, at=xbigs, lwd=0, line=-0.5)
mtext(side=1, line=1.6, text='PK ratio')
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, lwd=0, labels=percent(ybigs), las=2, line=-0.25)
axis(side=2, at=yats, tck=-0.025, labels=NA)
mtext(side=2, line=2.4, text='PD difference')
abline(h=0)
abline(v=1)
points(pd_pk_offset$pk_ratio, pd_pk_offset$pd_diff, pch=pd_pk_offset$pch, col='#00000088')
m = lm(pd_diff ~ logpkratio, data=pd_pk_offset %>% mutate(logpkratio = log(pk_ratio)))
pval = summary(m)$coefficients['logpkratio','Pr(>|t|)']
write_stats('Linear model: PD difference ~ log(PK) ratio, P = ',format_p(pval))
x = log(seq(0.01, 3.00, 0.01))
#x = seq(0.01, 3.00, 0.01)
y = predict(m, newdata=data.frame(logpkratio=x))
points(exp(x), y, type='l', col='red')
leg = animals %>% filter(treatment=='Malat1 ASO')
legend('topright', as.character(leg$animal), title='animal', pch=leg$pch, bty='n', cex=0.8)

silence = dev.off()





# FIGURE 2, S3 ####

if (file.exists('data/h5/')) {
  tell_user('done.\nCreating Figure 2 & S3...')
  
  
  panels = {}
  panel = 1
  supp_panels = {}
  supp_panel = 1
  
  for (this_region in regions$region) {
    
    disp_panel = (panel + 1) %/% 2
    
    single_label = c('oligodendrocyte')
    
    read_tsv('data/analytical/ana.tsv', col_types=cols()) %>%
      filter(region==this_region) %>%
      rename(umap_1 = umap1, umap_2 = umap2) %>%
      inner_join(types, by='cell_type') %>%
      select(cell_type, umap_1, umap_2, color) -> umapr
    umapr %>%
      filter(abs(umap_1) > 3 | abs(umap_2) > 3) %>%
      filter(!(cell_type %in% single_label)) %>%
      mutate(quadrant = paste0(umap_1 > 0, '-', umap_2 > 0)) %>%
      group_by(cell_type, quadrant) %>%
      summarize(.groups='keep',
                umap1_mean=mean(umap_1),
                umap2_mean=mean(umap_2),
                n_cells = n()) %>%
      ungroup() %>%
      filter(n_cells > 10) %>%
      filter(!(cell_type %in% single_label)) -> umap_centroids1
    umap_centroids1 %>%
      group_by(cell_type) %>%
      summarize(.groups='keep',
                n_quads = n(),
                dist=suppressWarnings(sqrt(max(diff(umap1_mean))^2 + max(diff(umap2_mean))^2))) -> quadcount
    umap_centroids1 %>%
      filter(cell_type %in% quadcount$cell_type[quadcount$n_quads > 1 & quadcount$dist > 5]) -> umap_centroids1
    
    umapr %>%
      filter(!(cell_type %in% quadcount$cell_type[quadcount$n_quads > 1 & quadcount$dist > 5])) %>%
      group_by(cell_type) %>%
      summarize(.groups='keep', umap1_mean=mean(umap_1), umap2_mean=mean(umap_2), n_cells=n()) %>%
      ungroup() %>% mutate(quadrant='') -> umap_centroids2
    rbind(umap_centroids1, umap_centroids2) -> umap_centroids
    
    
    this_umap_plot = ggplot(umapr, aes(umap_1, umap_2, colour=color)) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")) +
      geom_point(size=0.5, shape=20) +
      scale_color_identity() +
      labs(x = 'UMAP 1', y='UMAP 2') +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank()) +
      geom_text_repel(data=umap_centroids, point.size = NA, size=2, min.segment.length = 0,
                      aes(umap1_mean, umap2_mean, label=cell_type, color='black')) +
      ggtitle(this_region) +
      labs(tag = LETTERS[disp_panel])
    panels[[panel]] = this_umap_plot
    panel = panel + 1
    
    
    h5 = Read10X_h5(paste0('data/h5/',this_region,'_filtered_feature_bc_matrix.h5'),use.names=T)
    types = read_tsv('data/analytical/cell_types.tsv',col_types=cols())
    supplementary_gene_list = read_tsv('data/analytical/cell_types_with_supp_genes.tsv',col_types=cols()) %>%
      clean_names() %>%
      distinct(best_marker)
    ana = read_tsv('data/analytical/ana.tsv',col_types=cols()) %>%
      inner_join(types, by='cell_type') %>%
      filter(region==this_region)
    
    unique_levels = rev(types$cell_type[types$cell_type %in% ana$cell_type])
    ana$cell_type_fct = factor(ana$cell_type, levels = unique_levels, ordered = TRUE)
    ana %>%
      select(barcode, cell_type_fct) -> ana_idents
    
    seur = CreateSeuratObject(h5)
    
    Idents(seur) %>%
      as_tibble(rownames="barcode") %>%
      left_join(ana_idents, by="barcode") -> identities
    
    Idents(seur) = identities$cell_type_fct
    seur = subset(seur, idents = unique(ana$cell_type))
    
    this_dotplot = DotPlot(seur, features=types$best_marker[types$cell_type %in% ana$cell_type], dot.scale=4) +
      xlab('') + ylab('') +
      theme(axis.ticks = element_blank()) +
      theme(axis.text.y = element_text(size=8)) +
      theme(plot.margin = margin(l = 1, t=3, b=1, r=1)) +
      theme(axis.text.x = element_text(angle = 45, hjust=1, size=7, face='italic')) +
      theme(legend.position="none")
    
    supp_dotplot = DotPlot(seur, features=supplementary_gene_list$best_marker, dot.scale=3) +
      xlab('') + ylab('') +
      theme(axis.ticks = element_blank()) +
      theme(axis.text.y = element_text(size=8)) +
      theme(plot.margin = margin(l = 1, t=3, b=1, r=1)) +
      theme(axis.text.x = element_text(angle = 45, hjust=1, size=7, face='italic')) +
      theme(legend.position="none") +
      labs(tag = LETTERS[supp_panel]) +
      ggtitle(this_region)
    
    supp_panels[[supp_panel]] = supp_dotplot
    supp_panel = supp_panel + 1
    
    panels[[panel]] = this_dotplot
    panel = panel + 1
  }
  png(paste0('display_items/figure-2.png'),width=6.5*resx,height=8*resx,res=resx)
  grid.arrange(grobs=panels, ncol=2, widths=c(1,1))
  silence = dev.off()
  
  png(paste0('display_items/figure-s3.png'),width=6.5*resx,height=8*resx,res=resx)
  grid.arrange(grobs=supp_panels, ncol=1, widths=c(1))
  silence = dev.off()
} else {
  tell_user('done.\nSkipping Figure 2 & S3 (h5 files too large for GitHub)...')
}






# FIGURE 3 #### 

tell_user('done.\nCreating Figure 3...')


resx=300
png('display_items/figure-3.png',width=6.5*resx,height=7*resx,res=resx)

layout_matrix = matrix(c(1,1,1,1,
                         1,1,1,1,
                         2,2,2,2,
                         2,2,2,2,
                         3,3,3,3,
                         3,3,3,3), nrow=6, byrow=T)
layout(layout_matrix, 
       heights=c(1,1,1.1,1.1,1,1))

panel = 1
#layout.show(3)

ana = read_tsv('data/analytical/ana.tsv', col_types=cols())
types = read_tsv('data/analytical/cell_types.tsv',col_types=cols())

pd_raw %>%
  inner_join(animals, by='animal') %>%
  inner_join(region_names, by=c('region'='qpcr_name')) %>%
  group_by(meta_region, animal, pch) %>%
  summarize(.groups='keep', qpcr=mean(qpcr)) %>%
  ungroup() -> pd_all

bulk_sc <- ana %>%
  inner_join(animals, by='animal') %>%
  inner_join(treatments, by='treatment') %>%
  inner_join(types, by='cell_type') %>%
  group_by(animal, region, tx, treatment, pch) %>%
  summarize(.groups='keep',
            malat1_rpm = 1e6*mean(malat1_umis/total_umis)) %>%
  ungroup() %>%
  group_by(region) %>%
  mutate(malat1_rel = malat1_rpm / mean(malat1_rpm[tx =='control'])) %>%
  ungroup() %>%
  inner_join(pd_all %>% select(-pch), by=c('animal',c('region'='meta_region')))

bulk_sc %>%
  inner_join(treatments %>% select(-tx), by='treatment') %>%
  select(animal, region, treatment, color, pch, malat1_rel, qpcr) %>%
  pivot_longer(cols=c(malat1_rel, qpcr)) %>%
  mutate(ylab = case_when(name=='malat1_rel' ~ 'snRNAseq all',
                          name=='qpcr' ~ 'bulk qPCR')) -> bulk_sc_long


for (this_region in regions$region) {
  
  # tell_user('\rCreating Figure 3... in progress: ',this_region)
  
  this_model = paste0(this_region,'-Malat1ASO')
  
  aco = read_tsv('data/modeling_outputs/aco.tsv', col_types=cols()) %>%
    filter(model == this_model) %>%
    filter(shortname=='malat1aso')  %>% # only plot the fits for active
    mutate(y = rank(y))
  
  aco %>% distinct(cell_type, y) -> newy
  
  afi = read_tsv('data/modeling_outputs/afi.tsv', col_types=cols()) %>%
    filter(model == this_model) %>%
    select(-y) %>%
    inner_join(newy, by='cell_type') %>%
    mutate(animal = as.integer(substr(sample_id, 1, 4)),
           region = gsub('.*_','',sample_id)) %>%
    inner_join(animals %>% select(-treatment), by='animal')
  
  aco %>%
    distinct(y, category, cell_type) %>%
    arrange(y) -> yleg
  
  yleg %>%
    group_by(category) %>%
    summarize(.groups='keep',
              miny=min(y),
              maxy=max(y),
              midy = (min(y)+max(y))/2) -> tranches
  
  
  bulk_sc_long %>%
    filter(region==this_region) %>%
    mutate(y = case_when(ylab=='bulk qPCR' ~ max(newy$y) + 2,
                         ylab=='snRNAseq all' ~ max(newy$y) + 1)) -> bulk_sc_this
  
  bulk_sc_this %>%
    filter(treatment=='Malat1 ASO') %>%
    group_by(region, y, ylab, color) %>%
    summarize(.groups='keep',
              mean = mean(value),
              l95 = lower(value),
              u95 = upper(value)) %>%
    ungroup() -> bulk_sc_smry
  
  bulk_sc_this %>%
    distinct(y, ylab) -> moreleg
  
  if (panel==1) {
    par(mar=c(0,11,3,1))
  } else if (panel==2) {
    par(mar=c(1.5,11,1.5,1))
  } else if (panel==3) {
    par(mar=c(4,11,0,1))
  }
  
  xlims = c(0, 1.5)
  ylims = range(c(afi$y,bulk_sc_this$y)) + c(-0.5, 0.5)
  #ylims = range(afi$y) + c(-0.5, 0.5)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=0:6/4, labels=NA, tck=-0.04)
  abline(v=1, col=ctl_color, lwd=1, lty=3)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  mtext(side=2, at=yleg$y, text=yleg$cell_type, line=0.25, las=2, cex=0.75)
  mtext(side=2, at=moreleg$y, text=moreleg$ylab, line=0.25, las=2, cex=0.75)
  tranche_line = 8
  overhang = 0.3
  for (i in 1:nrow(tranches)) {
    axis(side=2, line=tranche_line, at=c(tranches$miny[i], tranches$maxy[i]) + c(-1,1)*overhang, tck=0.05, labels=NA)
    mtext(side=2, line=tranche_line + 0.25, at=tranches$midy[i], text=tranches$category[i], cex=0.8)
  }
  points(y=afi$y, x=afi$point_estimate, pch=afi$pch, cex=1.5, col=alpha(afi$color, ci_alpha))
  points(y=bulk_sc_this$y, x=bulk_sc_this$value, pch=bulk_sc_this$pch, cex=1.5, col=alpha(bulk_sc_this$color, ci_alpha))
  for (this_animal in unique(animals$animal[animals$treatment=='Malat1 ASO'])) {
    subs = afi %>% 
      filter(animal == this_animal) %>%
      arrange(desc(y))
    points(y=subs$y, x=subs$point_estimate, type='l', lwd=0.25, col=subs$color)
  }
  barwidth=0.4
  segments(y0=aco$y-barwidth, y1=aco$y+barwidth, x0=aco$normed, col=aco$color, lwd=2, lend=1)
  arrows(y0=aco$y, x0=aco$l95, x1=aco$u95, angle=90, length=0.03, code=3, col=aco$color, lwd=1.5)
  mtext(side=2, line=-0.2, adj=0, at=aco$y, text=percent(aco$normed), cex=0.5, las=2)
  segments(y0=bulk_sc_smry$y-barwidth, y1=bulk_sc_smry$y+barwidth, x0=bulk_sc_smry$mean, col=bulk_sc_smry$color, lwd=2, lend=1)
  arrows(y0=bulk_sc_smry$y, x0=bulk_sc_smry$l95, x1=bulk_sc_smry$u95, angle=90, length=0.03, code=3, col=bulk_sc_smry$color, lwd=1.5)
  mtext(side=2, line=-0.2, adj=0, at=bulk_sc_smry$y, text=percent(bulk_sc_smry$mean), cex=0.5, las=2)
  
  mtext(side=2, line=tranche_line+1.25, text=this_region, font=1, cex=1)
  
  if (panel == 1) {
    legend('topright', as.character(leg$animal), title='animal', pch=leg$pch, col=mal_color, bty='n', cex=0.8)
  }
  
  if (panel == 3) {
    axis(side=1, at=0:6/4, labels=percent(0:6/4), lwd=0, line=-0.2, cex.axis=1)
    mtext(side=1, line=2.2, text=substitute(paste('residual ',italic('Malat1'))), cex=1)
  }
  panel = panel + 1
}

silence = dev.off()



# FIGURE S5 ####

tell_user('done.\rCreating Figure S5...')

png('display_items/figure-s5.png',width=6.5*resx,height=8*resx,res=resx)
par(mfrow=c(8,4), mar=c(3,3,1,1))

aco = read_tsv('data/modeling_outputs/aco.tsv', col_types=cols())
aco %>%
  mutate(region = gsub('-.*','',model)) %>%
  distinct(region, cell_type) -> valid_region_types

ana = read_tsv('data/analytical/ana.tsv', col_types=cols())
ana %>%
  inner_join(animals, by='animal') -> ana_with

for (i in 1:nrow(valid_region_types)) {
  this_region = valid_region_types$region[i]
  this_cell_type = valid_region_types$cell_type[i]
  this_subtype_data = ana_with %>% filter(region==this_region & cell_type==this_cell_type)
  binsize = max(10^floor(log10(max(this_subtype_data$malat1_umis))-2),10)
  control_hist = hist(this_subtype_data$malat1_umis[this_subtype_data$treatment=='aCSF'], breaks=seq(0,1e5,binsize), plot=F)
  treatment_hist = hist(this_subtype_data$malat1_umis[this_subtype_data$treatment=='Malat1 ASO'], breaks=seq(0,1e5,binsize), plot=F)
  xmax = as.numeric(ceiling(quantile(this_subtype_data$malat1_umis[this_subtype_data$treatment=='aCSF'], .9)/binsize)*binsize)
  ymax = max(treatment_hist$counts)
  xlims = c(0,xmax*1.2)
  xbigs = seq(0,xmax,binsize*5)
  xats = seq(0,xmax,binsize)
  ylims = c(0,ymax)
  plot(NA, NA, xlim=xlims, ylim = ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  axis(side=2, at=ylims, labels=NA, lwd.ticks=0)
  rect(xleft=control_hist$mids-binsize/2, xright=control_hist$mids+binsize/2, ybottom=rep(0,length(control_hist)), ytop=control_hist$counts, col=alpha(ctl_color, ci_alpha), border=NA)
  rect(xleft=treatment_hist$mids-binsize/2, xright=treatment_hist$mids+binsize/2, ybottom=rep(0,length(treatment_hist)), ytop=treatment_hist$counts, col=alpha(mal_color, ci_alpha), border=NA)
  residual_fit = aco$normed[aco$model == paste0(this_region,'-Malat1ASO') & aco$cell_type==this_cell_type & aco$treatment=='Malat1 ASO']
  axis(side=1, labels=NA, tck=-0.025)
  axis(side=1, line=-1, lwd=0, cex.axis=0.8)
  mtext(side=1, line=1, text='Malat 1 UMIs', cex=0.6)
  axis(side=2, labels=NA, tck=-0.025)
  axis(side=2, line=-.75, lwd=0, las=2, cex.axis=0.8)
  mtext(side=2, line=1.25, text='N nuclei', cex=0.6)
  mtext(side=3, line=0, text=paste0(' ',this_region,' ',this_cell_type, ' ', percent(round(residual_fit*100)/100)), cex=0.6)
  axis(side=4, at=ylims, lwd.ticks=0, labels=NA)
  
}
silence = dev.off()

# FIGURE S4 #### 

tell_user('done.\nCreating Figure S4...')


resx=600
png('display_items/figure-s4.png',width=6.5*resx,height=7*resx,res=resx)

layout_matrix = matrix(1:21, nrow=7, byrow=F)
layout(layout_matrix, heights=c(0.3, rep(1,6)))


ana = read_tsv('data/analytical/ana.tsv', col_types=cols())
for (this_region in regions$region) {
  
  xlims = range(ana$umap1[ana$region==this_region])*1.2
  ylims = range(ana$umap2[ana$region==this_region])*1.2
  
  par(mar=c(0,1,0,1))
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i') 
  mtext(side=1,line=-1, text=this_region, cex=1.5)
  for (this_animal in animals$animal) {
    ana %>%
      inner_join(types %>% select(cell_type), by='cell_type') %>%
      filter(animal==this_animal) %>%
      filter(region==this_region) %>%
      inner_join(animals, by='animal') %>%
      inner_join(treatments, by='treatment') %>%
      select(animal, treatment, umap1, umap2, cell_type, color) -> this_umap
    
    par(mar=c(1,1,2,1))
    plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i') 
    axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
    mtext(side=1, line=0, text='UMAP 1', cex=0.5)
    axis(side=2, at=ylims, lwd.ticks=0, labels=NA)
    mtext(side=2, line=0, text='UMAP 2', cex=0.5)
    points(this_umap$umap1, this_umap$umap2, pch=20, cex=0.15, col=this_umap$color)
    mtext(side=3, text=this_animal, cex=0.8)
  }
}
silence = dev.off()








# OTHER #### 

tell_user('done.\nPreparing supplementary tables...')

pk %>%
  inner_join(animals, by=c('animal','pch')) %>%
  select(-pch, -y) -> pk_out


pd %>%
  inner_join(animals, by=c('animal','pch')) %>%
  select(-pch, -y) -> pd_out

pkpd %>%
  select(meta_region, n, mean_pd, l95_pd, u95_pd, mean_pk, l95_pk, u95_pk) -> pkpd_out

write_supp_table(pkpd_out, 'Pharmacokinetic and pharmacodynamic parameters across NHP brain regions - summary.')
write_supp_table(pk_out, 'Pharmacokinetic parameters across NHP brain regions - individual animals.')
write_supp_table(pd_out, 'Pharmacodynamic parameters across NHP brain regions - individual animals.')
  

ana %>%
  inner_join(types, by='cell_type') %>%
  inner_join(animals, by='animal') %>%
  inner_join(regions, by='region') %>%
  group_by(order, x, region, cell_type, animal, treatment) %>%
  summarize(.groups='keep',
            n_cells = n()) %>%
  ungroup() %>%
  select(-treatment) %>%
  pivot_wider(names_from=animal, values_from = n_cells, values_fill=0) %>%
  arrange(order, x) %>%
  select(-order, -x) -> counts_by_animal

write_supp_table(counts_by_animal, 'Nuclei counts per cell type in each animal.')

aco = read_tsv('data/modeling_outputs/aco.tsv', col_types=cols())
afi = read_tsv('data/modeling_outputs/afi.tsv', col_types=cols())
cell_type_composition = afi %>%
  mutate(region= sub(".*_", "", sample_id))%>%
  group_by(cell_type,category,region)%>%
  summarize(.groups='keep',
            n_cells=sum(n_cells),
            sum_umis= sum(total_umis)) %>%
  ungroup() %>%
  group_by(region)%>%
  mutate(
    total_cells = sum(n_cells),
    total_umi = sum(sum_umis),
    proportion_of_cells = n_cells / total_cells,
    proportion_of_umi = sum_umis / total_umi
  ) %>%
  ungroup()

write_supp_table(cell_type_composition, 'Cell type composition of each brain region by nuclei count and UMI count.')

write_supp_table(aco, 'Model fits of knockdown in each cell type + brain region combination.')
write_supp_table(afi, 'Point estimates for knockdown in each cell type + brain region + animal combination.')

bulk_sc %>%
  select(-tx, -pch) %>%
  rename(single_nucleus_residual = malat1_rel) %>%
  rename(bulk_qpcr_residual = qpcr) -> bulk_sc_out

bulk_sc_out %>%
  group_by(region, treatment) %>%
  summarize(.groups='keep',
            mean_snrnaseq = mean(single_nucleus_residual),
            l95_snrnaseq = lower(single_nucleus_residual),
            u95_snrnaseq = upper(single_nucleus_residual),
            mean_bulkqpcr = mean(bulk_qpcr_residual),
            l95_bulkqpcr = lower(bulk_qpcr_residual),
            u95_bulkqpcr = upper(bulk_qpcr_residual)) %>%
  ungroup() -> bulk_sc_smry_out

write_supp_table(bulk_sc_smry, 'Overall knockdown by bulk qPCR and by single nucleus sequencing - summary.')
write_supp_table(bulk_sc_out,  'Overall knockdown by bulk qPCR and by single nucleus sequencing - individual animal.')

# SUPPLEMENT #### 

tell_user('done.\nFinalizing supplementary tables...')

# write the supplement directory / table of contents
supplement_directory %>% rename(table_number = name, description=title) -> contents
addWorksheet(supplement,'contents')
bold_style = createStyle(textDecoration = "Bold")
writeData(supplement,'contents',contents,headerStyle=bold_style,withFilter=T)
freezePane(supplement,'contents',firstRow=T)
# move directory to the front
original_order = worksheetOrder(supplement)
n_sheets = length(original_order)
new_order = c(n_sheets, 1:(n_sheets-1))
worksheetOrder(supplement) = new_order
activeSheet(supplement) = 'contents'
# now save
saveWorkbook(supplement,supplement_path,overwrite = TRUE)

elapsed_time = Sys.time() - overall_start_time
cat(file=stderr(), paste0('done.\nAll tasks complete in ',round(as.numeric(elapsed_time),1),' ',units(elapsed_time),'.\n'))

