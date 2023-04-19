# R script to generate figures and related data for Hoffmann et al. 2023

# MIT License
# Copyright (C) 2023 Daniel Schmidt
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
  
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


### Preliminaries ####
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(colorspace)
library(ggpubr)
library(ComplexHeatmap)
library(rstatix)
library(NbClust)
library(bio3d)

# output directory for plots
pd = './plots/'

#source helper functions
source('./helper.R')

### Load AAV_alignment####
# this spreadsheet contain the alignment of the VP1 construct used in this study
# and various AAV serotypes and the AAV-DJ cryo-em structure (7KFR)
vp1_ref = openxlsx::read.xlsx('./input/VP1_reference.xlsx')

# load residues that are defined as external and internal to the capsid from
# the corresponding PDB files
vp1_ref$inside_res = rep(0, nrow(vp1_ref))
vp1_ref$outside_res = rep(0, nrow(vp1_ref))
inside_res = read.pdb('./input/inside_res.pdb')
outside_res = read.pdb(('./input/outside_res.pdb'))

inside_resno = unique(inside_res$atom$resno)
outside_resno = unique(outside_res$atom$resno)

idx = match(inside_resno, vp1_ref$resno_7kfr)
vp1_ref$inside_res[idx] = 1
vp1_ref$inside_res[vp1_ref$resno_vp1_ols < 218] = 1

idx = match(outside_resno, vp1_ref$resno_7kfr)
vp1_ref$outside_res[idx] = 1

#define boundaries of VP1u and PLA2 domain
vp1_ref$pla2 = ifelse(vp1_ref$resno_vp1_ols %in% c(52:97), 1, 0)
vp1_ref$vp1u = ifelse(vp1_ref$resno_vp1_ols %in% c(1:160), 1, 0)


### WT read count data ####
# from counting synonymous AAV-DJ mutation spiked into the input library
# order of data.frames:  input, packaging, pulldown, binding, uptake, infectivity
wt_1 = as.data.frame(t(c(51, 71, 104, 2682, 38, 55)))
wt_2 = as.data.frame(t(c(280, 289, 280, 192, 375, 236)))
colnames(wt_1) = c('input',
                   'packaging',
                   'pulldown',
                   'binding',
                   'uptake',
                   'infectivity')
colnames(wt_2) = c('input',
                   'packaging',
                   'pulldown',
                   'binding',
                   'uptake',
                   'infectivity')

wt = as.data.frame(t(rbind(wt_1, wt_2)))
colnames(wt) = c('rep1', 'rep2')

# wt confidence interval
# read counts (R1 & R2) for 10 synonymous AAV-DJ mutations
wt_s1_rep2 = c(125,
               159,
               137,
               137,
               139,
               131,
               177,
               137,
               162,
               159,
               240,
               214,
               175,
               132,
               158,
               144,
               135,
               108,
               20,
               8)
wt_s1_rep2_SE = qt(0.975, df = length(wt_s1_rep2) - 1 * sd(wt_s1_rep2) /
                     sqrt(length(wt_s1_rep2)))

### Load replicate 1 data ####
d1_FLAG = read.csv('./input/counts_FLAG_032822_cat.tsv',
                   sep = '\t',
                   header = T)
d1_gfp = read.csv('./input/counts_GFPnb_031622_cat.tsv',
                  sep = '\t',
                  header = T)
d1_mmobA = read.csv('./input/counts_mMobA_032822_cat.tsv',
                    sep = '\t',
                    header = T)
d1_DCV = read.csv('./input/counts_DCV_032822_cat.tsv',
                  sep = '\t',
                  header = T)
d1_WDV = read.csv('./input/counts_WDV_032822_cat.tsv',
                  sep = '\t',
                  header = T)
d1_SNAP = read.csv('./input/counts_SNAP_032822_cat.tsv',
                   sep = '\t',
                   header = T)
d1_SpyCatcher = read.csv('./input/counts_SpyCatcher_032822_cat.tsv',
                         sep = '\t',
                         header = T)
d1 = rbind(d1_FLAG,
           d1_gfp,
           d1_mmobA,
           d1_DCV,
           d1_WDV,
           d1_SNAP,
           d1_SpyCatcher)
d1 = d1[nchar(d1$counts) > 1, ]
d1$counts = as.list(strsplit(d1$counts, ','))
d1$target = as.factor(d1$target)
d1$domain = as.factor(d1$domain)
d1$cond = as.factor(d1$cond)
d1$rep = as.factor(rep('rep1', nrow(d1)))

### Load replicate 2 data ####
d2_FLAG = read.csv('./input/counts_FLAG_062722_cat.tsv',
                   sep = '\t',
                   header = T)
d2_gfp = read.csv('./input/counts_GFPnb_062722_cat.tsv',
                  sep = '\t',
                  header = T)
d2_mmobA = read.csv('./input/counts_mMobA_062722_cat.tsv',
                    sep = '\t',
                    header = T)
d2_DCV = read.csv('./input/counts_DCV_062722_cat.tsv',
                  sep = '\t',
                  header = T)
d2_WDV = read.csv('./input/counts_WDV_062722_cat.tsv',
                  sep = '\t',
                  header = T)
d2_SNAP = read.csv('./input/counts_SNAP_062722_cat.tsv',
                   sep = '\t',
                   header = T)
d2_SpyCatcher = read.csv('./input/counts_SpyCatcher_062722_cat.tsv',
                         sep = '\t',
                         header = T)
d2 = rbind(d2_FLAG,
           d2_gfp,
           d2_mmobA,
           d2_DCV,
           d2_WDV,
           d2_SNAP,
           d2_SpyCatcher)
d2 = d2[nchar(d2$counts) > 1, ]
d2$counts = as.list(strsplit(d2$counts, ','))
d2$target = as.factor(d2$target)
d2$domain = as.factor(d2$domain)
d2$cond = as.factor(d2$cond)
d2$rep = as.factor(rep('rep2', nrow(d2)))

### Combine and reorganize data####
d = rbind(d1, d2) #combined both replicates into one data frame

max_positions = max(unlist(lapply(d$counts, length)))
df = data.frame()
for (sdomain in levels(d$domain)) {
  for (scond in levels(d$cond)) {
    for (srep in levels(d$rep)) {
      subdf = droplevels(subset(d, cond == scond &
                                  domain == sdomain & rep == srep))
      domaindf = data.frame(
        'domain' = rep(sdomain, max_positions),
        'cond' = rep(scond, max_positions),
        'rep' = rep(srep, max_positions),
        'position' = seq(1, max_positions)
      )
      # If: there is more than one row for sample. Else: just unpack counts
      counts = as.numeric(unlist(subdf$counts))
      counts = c(counts, rep(0, max_positions - length(counts)))
      
      # add unpacked counts to domain dataframe
      domaindf = cbind(domaindf, counts) # raw counts
      # this method of finding a frequency for each position is for log-fold based enrichment
      domaindf = cbind(domaindf, tmp = (counts + 0.5) / sum(counts + 0.5)) # add 0.5 to all counts to assist positions with 0
      names(domaindf)[length(names(domaindf))] = 'counts_percent'
      df = rbind.fill(df, domaindf)
    }
  }
}

# remove data for extra measured phenotype assays not analyzed in this study
df = df[df$cond %!in% c('S7', 'S8', 'S9'), ]

# rename assay names to be human readable
df$cond = plyr::mapvalues(
  df$cond,
  from = c('S1', 'S2', 'S3', 'S4', 'S5', 'S6'),
  to = c(
    'input',
    'packaging',
    'pulldown',
    'binding',
    'uptake',
    'infectivity'
  )
)
df$cond = factor(
  df$cond,
  levels = c(
    'input',
    'packaging',
    'pulldown',
    'binding',
    'uptake',
    'infectivity'
  )
)

### Read depth ####
ggplot(df, aes(x = counts, color = domain)) +
  stat_ecdf(pad = FALSE) +
  scale_x_log10() +
  geom_vline(xintercept = 50,
             color = 'black',
             linetype = 'dashed') +
  facet_grid(cond ~ rep, scales = 'free') +
  ggtitle('Read Depth')
ggsave(paste0(pd,'read_depth.pdf'), device = 'pdf')

### Read count coverage by position and domain ####
ggplot(df, aes(x = position, y = counts)) +
  geom_col() +
  scale_y_log10() +
  facet_grid(cond ~ domain, scales = 'free') +
  theme_light()
ggsave(paste0(pd,'counts_coverage_by_pos_and_domain.pdf'), device = 'pdf')

### Replicate summary read stats ####
mapped_reads = df %>%
  group_by(cond, rep, domain) %>%
  summarise(mapped_reads = sum(counts))
write.csv(mapped_reads, './output/mapped_reads.csv')

### Data completenes s####
# mark positions with more than 20 reads
df$is.read = df$counts > 20

# completeness
completeness = df[df$position %!in% c(428:445), ] %>%
  group_by(cond, domain, rep) %>%
  summarise(fract = 100 - sum(!is.read) / sum(is.read) * 100)
paste0('median completeness: ', median(completeness$fract), '%')

ggplot(df, aes(x = cond, y = 1, fill = is.read)) +
  geom_bar(position = "fill", stat = "identity") +
  ylab('Percent of Positions') +
  facet_grid(rep ~ domain) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.text = element_text(size = 7)
  )
ggsave(paste0(pd,'read_threshold.pdf'), device = 'pdf')

# use only positions with reads above threshold (20)
# do not consider missing position 428-445 (library construction error)
df = df[df$is.read & df$position %!in% c(428:445), ]

## Median read count ####
med_readcounts_rep1 = df[df$rep == 'rep1', ] %>%
  group_by(cond, domain) %>%
  summarise(med_read = median(counts, na.rm = T),
            med_sd = sd(counts, na.rm = T))
range(med_readcounts_rep1$med_read)
print(med_readcounts_rep1, n = 42)
write.csv(med_readcounts_rep1, './output/med_readcounts_rep1_bydomain.csv')

med_readcounts_rep2 = df[df$rep == 'rep2', ] %>%
  group_by(cond, domain) %>%
  summarise(med_read = median(counts, na.rm = T),
            med_sd = sd(counts, na.rm = T))
range(med_readcounts_rep2$med_read)
print(med_readcounts_rep2, n = 42)
write.csv(med_readcounts_rep2, './output/med_readcounts_rep2_bydomain.csv')

### Replicate reproducibility ####

# 90% confidence interval for read count data based on Poisson distribution with 0.05% tail prob.
df[,c('counts_low','counts_high')] = t(sapply(df$counts, function(x) poisson.test(x, conf.level = 0.95)$conf.int[1:2]))
df$conf = (df$counts_high - df$counts_low) / df$counts

reprod_pl_data_p1 = dcast(df,
                          position + domain + cond ~ rep, value.var = 'counts_percent')
colnames(reprod_pl_data_p1) = c('position',
                                'domain',
                                'cond',
                                'counts_percent_rep1',
                                'counts_percent_rep2')
reprod_pl_data_p2 = dcast(df,
                          position + domain + cond ~ rep, value.var = 'conf')
colnames(reprod_pl_data_p2) = c('position', 'domain', 'cond', 'conf_rep1', 'conf_rep2')
reprod_pl_data = merge(reprod_pl_data_p1, reprod_pl_data_p2)
reprod_pl_data$conf = apply(reprod_pl_data[, c('conf_rep1', 'conf_rep2')], 1, var)

ggplot(reprod_pl_data,
       aes(x = counts_percent_rep1, y = counts_percent_rep2)) +
  geom_abline(slope = 1,
              intercept = 0,
              color = 'grey50') +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor() +
  geom_point(size = .5, alpha = 0.2) +
  scale_color_viridis_b() +
  facet_grid(cond ~ domain) +
  xlab('Replicate 1 Count Frequency') +
  ylab('Replicate 2 Count Frequency') +
  coord_fixed() +
  theme_light() +
  rotate_x_text() +
  xlim(0, 0.01) + ylim(0, 0.01) +
  ggtitle('replicate reproducibility')
ggsave(paste0(pd,'replicate_reproducibility.pdf'), device = 'pdf')

### Fitness calculation ####
# reorganize data.frame into a list
dfc = list(
  "FLAG" = dcast(
    df[df$domain == "FLAG" &
         df$rep == 'rep2', ],
    position ~ cond + rep,
    value.var = 'counts',
    fun.aggregate = sum
  ),
  "GFPnb" = dcast(
    df[df$domain == "GFPnb" &
         df$rep == 'rep2', ],
    position ~ cond + rep,
    value.var = 'counts',
    fun.aggregate = sum
  ),
  "SNAP" = dcast(
    df[df$domain == "SNAP" &
         df$rep == 'rep2', ],
    position ~ cond + rep,
    value.var = 'counts',
    fun.aggregate = sum
  ),
  "SpyCatcher" = dcast(
    df[df$domain == "SpyCatcher" &
         df$rep == 'rep2', ],
    position ~ cond + rep,
    value.var = 'counts',
    fun.aggregate = sum
  ),
  "WDV" = dcast(
    df[df$domain == "WDV" &
         df$rep == 'rep2', ],
    position ~ cond + rep,
    value.var = 'counts',
    fun.aggregate = sum
  ),
  "DCV" = dcast(
    df[df$domain == "DCV" &
         df$rep == 'rep2', ],
    position ~ cond + rep,
    value.var = 'counts',
    fun.aggregate = sum
  ),
  "mMobA" = dcast(
    df[df$domain == "mMobA" &
         df$rep == 'rep2', ],
    position ~ cond + rep,
    value.var = 'counts',
    fun.aggregate = sum
  )
)

# allocate the list that hold fitness seggregated by doman
dfc_enr = list()

for (i in names(dfc)) { # loop over all inserted domains
  tryCatch({ # loop error handling
    dfc_enr[[i]] = dfc[[i]] # copy over pre-fitness calculation data frame
    
    # packaging fitness and standard error
    dfc_enr[[i]]$W_packaging = log(dfc[[i]]$packaging_rep2 / dfc[[i]]$input_rep2 *
                                     wt_2$input / wt_2$packaging)
    dfc_enr[[i]]$W_packaging_SE = sqrt(
      1 / (dfc[[i]]$input_rep2 + 0.5) + 1 / (wt_2$input + 0.5) + 1 / (dfc[[i]]$packaging_rep2 + 0.5) + 1 /
        (wt_2$packaging + 0.5)
    )
    
    # pulldow and standard error (relative to packaged AAV)
    dfc_enr[[i]]$W_pulldown = log(dfc[[i]]$pulldown_rep2 / dfc[[i]]$packaging_rep2 *
                                    wt_2$packaging / wt_2$pulldown)
    dfc_enr[[i]]$W_pulldown_SE = sqrt(
      1 / (dfc[[i]]$packaging_rep2 + 0.5) + 1 / (wt_2$packaging + 0.5) + 1 / (dfc[[i]]$pulldown_rep2 + 0.5) + 1 /
        (wt_2$pulldown + 0.5)
    )
    
    # binding and standard error (relative to packaged AAV)
    dfc_enr[[i]]$W_binding = log(dfc[[i]]$binding_rep2 / dfc[[i]]$packaging_rep2 *
                                   wt_2$packaging / wt_2$binding)
    dfc_enr[[i]]$W_binding_SE = sqrt(
      1 / (dfc[[i]]$packaging_rep2 + 0.5) + 1 / (wt_2$packaging + 0.5) + 1 / (dfc[[i]]$binding_rep2 + 0.5) + 1 /
        (wt_2$binding + 0.5)
    )
    
    # uptake and standard error (relative to packaged AAV)
    dfc_enr[[i]]$W_uptake = log(dfc[[i]]$uptake_rep2 / dfc[[i]]$packaging_rep2 *
                                  wt_2$packaging / wt_2$uptake)
    dfc_enr[[i]]$W_uptake_SE = sqrt(
      1 / (dfc[[i]]$packaging_rep2 + 0.5) + 1 / (wt_2$packaging + 0.5) + 1 / (dfc[[i]]$uptake_rep2 + 0.5) + 1 /
        (wt_2$uptake + 0.5)
    )
    
    # infectivity and standard error(relative to packaged AAV)
    dfc_enr[[i]]$W_infectivity = log(dfc[[i]]$infectivity_rep2 / dfc[[i]]$packaging_rep2 *
                                       wt_2$packaging / wt_2$infectivity)
    dfc_enr[[i]]$W_infectivity_SE = sqrt(
      1 / (dfc[[i]]$packaging_rep2 + 0.5) + 1 / (wt_2$packaging + 0.5) + 1 / (dfc[[i]]$infectivity_rep2 + 0.5) + 1 /
        (wt_2$infectivity + 0.5)
    )
    
    # data polishing; replace infinity with NA, replace NaN with NA
    dfc_enr[[i]][dfc_enr[[i]] == Inf] = NA
    dfc_enr[[i]][dfc_enr[[i]] == -Inf] = NA
    dfc_enr[[i]][dfc_enr[[i]] == NaN] = NA
    dfc_enr[[i]]$gene = rep(i, nrow(dfc_enr[[i]]))
    
    # merge with other meta data from VP1 alignment
    dfc_enr[[i]] = left_join(dfc_enr[[i]], vp1_ref, by = c('position' = 'resno_vp1_ols'))
    
    print(paste0('Done with ', i))
  }, error = function(e) {
  })
}

# correction factor for binding, uptake, and infectivity fitness
# median binding fitness of FLAG insertion library (peaked around wildtype fitness = 0)
# this is necessitated by the HA tag inserted in to silent mutation DJ VP2+3,
# which slightly decreases binding fitness of silent mutation AAV-DJ compared to wt AAV-DJ
corr_factor = median(dfc_enr$FLAG$W_binding, na.rm = T)

for (i in names(dfc)) {
  dfc_enr[[i]]$W_binding = dfc_enr[[i]]$W_binding - corr_factor
  dfc_enr[[i]]$W_uptake = dfc_enr[[i]]$W_uptake - corr_factor
  dfc_enr[[i]]$W_infectivity = dfc_enr[[i]]$W_infectivity - corr_factor
}

### Map fitness onto AAV structure ####
s = Rpdb::read.pdb('./input/aav_map_target.pdb')

for (i in names(dfc_enr)) {
  idx = match(s$atoms$resid, dfc_enr[[i]]$resno_7kfr)
  
  # packaging fitness AAV_mapping
  s$atoms$temp = dfc_enr[[i]]$W_packaging[idx]
  s$atoms$temp[is.na(s$atoms$temp)] = 99 # B-factor set to this value for all NA fields
  Rpdb::write.pdb(s, file = paste0('./pdbs/', i, '_packaging.pdb'))
  
  # pulldown fitness AAV_mapping
  s$atoms$temp = dfc_enr[[i]]$W_pulldown[idx]
  s$atoms$temp[is.na(s$atoms$temp)] = 99
  Rpdb::write.pdb(s, file = paste0('./pdbs/', i, '_pulldown.pdb'))
  
  # binding fitness AAV_mapping
  s$atoms$temp = dfc_enr[[i]]$W_binding[idx]
  s$atoms$temp[is.na(s$atoms$temp)] = 99
  Rpdb::write.pdb(s, file = paste0('./pdbs/', i, '_binding.pdb'))
  
  # uptake fitness AAV_mapping
  s$atoms$temp = dfc_enr[[i]]$W_uptake[idx]
  s$atoms$temp[is.na(s$atoms$temp)] = 99
  Rpdb::write.pdb(s, file = paste0('./pdbs/', i, '_uptake.pdb'))
  
  # infectivity fitness AAV_mapping
  s$atoms$temp = dfc_enr[[i]]$W_infectivity[idx]
  s$atoms$temp[is.na(s$atoms$temp)] = 99
  Rpdb::write.pdb(s, file = paste0('./pdbs/', i, '_infectivity.pdb'))
  
}

### Insertion library titer ####
lib_titer = read.csv('./input/lib_titer.csv')
lib_titer = melt(
  lib_titer,
  id.vars = 'library',
  variable.name = 'replicate',
  value.name = 'titer'
)
lib_titer$library = as.factor(lib_titer$library)

ggboxplot(lib_titer, x = 'library', y = 'titer') +
  rotate_x_text()
ggsave(paste0(pd,'library_titer.pdf'), device = 'pdf')

# one-way ANOVA
lib_titer_model = aov(titer ~ library, data = lib_titer)
summary(lib_titer_model)
DescTools::DunnettTest(x = lib_titer$titer, g = lib_titer$library)

### Insertion library VP1 content ####
lib_vp1 = read.csv(
  './input/lib_vp1_content.csv',
  colClasses = c('factor', 'factor', 'factor', 'factor', 'numeric')
)

# two-way anova repeated measured ANOVA
lib_vp1.aov <- anova_test(
  data = lib_vp1,
  dv = content,
  wid = id,
  within = vp,
  between = library
)
get_anova_table(lib_vp1.aov)


# one way ANOVA testing for significant differences in VP1 content
lib_vp1.aov2 = aov(content ~ library, data = lib_vp1[lib_vp1$vp == 'VP1', ])
summary(lib_vp1.aov2)

lib_vp1.pw = lib_vp1 %>%
  group_by(vp) %>%
  pairwise_t_test(content ~ library,
                  paired = F,
                  p.adjust.method = "bonferroni")
print(lib_vp1.pw, n = 100)

DescTools::DunnettTest(x = lib_vp1[lib_vp1$vp == 'VP1', 'content'],
                       g = lib_vp1[lib_vp1$vp == 'VP1', 'library'],
                       control = 'wtDJ')

### Cysteine Mutants (all VPs) ####
cys_titer = read.csv('./input/cys_titer.csv')
cys_titer$library = as.factor(cys_titer$library)
cys_titer$replicate = as.factor(cys_titer$replicate)

ggboxplot(cys_titer, x = 'library', y = 'titer') +
  rotate_x_text()
ggsave(paste0(pd,'cysteine_mutant_titer.pdf'), device = 'pdf')

# one way ANOVA testing for significant differences in VP1 content
cys_titer_model = aov(titer ~ library, data = cys_titer)
summary(cys_titer_model)
DescTools::DunnettTest(x = cys_titer$titer,
                       g = cys_titer$library,
                       control = 'wt DJ')

### WDV crude lysate titer ####
wdv_titer = read.csv('./input/wdv_titer.csv')
wdv_titer$gene = as.factor(wdv_titer$gene)
wdv_titer$replicate = as.factor(wdv_titer$replicate)

ggboxplot(wdv_titer, x = 'gene', y = 'pur_titer') +
  rotate_x_text()
ggsave(paste0(pd,'WDV_titer_crude_lysate.pdf'), device = 'pdf')

wdv_titer_model = aov(pur_titer ~ gene, data = wdv_titer)
summary(wdv_titer_model)
DescTools::DunnettTest(x = wdv_titer$pur_titer,
                       g = wdv_titer$gene,
                       control = 'wtDJ')

### WDV purified titer ####
wdv_ptiter = read.csv('./input/wdv_pur_titer.csv')
wdv_ptiter$gene = as.factor(wdv_ptiter$gene)
wdv_ptiter$replicate = as.factor(wdv_ptiter$replicate)

ggboxplot(wdv_ptiter, x = 'gene', y = 'pur_titer')  +
  rotate_x_text()

wdv_ptiter_model = aov(pur_titer ~ gene, data = wdv_ptiter)
summary(wdv_ptiter_model)
DescTools::DunnettTest(x = wdv_ptiter$pur_titer,
                       g = wdv_ptiter$gene,
                       control = 'wtDJ')


### Heat maps of fitness scores (and SE) ####

# fitness values select for plotting
col_select_1 = c(
  'gene',
  'position',
  'W_packaging',
  'W_pulldown',
  'W_binding',
  'W_uptake',
  'W_infectivity'
)

dfc_enr_pivot = do.call(rbind, dfc_enr)[, col_select_1]
dfc_enr_pivot_melt = melt(dfc_enr_pivot, id.vars = c('gene', 'position'))
dfc_enr_pivot_melt$org_var = interaction(dfc_enr_pivot_melt$gene, dfc_enr_pivot_melt$variable)
dfc_enr_pivot = as.matrix(dcast(dfc_enr_pivot_melt, position ~ org_var, fun.aggregate = sum))

# fitness standard errors  select for plotting
col_select_2 = c(
  'gene',
  'position',
  'W_packaging_SE',
  'W_pulldown_SE',
  'W_binding_SE',
  'W_uptake_SE',
  'W_infectivity_SE'
)

dfc_se_pivot = do.call(rbind, dfc_enr)[, col_select_2]
dfc_se_pivot_melt = melt(dfc_se_pivot, id.vars = c('gene', 'position'))
dfc_se_pivot_melt$org_var = interaction(dfc_se_pivot_melt$gene, dfc_se_pivot_melt$variable)
dfc_se_pivot = as.matrix(dcast(dfc_se_pivot_melt, position ~ org_var, fun.aggregate = sum))

# raw ead counts selected for plotting
col_select_3 = c(
  'gene',
  'position',
  'input_rep2',
  'packaging_rep2',
  'pulldown_rep2',
  'binding_rep2',
  'uptake_rep2',
  'infectivity_rep2'
)

dfc_rc_pivot = do.call(rbind, dfc_enr)[, col_select_3]
dfc_rc_pivot_melt = melt(dfc_rc_pivot, id.vars = c('gene', 'position'))
dfc_rc_pivot_melt$org_var = interaction(dfc_rc_pivot_melt$gene, dfc_rc_pivot_melt$variable)
dfc_rc_pivot = as.matrix(dcast(dfc_rc_pivot_melt, position ~ org_var, fun.aggregate = sum))

rownames(dfc_enr_pivot) = dfc_enr_pivot[, 1]
rownames(dfc_se_pivot) = dfc_se_pivot[, 1]
rownames(dfc_rc_pivot) = dfc_rc_pivot[, 1]
dfc_enr_pivot = dfc_enr_pivot[, -1]
dfc_se_pivot = dfc_se_pivot[, -1]
dfc_rc_pivot = dfc_rc_pivot[, -1]
dfc_rc_pivot = log10(dfc_rc_pivot)
dfc_rc_pivot[dfc_rc_pivot == Inf |
               dfc_rc_pivot == -Inf | dfc_rc_pivot == NaN] = NA

dfc_enr_pivot[dfc_se_pivot > 0.5] = NA

col_fun1 = circlize::colorRamp2(c(-1,-0.75,-0.5, 0, 0.5, 0.75, 1),
                                divergingx_hcl(7, palette = 'PiYG'))
col_fun2 = circlize::colorRamp2(seq(0, 0.5, by = 0.1),
                                sequential_hcl(6, palette = 'OrRd', rev = T))
col_fun3 = circlize::colorRamp2(seq(0, 4, by = 0.5),
                                sequential_hcl(9, palette = 'OrRd', rev = T))

##### median poisson error ####

# overall
dfc_se_pivot_melt %>%
  summarise(median_se = median(value, na.rm = T))

# by assay
dfc_se_pivot_melt %>%
  group_by(variable) %>%
  summarise(median_se = median(value, na.rm = T))

# Fitness heatmap
Heatmap(
  t(dfc_enr_pivot),
  cluster_rows = F,
  cluster_columns = F,
  na_col = 'yellow',
  col = col_fun1,
  heatmap_legend_param = list(at = c(-1, 0, 1),
                              title = 'Fitness'),
  column_names_gp = grid::gpar(fontsize = 1),
  row_names_gp = grid::gpar(fontsize = 4)
)

# Fitness standard error heatmap
Heatmap(
  t(dfc_se_pivot),
  cluster_rows = F,
  cluster_columns = F,
  na_col = 'yellow',
  col = col_fun2,
  heatmap_legend_param = list(at = seq(0, 0.5, by = 0.1),
                              title = 'Fitness Error'),
  column_names_gp = grid::gpar(fontsize = 1),
  row_names_gp = grid::gpar(fontsize = 4)
)

# Read count (log10) heatmap
Heatmap(
  t(dfc_rc_pivot),
  cluster_rows = F,
  cluster_columns = F,
  na_col = 'yellow',
  col = col_fun3,
  heatmap_legend_param = list(at = seq(0, 4, by = 0.5),
                              title = 'Log10 Read Count'),
  column_names_gp = grid::gpar(fontsize = 1),
  row_names_gp = grid::gpar(fontsize = 4)
)

## read count vs. fitness vs. position count (to establish dynamic range of the assays)
dfc_enr_ul = do.call(rbind, dfc_enr)
dfc_enr_ul_melt = melt(dfc_enr_ul[, col_select_1], id.vars = c('gene', 'position'))

ggplot(dfc_enr_ul_melt, aes(x = value, color = variable)) +
  stat_ecdf()

ggplot(dfc_enr_ul_melt, aes(x = value)) +
  geom_density(aes(color = variable), alpha = 0.25) +
  xlim(-2, 1) +
  facet_grid(variable ~ .) +
  xlab('Fitness') +
  guides(color=guide_legend(title='Phenotype Assay'))
ggsave(paste0(pd,'fitness_distribution.pdf'), device = 'pdf')


### Pulldown fitness for VP1u residues ####
ggplot(dfc_enr_ul,
       aes(x = W_pulldown, color = interaction(outside_res, vp1u))) +
  stat_ecdf() +
  facet_grid(gene ~ .) +
  xlim(-0.5, 0.5) +
  theme_classic2() +
  xlab('Pulldown Fitness') +
  guides(color=guide_legend(title='[outside residue].[within VP1u] [0/1]'))
ggsave(paste0(pd,'vp1u_vs_inside_pd.pdf'), device = 'pdf')

# test significance for distribution differences between inside residues 
# (that are not VP1u) and VP1
for (i in levels(as.factor(dfc_enr_ul$gene))) {
  tmp = ks.test(dfc_enr_ul$W_pulldown[dfc_enr_ul$outside_res == 0 &
                                        dfc_enr_ul$gene == i],
                dfc_enr_ul$W_pulldown[dfc_enr_ul$vp1u == 1 &
                                        dfc_enr_ul$gene == i],
                exact = T)$p.value
  print(paste0('Gene: ', i, ' ; inside vs vp1u p.value: ', tmp))
}

### Fitness at 2-,3-, and 5-fold interfaces ####
dfc_enr_ul_d = melt(
  dfc_enr_ul[, c(
    'position',
    'W_packaging',
    'W_pulldown',
    'W_binding',
    'W_uptake',
    'gene',
    'interface_2fold',
    'interface_3fold',
    'interface_5fold',
    'inside_res',
    'outside_res'
  )],
  id.vars = c(
    'position',
    'gene',
    'interface_2fold',
    'interface_3fold',
    'interface_5fold',
    'inside_res',
    'outside_res'
  ),
  variable.name = 'phenotype'
)

# interface membership is estblished using the pymol interfaceResidue script
# https://pymolwiki.org/index.php/InterfaceResidues
dfc_enr_ul_d$interface_2fold = as.factor(ifelse(is.na(dfc_enr_ul_d$interface_2fold), 0, 1))
dfc_enr_ul_d$interface_3fold = as.factor(ifelse(is.na(dfc_enr_ul_d$interface_3fold), 0, 1))
dfc_enr_ul_d$interface_5fold = as.factor(ifelse(is.na(dfc_enr_ul_d$interface_5fold), 0, 1))
dfc_enr_ul_d$resno_AAVDJwt = vp1_ref$resno_AAVDJwt[match(dfc_enr_ul_d$pos, vp1_ref$resno_vp1_ols)]

# test for significance of distribution difference (in interface vs outside interface)
stat.test_2f = dfc_enr_ul_d %>%
  group_by(phenotype, gene) %>%
  wilcox_test(value ~ interface_2fold) %>%
  adjust_pvalue() %>%
  add_significance()
stat.test_3f = dfc_enr_ul_d %>%
  group_by(phenotype, gene) %>%
  wilcox_test(value ~ interface_3fold) %>%
  adjust_pvalue() %>%
  add_significance()
stat.test_5f = dfc_enr_ul_d %>%
  group_by(phenotype, gene) %>%
  wilcox_test(value ~ interface_5fold) %>%
  adjust_pvalue() %>%
  add_significance()

stat.test_2f = stat.test_2f %>% add_xy_position(x = "interface_2fold")
stat.test_3f = stat.test_3f %>% add_xy_position(x = "interface_3fold")
stat.test_5f = stat.test_5f %>% add_xy_position(x = "interface_5fold")

# cumulative density plot of fitness within interface vs outside the interface

ggecdf(
  dfc_enr_ul_d,
  x = 'value',
  color = 'interface_2fold',
  facet.by = c('gene', 'phenotype')
) +
  xlim(-1, 1) +
  stat_pvalue_manual(stat.test_2f, label = 'p.adj.signif', y.position = 0.5) +
  xlab('Fitness') +
  guides(color=guide_legend(title='within 2-fold interface [0/1]'))
ggsave(paste0(pd,'packaging_fitness_2fold_interface.pdf'), device = 'pdf')

ggecdf(
  dfc_enr_ul_d,
  x = 'value',
  color = 'interface_3fold',
  facet.by = c('gene', 'phenotype')
) +
  xlim(-1, 1) +
  stat_pvalue_manual(stat.test_3f, label = 'p.adj.signif', y.position = 0.5) +
  xlab('Fitness') +
  guides(color=guide_legend(title='within 3-fold interface [0/1]'))
ggsave(paste0(pd,'packaging_fitness_3fold_interface.pdf'), device = 'pdf')

ggecdf(
  dfc_enr_ul_d,
  x = 'value',
  color = 'interface_5fold',
  facet.by = c('gene', 'phenotype')
) +
  xlim(-1, 1) +
  stat_pvalue_manual(stat.test_5f, label = 'p.adj.signif', y.position = 0.5) +
  xlab('Fitness') +
  guides(color=guide_legend(title='within 5-fold interface [0/1]'))
ggsave(paste0(pd,'packaging_fitness_5fold_interface.pdf'), device = 'pdf')

# interface membership overlap
interface_list = list('2fold' = dfc_enr_ul[dfc_enr_ul_d$interface_2fold == '1', 'resno_AAVDJwt'],
                      '3fold' = dfc_enr_ul[dfc_enr_ul_d$interface_3fold == '1', 'resno_AAVDJwt'],
                      '5fold' = dfc_enr_ul[dfc_enr_ul_d$interface_5fold == '1', 'resno_AAVDJwt'])

ggVennDiagram::ggVennDiagram(interface_list) +
  ggplot2::scale_fill_gradient(low = '#f2f0f7', high = '#6a51a3') +
  ggplot2::ggtitle('Interface Residues Membership Overlap')
ggsave(paste0(pd,'venn_interfaces.pdf'), device = 'pdf')

### UMAP clustering####

# select phenotype assays used for UMAP clustering
umap_select = c('gene',
                'position',
                'W_packaging',
                'W_pulldown',
                'W_binding',
                'W_uptake')
dfc_enr_pivot = do.call(rbind, dfc_enr)[, umap_select]
dfc_enr_pivot_melt = melt(dfc_enr_pivot, id.vars = c('gene', 'position'))
dfc_enr_pivot_melt$org_var = interaction(dfc_enr_pivot_melt$gene, dfc_enr_pivot_melt$variable)
dfc_enr_pivot = as.matrix(dcast(dfc_enr_pivot_melt, position ~ org_var, fun.aggregate = sum))

rownames(dfc_enr_pivot) = dfc_enr_pivot[, 1]
dfc_enr_pivot = dfc_enr_pivot[, -1]

W = t(dfc_enr_pivot)
W[is.na(W)] = 0
W[W == 'Inf' | W == '-Inf'] = 0
W = data.matrix(W)
umap_W = t(W)

set.seed(220)
dfc.umap = uwot::umap(
  umap_W,
  metric = 'cosine',
  n_epochs = 500,
  n_neighbors = 20,
  nn_method = 'annoy',
  n_trees = 50
)


dfc.umap.scores = data.frame(dfc.umap) # PC score matrix
# find ideal cluster number
dfc.umap.scores$group = NbClust::NbClust(dfc.umap.scores, method = 'ward.D2', min.nc = 3)$Best.partition
dfc.umap.scores$position = as.numeric(colnames(W))

ggplot(dfc.umap.scores, aes(X1, X2, color = as.factor(group))) +
  geom_point() +
  scale_color_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  coord_fixed() +
  xlab('UMAP1') +
  ylab('UMAP2') +
  guides(color=guide_legend(title='cluster')) +
  ggtitle('UMAP clustering of fitness')
ggsave(paste0(pd,'umap_clusters.pdf'), device = 'pdf')


ggplot(dfc.umap.scores, aes(x = position, y = 1, fill = as.factor(group))) +
  geom_tile() +
  scale_x_continuous(breaks = seq(0, max(dfc.umap.scores$position), by = 20)) + ggtitle('all phenotypes') +
  theme_classic2() +
  #  facet_grid(as.factor(group)~.) +
  scale_fill_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  theme(axis.text.x = element_text(angle = 90)) +
  guides(color=guide_legend(title='cluster'))
ggsave(paste0(pd,'umap_cluster_by_pos.pdf'), device = 'pdf')

# map UMAP cluster on structure
dfc.umap.scores$pos_7kfr = vp1_ref$resno_AAVDJwt[match(dfc.umap.scores$position, vp1_ref$resno_vp1_ols)]
wcpdb(data_source = 'dfc.umap.scores', data_col = 'group',
      file_in = './input/aav_map_target.pdb',
      file_out = './pdbs/umap_clusters.pdb')

dfc_enr_pivot_melt = merge(dfc_enr_pivot_melt, dfc.umap.scores[, c('position', 'pos_7kfr', 'group')], all.x = T)

ggplot(dfc_enr_pivot_melt,
       aes(x = value)) +
  geom_histogram(
    aes(fill = as.factor(group), y = after_stat(ndensity)),
    position = 'identity',
    bins = 100,
    alpha = 0.7
  ) +
  scale_fill_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  geom_vline(xintercept = 0) +
  facet_grid(group ~ variable) +
  guides(fill=guide_legend(title='cluster')) +
  xlab('Fitness')
ggsave(paste0(pd,'umap_cluster_fitness_distribution.pdf'), device = 'pdf')


### Fitness Variance ####

# select phenotype assays used for fitness variance mapping
col_select_var = c('gene',
                   'position',
                   'W_packaging',
                   'W_pulldown',
                   'W_binding',
                   'W_uptake')

dfc_enr_pivot = do.call(rbind, dfc_enr)[, col_select_var]
rownames(dfc_enr_pivot) = NULL

#### Packaging Variance ####
packaging_var = dfc_enr_pivot %>%
  group_by(position) %>%
  summarize(
    W_packaging_mean = mean(W_packaging, na.rm = T) + 0.001,
    W_packaging_sd = sd(W_packaging, na.rm = T) + 0.001,
    W_packaging_var = var(W_packaging, na.rm = T) + 0.001
  )

# rolling mean to smooth out the data (window size 5 residues)
packaging_var$pos_7kfr = vp1_ref$resno_AAVDJwt[match(packaging_var$position, vp1_ref$resno_vp1_ols)]
packaging_var$W_packaging_mean_roll = zoo::rollmean(packaging_var$W_packaging_mean, 5, na.pad = T)
packaging_var$W_packaging_var_roll = zoo::rollmean(packaging_var$W_packaging_var, 5, na.pad = T)

# merge UMAP assignment into dataframe
packaging_var = merge(packaging_var, dfc.umap.scores[, c('position', 'group')], all.x = T)

p1 = ggplot(packaging_var,
            aes(
              x = position,
              y = W_packaging_mean_roll,
              color = as.factor(group)
            )) +
  geom_point() +
  scale_color_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  ylab('Packaging Fitness, rolling mean')
  theme_classic2()

p2 = ggplot(packaging_var, aes(x = W_packaging_mean_roll, fill = as.factor(group))) +
  geom_histogram(aes(y = after_stat(ndensity)), position = 'identity', bins = 50) +
  scale_fill_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  facet_grid(group ~ .) +
  xlab('Packaging Fitness, rolling mean') +
  theme_classic2()

p3 = ggplot(packaging_var,
            aes(
              x = position,
              y = W_packaging_var_roll,
              color = as.factor(group)
            )) +
  geom_point() +
  scale_color_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  ylab('Packaging Fitness, rolling variance')
  theme_classic2()

p4 = ggplot(packaging_var, aes(x = W_packaging_var_roll, fill = as.factor(group))) +
  geom_histogram(aes(y = after_stat(ndensity)), position = 'identity', bins = 50) +
  scale_fill_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  facet_grid(group ~ .) +
  xlab('Packaging Fitness, rolling variance') +
  theme_classic2()

ggarrange(p1,
          p2,
          p3,
          p4,
          ncol = 2,
          nrow = 2,
          common.legend = T)
ggsave(paste0(pd,'packaging_var_by_pos_and_umap_group.pdf'), device = 'pdf')

# Map packaging variance onto AAV structure
wcpdb(data_source = 'packaging_var', data_col = 'W_packaging_var_roll', na_val = 0,
      file_in = './input/aav_map_target.pdb',
      file_out = './pdbs/packaging_var_roll.pdb')

#### Pulldown Variance ####
pulldown_var = dfc_enr_pivot %>%
  group_by(position) %>%
  summarize(
    W_pulldown_mean = mean(W_pulldown, na.rm = T) + 0.001,
    W_pulldown_sd = sd(W_pulldown, na.rm = T) + 0.001,
    W_pulldown_var = var(W_pulldown, na.rm = T) + 0.001
  )

# rolling mean to smooth out the data (window size 5 residues)
pulldown_var$pos_7kfr = vp1_ref$resno_AAVDJwt[match(pulldown_var$position, vp1_ref$resno_vp1_ols)]
pulldown_var$W_pulldown_mean_roll = zoo::rollmean(pulldown_var$W_pulldown_mean, 5, na.pad = T)
pulldown_var$W_pulldown_var_roll = zoo::rollmean(pulldown_var$W_pulldown_var, 5, na.pad = T)

# merge UMAP assignment into dataframe
pulldown_var = merge(pulldown_var, dfc.umap.scores[, c('position', 'group')], all.x = T)

p1 = ggplot(pulldown_var,
            aes(
              x = position,
              y = W_pulldown_mean_roll,
              color = as.factor(group)
            )) +
  geom_point() +
  scale_color_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  ylab('Pulldown Fitness, rolling mean')
  theme_classic2()

p2 = ggplot(pulldown_var, aes(x = W_pulldown_mean_roll, fill = as.factor(group))) +
  geom_histogram(aes(y = after_stat(ndensity)), position = 'identity', bins = 50) +
  scale_fill_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  facet_grid(group ~ .) +
  xlab('Pulldown Fitness, rolling mean') +
  theme_classic2()

p3 = ggplot(pulldown_var,
            aes(
              x = position,
              y = W_pulldown_var_roll,
              color = as.factor(group)
            )) +
  geom_point() +
  scale_color_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  ylab('Pulldown Fitness, rolling variance')
  theme_classic2()

p4 = ggplot(pulldown_var, aes(x = W_pulldown_var_roll, fill = as.factor(group))) +
  geom_histogram(aes(y = after_stat(ndensity)), position = 'identity', bins = 50) +
  scale_fill_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  facet_grid(group ~ .) +
  xlab('Pulldown Fitness, rolling variance') +
  theme_classic2()

ggarrange(p1,
          p2,
          p3,
          p4,
          ncol = 2,
          nrow = 2,
          common.legend = T)
ggsave(paste0(pd,'pulldown_var_by_pos_and_umap_group.pdf'), device = 'pdf')

#Map pulldown variance onto AAV structure
wcpdb(data_source = 'pulldown_var', data_col = 'W_pulldown_var_roll', na_val = 0,
      file_in = './input/aav_map_target.pdb',
      file_out = './pdbs/pulldown_var_roll.pdb')

#### Binding Variance ####
binding_var = dfc_enr_pivot %>%
  group_by(position) %>%
  summarize(
    W_binding_mean = mean(W_binding, na.rm = T) + 0.001,
    W_binding_sd = sd(W_binding, na.rm = T) + 0.001,
    W_binding_var = var(W_binding, na.rm = T) + 0.001
  )

# rolling mean to smooth out the data (window size 5 residues)
binding_var$pos_7kfr = vp1_ref$resno_AAVDJwt[match(binding_var$position, vp1_ref$resno_vp1_ols)]
binding_var$W_binding_mean_roll = zoo::rollmean(binding_var$W_binding_mean, 5, na.pad = T)
binding_var$W_binding_var_roll = zoo::rollmean(binding_var$W_binding_var, 5, na.pad = T)

# merge UMAP assignment into dataframe
binding_var = merge(binding_var, dfc.umap.scores[, c('position', 'group')], all.x = T)

p1 = ggplot(binding_var,
            aes(
              x = position,
              y = W_binding_mean_roll,
              color = as.factor(group)
            )) +
  geom_point() +
  scale_color_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  ylab('Binding Fitness, rolling mean')
  theme_classic2()

p2 = ggplot(binding_var, aes(x = W_binding_mean_roll, fill = as.factor(group))) +
  geom_histogram(aes(y = after_stat(ndensity)), position = 'identity', bins = 50) +
  scale_fill_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  facet_grid(group ~ .) +
  xlab('Binding Fitness, rolling mean') +
  theme_classic2()

p3 = ggplot(binding_var,
            aes(
              x = position,
              y = W_binding_var_roll,
              color = as.factor(group)
            )) +
  geom_point() +
  scale_color_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  ylab('Binding Fitness, rolling variance')
  theme_classic2()

p4 = ggplot(binding_var, aes(x = W_binding_var_roll, fill = as.factor(group))) +
  geom_histogram(aes(y = after_stat(ndensity)), position = 'identity', bins = 50) +
  scale_fill_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  facet_grid(group ~ .) +
  xlab('Binding Fitness, rolling variance') +
  theme_classic2()

ggarrange(p1,
          p2,
          p3,
          p4,
          ncol = 2,
          nrow = 2,
          common.legend = T)
ggsave(paste0(pd,'binding_var_by_pos_and_umap_group.pdf'), device = 'pdf')

# Map binding variance onto AAV structure
wcpdb(data_source = 'binding_var', data_col = 'W_binding_var_roll', na_val = 0,
      file_in = './input/aav_map_target.pdb',
      file_out = './pdbs/binding_var_roll.pdb')

#### Uptake Variance ####
uptake_var = dfc_enr_pivot %>%
  group_by(position) %>%
  summarize(
    W_uptake_mean = mean(W_uptake, na.rm = T) + 0.001,
    W_uptake_sd = sd(W_uptake, na.rm = T) + 0.001,
    W_uptake_var = var(W_uptake, na.rm = T) + 0.001
  )

# rolling mean to smooth out the data (window size 5 residues)
uptake_var$pos_7kfr = vp1_ref$resno_AAVDJwt[match(uptake_var$position, vp1_ref$resno_vp1_ols)]
uptake_var$W_uptake_mean_roll = zoo::rollmean(uptake_var$W_uptake_mean, 5, na.pad = T)
uptake_var$W_uptake_var_roll = zoo::rollmean(uptake_var$W_uptake_var, 5, na.pad = T)

# merge UMAP assignment into dataframe
uptake_var = merge(uptake_var, dfc.umap.scores[, c('position', 'group')], all.x = T)

p1 = ggplot(uptake_var,
            aes(
              x = position,
              y = W_uptake_mean_roll,
              color = as.factor(group)
            )) +
  geom_point() +
  scale_color_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  ylab('Uptake Fitness, rolling mean')
  theme_classic2()

p2 = ggplot(uptake_var, aes(x = W_uptake_mean_roll, fill = as.factor(group))) +
  geom_histogram(aes(y = after_stat(ndensity)), position = 'identity', bins = 50) +
  scale_fill_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  facet_grid(group ~ .) +
  xlab('Uptake Fitness, rolling mean') +
  theme_classic2()

p3 = ggplot(uptake_var,
            aes(
              x = position,
              y = W_uptake_var_roll,
              color = as.factor(group)
            )) +
  geom_point() +
  scale_color_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  ylab('Uptake Fitness, rolling variance')
  theme_classic2()

p4 = ggplot(uptake_var, aes(x = W_uptake_var_roll, fill = as.factor(group))) +
  geom_histogram(aes(y = after_stat(ndensity)), position = 'identity', bins = 50) +
  scale_fill_manual(values = c('#FB4B5A', '#FFB347', '#0FBD9A', '#D3BFD9', '#0F497E')) +
  facet_grid(group ~ .) +
  xlab('Uptake Fitness, rolling variance') +
  theme_classic2()

ggarrange(p1,
          p2,
          p3,
          p4,
          ncol = 2,
          nrow = 2,
          common.legend = T)
ggsave(paste0(pd,'uptake_var_by_pos_and_umap_group.pdf'), device = 'pdf')

# Map uptake variance onto AAV structure
wcpdb(data_source = 'uptake_var', data_col = 'W_uptake_var_roll', na_val = 0,
      file_in = './input/aav_map_target.pdb',
      file_out = './pdbs/uptake_var_roll.pdb')

var_list =
  list(packaging_var, pulldown_var, binding_var, uptake_var)
var_df = var_list %>% purrr::reduce(full_join, by = c('position', 'pos_7kfr', 'group'))

ggplot(var_df, aes(x = position)) +
  #geom_col(aes(y = 0.25, fill = as.factor(group))) +
  geom_line(aes(y = W_packaging_mean_roll), color = '#D9496F') +
  geom_line(aes(y = W_pulldown_mean_roll), color = '#66997D') +
  geom_line(aes(y = W_binding_mean_roll), color = '#189BDC') +
  geom_line(aes(y = W_uptake_mean_roll), color = '#946699') +
  ylab('Fitness, rolling mean')
  theme_pubclean()
ggsave(paste0(pd,'mean_fitnness_by_pos.pdf'), device = 'pdf')


### WDV N664 & K708 infectivity ####

# take in infectivity data
d = read.csv('./input/wdv_anova_v1.csv',
             colClasses = c('factor',
                            'factor',
                            'factor',
                            'numeric'))

# two-way repeated measure ANOVA for conditions without GFP-GPI expression
res.aov_woGFP = anova_test(
  data = d[d$cond %in% c('4', '4a'), ],
  dv = tdt_pos,
  wid = id,
  within = cond,
  between = gene
)
get_anova_table(res.aov_woGFP)

# Bonferroni post test
pwc_woGFP = d %>%
  filter(cond %in% c('4', '4a')) %>%
  group_by(gene) %>%
  pairwise_t_test(tdt_pos ~ cond, paired = TRUE,
                  p.adjust.method = "bonferroni")
pwc_woGFP

# two-way repeated measure ANOVA for conditions with GFP-GPI expression
res.aov_wGFP <- anova_test(
  data = d[d$cond %in% c('4g', '4ga'), ],
  dv = tdt_pos,
  wid = id,
  within = cond,
  between = gene
)
get_anova_table(res.aov_wGFP)

# Bonferroni post test
pwc_wGFP = d %>%
  filter(cond %in% c('4g', '4ga')) %>%
  group_by(gene) %>%
  pairwise_t_test(tdt_pos ~ cond, paired = TRUE,
                  p.adjust.method = "bonferroni")
pwc_wGFP


p1 = ggpaired(
  d[d$gene == 'wtDJ', ],
  x = "cond",
  y = "tdt_pos",
  color = "cond",
  line.color = "gray",
  line.size = 0.4,
  #facet.by = 'gene',
  palette = "jco"
) +
  rotate_x_text() +
  stat_compare_means(
    method = 't.test',
    paired = T,
    comparisons = list(c('4', '4a'), c('4g', '4ga')),
    method.args = list(p.adjust.method = "bonferroni")
  ) +
  ylab('% tdtomota positive') +
  ylim(0, 100) +
  ggtitle('WT')

p2 = ggpaired(
  d[d$gene == 'N664', ],
  x = "cond",
  y = "tdt_pos",
  color = "cond",
  line.color = "gray",
  line.size = 0.4,
  #facet.by = 'gene',
  palette = "jco"
) +
  rotate_x_text() +
  stat_compare_means(
    method = 't.test',
    paired = T,
    comparisons = list(c('4', '4a'), c('4g', '4ga')),
    method.args = list(p.adjust.method = "bonferroni")
  ) +
  ylab('% tdtomota positive') +
  ylim(0, 100) +
  ggtitle('N664')

p3 = ggpaired(
  d[d$gene == 'K708', ],
  x = "cond",
  y = "tdt_pos",
  color = "cond",
  line.color = "gray",
  line.size = 0.4,
  #facet.by = 'gene',
  palette = "jco"
) +
  rotate_x_text() +
  stat_compare_means(
    method = 't.test',
    paired = T,
    comparisons = list(c('4', '4a'), c('4g', '4ga')),
    method.args = list(p.adjust.method = "bonferroni")
  ) +
  ylab('% tdtomota positive') +
  ylim(0, 100) +
  ggtitle('K708')


ggarrange(p1, p2, p3, nrow = 1, common.legend = T)
ggsave(paste0(pd,'wdv_inf_comp.pdf'), device = 'pdf')


wdv_titer = read.csv('./input/wdv_variant_titers.csv',
                     colClasses = c('factor',
                                    'factor',
                                    'numeric'))
wdv_titer$gene = factor(
  wdv_titer$gene,
  levels = c(
    'wtDJ',
    'S244',
    'L256',
    'S268',
    'K311',
    'T331',
    'E349',
    'T414',
    'N664',
    'Y702',
    'K708'
  ),
  ordered = T
)

ggboxplot(
  wdv_titer[wdv_titer$replicate != 'rep3', ],
  x = 'gene',
  y = 'titer',
  add = 'point',
  add.params = list(color = 'grey50',
                    shape = 1)
) +
  rotate_x_text()
ggsave(paste0(pd,'wdv_crude_lysate_titer.pdf'), device = 'pdf')

wdv_inf = read.csv(
  './input/wdv_variant_inf.csv',
  colClasses = c('factor',
                 'factor',
                 'factor',
                 'numeric')
)
wdv_inf$gene = factor(
  wdv_inf$gene,
  levels = c(
    'wtDJ',
    'S244',
    'L256',
    'S268',
    'K311',
    'T331',
    'E349',
    'T414',
    'N664',
    'Y702',
    'K708'
  ),
  ordered = T
)

ggline(
  wdv_inf,
  x = 'moi',
  y = 'tdt_pos',
  color = 'gene',
  add = c('mean_se', 'point'),
  add.params = list(size = 0.5),
  facet.by = 'gene',
  nrow = 1,
) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(paste0(pd,'wdv_crude_lysate_inf.pdf'), device = 'pdf')
