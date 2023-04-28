##### look at luciferase SNPs -------
library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)

#### read data -----
### read differential results ####
CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.CM.new_back.005.txt', sep='\t', header = T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.VSMC.new_back.005.txt', sep='\t', header = T)

### index 
index <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/Hypertension__pooled.index.txt', sep = '\t', header = T)
head(index)

#### merge results with index ####
index['dupe_info'] <- gsub('\\..*','',index$tile_id)

CM_results = merge(CM_results,index[index$tile_type != 'RANDOM',c("name","dupe_info")], by="dupe_info")
VSMC_results = merge(VSMC_results,index[index$tile_type != 'RANDOM',c("name","dupe_info")], by="dupe_info")

CM_results <- CM_results %>% distinct()
VSMC_results <- VSMC_results %>% distinct()

CM_results
VSMC_results

CM_results <- CM_results %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(CM_results)
VSMC_results <- VSMC_results %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(VSMC_results)

### luciferase snps ----
luciferase <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/luciferase_snps.txt',
                         sep = '\t', header = F)
colnames(luciferase) <- c('snp_info','sentinel','gene')

##### merge info ----
luci_stats <- merge(luciferase, CM_results[,c('fdr_comp','logFC_comp','snp_info')])
luci_stats <- merge(luci_stats, VSMC_results[,c('fdr_comp','logFC_comp','snp_info')], by = 'snp_info')

### swipe logfc of reversed snps -----
luci_stats[luci_stats$snp_info == 'rs319686',c('logFC_comp.x','logFC_comp.y')] <- -luci_stats[luci_stats$snp_info == 'rs319686',c('logFC_comp.x','logFC_comp.y')]
luci_stats[luci_stats$snp_info == 'rs6773208',c('logFC_comp.x','logFC_comp.y')] <- -luci_stats[luci_stats$snp_info == 'rs319686',c('logFC_comp.x','logFC_comp.y')]
luci_stats[luci_stats$snp_info == 'rs11712464',c('logFC_comp.x','logFC_comp.y')] <- -luci_stats[luci_stats$snp_info == 'rs319686',c('logFC_comp.x','logFC_comp.y')]
luci_stats[luci_stats$snp_info == 'rs12330363',c('logFC_comp.x','logFC_comp.y')] <- -luci_stats[luci_stats$snp_info == 'rs319686',c('logFC_comp.x','logFC_comp.y')]

## bubble 
# cm 662b8fff
# vsmc f05928ff
luci_stats <- luci_stats[order(luci_stats$gene),]
colnames(luci_stats) <- c('snp_info','sentinel','gene', 'CM_fdr_comp','CM_logFC_comp', 'VSMC_fdr_comp','VSMC_logFC_comp')
library(reshape2)
luci_stats_fdr <- luci_stats[,c('snp_info','sentinel','gene', 'CM_fdr_comp','VSMC_fdr_comp')]
luci_stats_logfc <- luci_stats[,c('snp_info','sentinel','gene', 'CM_logFC_comp','VSMC_logFC_comp')]

colnames(luci_stats_fdr) <- c('snp_info','sentinel','gene','CM','VSMC')
colnames(luci_stats_logfc) <- c('snp_info','sentinel','gene','CM','VSMC')

luci_stats_fdr <- melt(luci_stats_fdr, id.vars = c('snp_info','sentinel','gene'))
luci_stats_logfc <- melt(luci_stats_logfc, id.vars = c('snp_info','sentinel','gene'))

colnames(luci_stats_fdr) <- c('snp_info','sentinel','gene','cell','FDR')
colnames(luci_stats_logfc) <- c('snp_info','sentinel','gene','cell','logFC')

all <- merge(luci_stats_fdr, luci_stats_logfc)
all <- all[order(all$gene, decreasing = T),]

library(scales)
ggplot(all, aes(y=snp_info, x=cell, size=-log10(FDR), fill=logFC)) +
  geom_point(alpha=0.8, shape=21, color="black") +
  scale_size(range = c(0.5, 12), name="FDR", limits = c(0.05,103)) +
  #scale_fill_gradient2(low=muted("dark blue"), high=muted("red"), name='logFC', breaks = c(-2,0,2))+
  #theme_ipsum() +
  scale_fill_gradientn(values=c(1, .63, .53, .43, 0), colours=c(muted("dark blue"), muted("dark blue"), "white", muted("red"), muted("red")))+
  theme(legend.position="right") +
  ylab("") +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        panel.grid = element_line(colour = 'light grey'))

library(forcats)
all %>%
  mutate(snp_info = fct_reorder(snp_info, desc(gene))) %>%
  ggplot( aes(y=snp_info, x=cell, size=-log10(FDR), fill=logFC)) +
  geom_point(alpha=0.8, shape=21, color="black") +
  scale_size(range = c(0.5, 12), name="FDR", limits = c(0.05,103)) +
  #scale_fill_gradient2(low=muted("dark blue"), high=muted("red"), name='logFC', breaks = c(-2,0,2))+
  #theme_ipsum() +
  scale_fill_gradientn(values=c(1, .63, .53, .43, 0), colours=c(muted("dark blue"), muted("dark blue"), "white", muted("red"), muted("red")))+
  theme(legend.position="right") +
  ylab("") +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        panel.grid = element_line(colour = 'light grey'))

### test significance ----
# one-proportion test
test <- prop.test(
  x = 44, # number of successes
  n = 62, # total number of trials
  p = 0.5 # we test for equal proportion so prob = 0.5 in each group
)

test

### All regulatory ranked by clustering -----
sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

#### read length block -----
length <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_snps_length_block.txt',
                     sep = '\t', header = T)
#1,000,000

number_snps_all <- as.data.frame(table(sentinel_snps$sentinel))
number_snps_all <- merge(number_snps_all, length, by.x='Var1', by.y='sentinel')
number_snps_all$clust <- number_snps_all$Freq/(number_snps_all$length/1000000)

CM_results <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/CM_eqtls_rep_cons_dhs_tbs_scores.txt')
VSMC_results <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/VSMC_eqtls_rep_cons_dhs_tbs_scores.txt')
  
#merge results with clustering ------
CM_results <- merge(CM_results, sentinel_snps, by.x = 'snp_info', by.y = 'snp')
VSMC_results <- merge(VSMC_results, sentinel_snps, by.x = 'snp_info', by.y = 'snp')

number_snps_CM <- merge(CM_results, number_snps_all, by.x='sentinel', by.y='Var1')
number_snps_VSMC <- merge(VSMC_results, number_snps_all, by.x='sentinel', by.y='Var1')

number_snps_CM <- number_snps_CM[number_snps_CM$fdr_comp<0.05,]
number_snps_VSMC <- number_snps_VSMC[number_snps_VSMC$fdr_comp<0.05,]

##### merge info ----
stats <- merge(number_snps_CM[,c('fdr_comp','logFC_comp','snp_info','clust')], number_snps_VSMC[,c('fdr_comp','logFC_comp','snp_info','clust')], by = c('snp_info', 'clust'), all.x = T, all.y = T)

## bubble 
# cm 662b8fff
# vsmc f05928ff
colnames(stats) <- c('snp_info','clust', 'CM_fdr_comp','CM_logFC_comp', 'VSMC_fdr_comp','VSMC_logFC_comp')
library(reshape2)
luci_stats_fdr <- stats[,c('snp_info','clust', 'CM_fdr_comp','VSMC_fdr_comp')]
luci_stats_logfc <- stats[,c('snp_info','clust', 'CM_logFC_comp','VSMC_logFC_comp')]

colnames(luci_stats_fdr) <- c('snp_info','clust','CM','VSMC')
colnames(luci_stats_logfc) <- c('snp_info','clust','CM','VSMC')

luci_stats_fdr <- melt(luci_stats_fdr, id.vars = c('snp_info','clust'))
luci_stats_logfc <- melt(luci_stats_logfc, id.vars = c('snp_info','clust'))

colnames(luci_stats_fdr) <- c('snp_info','clust','cell','FDR')
colnames(luci_stats_logfc) <- c('snp_info','clust','cell','logFC')

all <- merge(luci_stats_fdr, luci_stats_logfc)

library(forcats)
all %>%
  mutate(snp_info = fct_reorder(snp_info, desc(clust))) %>%
  ggplot( aes(y=snp_info, x=cell, fill=logFC)) +
  geom_point(alpha=0.8, shape=21, color="black") +
  #scale_size(range = c(0.5, 12), name="FDR", limits = c(0.05,103)) +
  #scale_fill_gradient2(low=muted("dark blue"), high=muted("red"), name='logFC', breaks = c(-2,0,2))+
  #theme_ipsum() +
  scale_fill_gradientn(values=c(1, .63, .53, .43, 0), colours=c(muted("red"), muted("red"), "white", muted("dark blue"), muted("dark blue")))+
  theme(legend.position="right") +
  ylab("") +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        panel.grid = element_line(colour = 'light grey'), axis.text.y = element_blank(),
        axis.ticks = element_blank())

sum(VSMC_results$logFC_comp>0)
nrow(VSMC_results)
sum(CM_results$logFC_comp>0)
nrow(CM_results)

test <- prop.test(
  x = 866, # number of successes 1172
  n = 1788, # total number of trials 2179
  p = 0.5 # we test for equal proportion so prob = 0.5 in each group
)

test

binom.test(866, 1788, 0.5, alternative="greater")

stats <- stats %>%
  mutate(snp_info = fct_reorder(snp_info, desc(clust)))  
stats <- stats[order(stats$clust, decreasing = T),]

stats_filt <- stats[c(1:500),]

pheatmap::pheatmap(as.matrix(stats_filt[,c('CM_logFC_comp','VSMC_logFC_comp')]), 
                     scale = 'none', cluster_rows = F, na_col = 'grey', 
                   color = c(("red"), ("blue")), show_rownames = FALSE, 
                   breaks = c(-4, 0, 4))

#### TFBS on luciferase SNPs -----
### read data 
### read fimo results ####
tile_motifs = read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/tiles_tfbs.best_motifs.csv', sep='\t', header = T)

head(tile_motifs)
tile_motifs <- tile_motifs %>%
  separate(sequence_name, c("pos", "snp_info","info"), "__")
head(tile_motifs)
motifs_wildtype <- tile_motifs[tile_motifs$tile_type == 'WILDTYPE_BUT',]
motifs_snps <- tile_motifs[tile_motifs$tile_type == 'WILDTYPE_SNP',]

motifs_luciferase_wt <- motifs_wildtype[motifs_wildtype$snp_info %in% luciferase$snp_info,]
motifs_luciferase_alt <- motifs_snps[motifs_snps$snp_info %in% luciferase$snp_info,]

library(ComplexHeatmap)
lt = list(CPEB4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'CPEB4'],"motif_alt_id"]),
          ESR1_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ESR1'],"motif_alt_id"]),
          INSR_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'INSR'],"motif_alt_id"]),
          MAP4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'MAP4'],"motif_alt_id"]),
          PDE5A_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'PDE5A'],"motif_alt_id"]),
          SMARCC1_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'SMARCC1'],"motif_alt_id"]),
          ST7L_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ST7L'],"motif_alt_id"]),
          ULK4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ULK4'],"motif_alt_id"]),
          CPEB4_alt = unique(motifs_luciferase_alt[motifs_luciferase_alt$snp_info %in% luciferase$snp_info[luciferase$gene == 'CPEB4'],"motif_alt_id"]),
          ESR1_alt = unique(motifs_luciferase_alt[motifs_luciferase_alt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ESR1'],"motif_alt_id"]),
          INSR_alt = unique(motifs_luciferase_alt[motifs_luciferase_alt$snp_info %in% luciferase$snp_info[luciferase$gene == 'INSR'],"motif_alt_id"]),
          MAP4_alt = unique(motifs_luciferase_alt[motifs_luciferase_alt$snp_info %in% luciferase$snp_info[luciferase$gene == 'MAP4'],"motif_alt_id"]),
          PDE5A_alt = unique(motifs_luciferase_alt[motifs_luciferase_alt$snp_info %in% luciferase$snp_info[luciferase$gene == 'PDE5A'],"motif_alt_id"]),
          SMARCC1_alt = unique(motifs_luciferase_alt[motifs_luciferase_alt$snp_info %in% luciferase$snp_info[luciferase$gene == 'SMARCC1'],"motif_alt_id"]),
          ST7L_alt = unique(motifs_luciferase_alt[motifs_luciferase_alt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ST7L'],"motif_alt_id"]),
          ULK4_alt = unique(motifs_luciferase_alt[motifs_luciferase_alt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ULK4'],"motif_alt_id"]))

lt = list(CPEB4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'CPEB4'],"motif_alt_id"]),
          ESR1_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ESR1'],"motif_alt_id"]),
          INSR_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'INSR'],"motif_alt_id"]),
          MAP4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'MAP4'],"motif_alt_id"]),
          PDE5A_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'PDE5A'],"motif_alt_id"]),
          SMARCC1_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'SMARCC1'],"motif_alt_id"]),
          ST7L_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ST7L'],"motif_alt_id"]),
          ULK4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ULK4'],"motif_alt_id"]))

lt = list(
          MAP4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'MAP4'],"motif_alt_id"]),
          ULK4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ULK4'],"motif_alt_id"]),
          MAP4_alt = unique(motifs_luciferase_alt[motifs_luciferase_alt$snp_info %in% luciferase$snp_info[luciferase$gene == 'MAP4'],"motif_alt_id"]),
          ULK4_alt = unique(motifs_luciferase_alt[motifs_luciferase_alt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ULK4'],"motif_alt_id"]),
          unique = wt)

m = make_comb_mat(lt)
UpSet(m)

ht = draw(UpSet(m, 
                comb_order = order(comb_size(m), decreasing = T),
                right_annotation = upset_right_annotation(m,
                                                          add_numbers = TRUE),
                top_annotation = upset_top_annotation(m, add_numbers = TRUE)))


extract_comb(m, "10010001")
wt <- extract_comb(m, "00010001")

#### compare single luciferase FIMO with haplotype ------
### New TFBS created?
tfbs_luci <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/luciferase/mapped_best_motifs_luciferase.txt',
                        sep = '\t', header = F)
head(tfbs_luci)
# 135 bp
# ULK = 10 snps
# CPEB4 = 3
# ESR1 = 2
# INSR = 1
# MAP4 = 6
# PDE5A = 2
# SMARCC1 = 4
# ST7L = 2

ulk <- list()
for (snp in c(1:10)) {
  tfbs <- tfbs_luci[tfbs_luci$V4<(135*snp) & tfbs_luci$V5>(135*snp) & tfbs_luci$V3 == 'ULK',]
  ulk[[snp]] <- tfbs
}
cpeb4 = list()
for (snp in c(1:3)) {
  tfbs <- tfbs_luci[tfbs_luci$V4<(135*snp) & tfbs_luci$V5>(135*snp) & tfbs_luci$V3 == 'CPEB4',]
  cpeb4[[snp]] <- tfbs
}
ESR1 = list()
for (snp in c(1:2)) {
  tfbs <- tfbs_luci[tfbs_luci$V4<(135*snp) & tfbs_luci$V5>(135*snp) & tfbs_luci$V3 == 'ESR1',]
  ESR1[[snp]] <- tfbs
}
MAP4 = list()
for (snp in c(1:6)) {
  tfbs <- tfbs_luci[tfbs_luci$V4<(135*snp) & tfbs_luci$V5>(135*snp) & tfbs_luci$V3 == 'MAP4',]
  MAP4[[snp]] <- tfbs
}
PDE5A = list()
for (snp in c(1:2)) {
  tfbs <- tfbs_luci[tfbs_luci$V4<(135*snp) & tfbs_luci$V5>(135*snp) & tfbs_luci$V3 == 'PDE5A',]
  PDE5A[[snp]] <- tfbs
}
SMARCC1 = list()
for (snp in c(1:4)) {
  tfbs <- tfbs_luci[tfbs_luci$V4<(135*snp) & tfbs_luci$V5>(135*snp) & tfbs_luci$V3 == 'SMARCC1',]
  SMARCC1[[snp]] <- tfbs
}
ST7L = list()
for (snp in c(1:2)) {
  tfbs <- tfbs_luci[tfbs_luci$V4<(135*snp) & tfbs_luci$V5>(135*snp) & tfbs_luci$V3 == 'ST7L',]
  ST7L[[snp]] <- tfbs
}

## now look at unique TF ----
lt = list(CPEB4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'CPEB4'],"motif_alt_id"]),
          ESR1_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ESR1'],"motif_alt_id"]),
          INSR_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'INSR'],"motif_alt_id"]),
          MAP4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'MAP4'],"motif_alt_id"]),
          PDE5A_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'PDE5A'],"motif_alt_id"]),
          SMARCC1_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'SMARCC1'],"motif_alt_id"]),
          ST7L_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ST7L'],"motif_alt_id"]),
          ULK4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ULK4'],"motif_alt_id"]))


lt_oligos_luci = list(CPEB4_luc = unique(tfbs_luci[tfbs_luci$V3 == 'CPEB4',"V2"]),
                      ESR1_luc = unique(tfbs_luci[tfbs_luci$V3 == 'ESR1',"V2"]),
                      INSR_luc = unique(tfbs_luci[tfbs_luci$V3 == 'INSR',"V2"]),
                      MAP4_luc = unique(tfbs_luci[tfbs_luci$V3 == 'MAP4',"V2"]),
                      PDE5A_luc = unique(tfbs_luci[tfbs_luci$V3 == 'PDE5A',"V2"]),
                      SMARCC1_luc = unique(tfbs_luci[tfbs_luci$V3 == 'SMARCC1',"V2"]),
                      ST7L_luc = unique(tfbs_luci[tfbs_luci$V3 == 'ST7L',"V2"]),
                      ULK4_luc = unique(tfbs_luci[tfbs_luci$V3 == 'ULK4',"V2"]),
                      CPEB4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'CPEB4'],"motif_alt_id"]),
                      ESR1_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ESR1'],"motif_alt_id"]),
                      INSR_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'INSR'],"motif_alt_id"]),
                      MAP4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'MAP4'],"motif_alt_id"]),
                      PDE5A_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'PDE5A'],"motif_alt_id"]),
                      SMARCC1_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'SMARCC1'],"motif_alt_id"]),
                      ST7L_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ST7L'],"motif_alt_id"]),
                      ULK4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ULK4'],"motif_alt_id"]))


lt_oligos_luci = list(ULK4_luc = unique(tfbs_luci[tfbs_luci$V3 == 'ULK4',"V2"]),
                      ULK4_wt = unique(motifs_luciferase_wt[motifs_luciferase_wt$snp_info %in% luciferase$snp_info[luciferase$gene == 'ULK4'],"motif_alt_id"]))

m = make_comb_mat(lt_oligos_luci)
UpSet(m)

ht = draw(UpSet(m, 
                comb_order = order(comb_size(m), decreasing = T),
                right_annotation = upset_right_annotation(m,
                                                          add_numbers = TRUE),
                top_annotation = upset_top_annotation(m, add_numbers = TRUE)))


extract_comb(m, "10")
wt <- extract_comb(m, "00010001")

### plot luciferase counts -----

counts_dir = "~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/01_counts/"
mpranalyze_dir = paste0(counts_dir,"/mpranalyze_files")

### VSMC data ####
VSMC_data = read.table(paste0(counts_dir,"/VSMC__all_counts_final.4.txt"), sep="\t", header = T)
head(VSMC_data)

### CM data ####
CM_data = read.table(paste0(counts_dir,"/CM__all_counts_final.5.txt"), sep="\t", header = T)
head(CM_data)

colnames(VSMC_data) = c("barcode", "dna_1", "VSMC_rep1", "VSMC_rep2", "VSMC_rep3","VSMC_rep6")
colnames(CM_data) = c("barcode", "dna_1", "CM_rep1", "CM_rep2", "CM_rep3", "CM_rep4","CM_rep5")

df_VSMC = merge(VSMC_data, index[,c("barcode", "element", "tile_type", "tile_id","snp","name","dupe_info")], by="barcode")
df_CM = merge(CM_data, index[,c("barcode", "element", "tile_type", "tile_id","snp","name","dupe_info")], by="barcode")
#df["barc_id"] = df.apply(get_barc_id, axis=1).astype(int)
head(df_VSMC)

#positive_controls_VSMC <- df_VSMC[df_VSMC$tile_type %in% c('CONTROL_BUT_HAS_SNP','CONTROL_SNP_INDIV'),]
#positive_controls_CM <- df_CM[df_CM$tile_type %in% c('CONTROL_BUT_HAS_SNP','CONTROL_SNP_INDIV'),]

df_CM$median_rep = rowMeans(df_CM[,c("CM_rep1", "CM_rep2", "CM_rep3", "CM_rep4","CM_rep5")], na.rm = T)
df_VSMC$median_rep = rowMeans(df_VSMC[,c("VSMC_rep1", "VSMC_rep2", "VSMC_rep3","VSMC_rep6")], na.rm = T)

dupe_snp <- CM_results[CM_results$snp_info %in% luciferase$snp_info,c('dupe_info','snp_info')]
dupe_snp <- dupe_snp[dupe_snp$snp!='none',]
dupe_snp <- dupe_snp %>% distinct()

positive_controls_VSMC <- df_VSMC[df_VSMC$dupe_info %in% dupe_snp$dupe_info,]
positive_controls_CM <- df_CM[df_CM$dupe_info %in% dupe_snp$dupe_info,]

positive_controls_CM$activity = positive_controls_CM$median_rep/positive_controls_CM$dna_1
positive_controls_VSMC$activity = positive_controls_VSMC$median_rep/positive_controls_VSMC$dna_1

ggplot(positive_controls_CM, aes(x = dupe_info, y = activity, fill=as.factor(tile_type))) +
  #geom_jitter(aes(color=tile_type),alpha=0.4)+
  geom_point(aes(color=as.factor(tile_type)),pch = 21, alpha=0.3, position = position_jitterdodge())+
  geom_boxplot(outlier.shape = NA)+
  ylab("median counts") + 
  xlab("") + 
  theme(legend.position="top", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(color="", fill = "")+
  scale_color_discrete(name = "", labels = c("REF", "ALT"))+
  scale_fill_discrete(name = "", labels = c("REF", "ALT"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  #scale_color_brewer(palette="Set1") + 
  scale_x_discrete(breaks = dupe_snp$dupe_info, labels = dupe_snp$snp_info)+
  #stat_compare_means(aes(group = tile_type), label = "p.signif", hide.ns = T)+
  ylim(c(0,5))

ggplot(positive_controls_CM, aes(x = tile_type, y = activity)) +
  #geom_jitter(aes(color=tile_type),alpha=0.4)+
  geom_point(aes(color=tile_type),pch = 21, alpha=0.3, position = position_jitterdodge())+
  geom_violin(outlier.shape = NA, alpha=0.5)+
  ylab("Ratio RNA/DNA") + 
  xlab("") + 
  scale_color_discrete(name = "", labels = c("REF", "ALT"))+
  scale_fill_discrete(name = "", labels = c("REF", "ALT"))+
  theme(legend.position="top", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))
  ylim(0,10)
#scale_color_brewer(palette="Set1") + 
  
  
