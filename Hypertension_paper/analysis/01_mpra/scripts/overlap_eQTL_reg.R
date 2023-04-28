#### Overlap GTEx eQTLs all tissues with regulatory variants #####
### GTEx eQTLs ####
### these data comes from GTEx portal ##
library(dplyr)
library(tidyr)

gtex_path <- '~/marenostrum/Projects/GTEx_v8/Manuscript/'
TissueInfoFile <- "TissuesInfo.rds" 

TissueInfo <- readRDS(paste0(gtex_path,TissueInfoFile))
tissues <- as.character(TissueInfo$Tissue_id)
##  kidneys, ureters, bladder, and the urethra.
## Nerve Brain
## pancreas, pituitary, thyroid, Glands, Ovary, Testis
subset <- tissues[grep('Arter|Heart|Brain|Pitu|Gland|Thyr|Panc|Nerve|Kidn|Uret|Blad',tissues)]

eQTLs_path <- '~/marenostrum_scratch/GTEx/v8/cis_QTLs/cis_eQTLs/GTEx_Analysis_v8_eQTL/'

eQTLs <- lapply(subset, function(tissue) read.table(paste0(eQTLs_path,tissue,".v8.signif_variant_gene_pairs.txt.gz"), sep='\t', header=T) )
names(eQTLs) <- subset

#whole_blood <- 

eQTLs <- lapply(subset, function(tissue) eQTLs[[tissue]] %>%
                  separate(variant_id, c("chrom", "start","ref","alt",'genome'), "_"))
names(eQTLs) <- subset
eQTLs_n <- lapply(subset, function(tissue) eQTLs[[tissue]]$variant_id <- paste0(eQTLs[[tissue]]$chrom,
                                                                                '_',eQTLs[[tissue]]$start,
                                                                                '_',eQTLs[[tissue]]$ref,
                                                                                '_',eQTLs[[tissue]]$alt,
                                                                                '_',eQTLs[[tissue]]$genome))
names(eQTLs_n) <- subset

liftover_coord <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/snp_coord_master.hg38.simple.chr.sorted.bed', sep='\t', header=F)
colnames(liftover_coord) <- c('chr_b38', 'start_b38', 'end_b38','info_snp')

sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]
liftover_coord <- merge(liftover_coord, sentinel_snps, by.x='info_snp', by.y='change')
liftover_coord <- liftover_coord %>% distinct()

liftover_coord$variant_id <- paste0(liftover_coord$chr_b38,'_',as.numeric(liftover_coord$start_b38)+1,'_',gsub(':.*','',liftover_coord$E),'_',gsub('.*:','',liftover_coord$E),'_b38')
head(liftover_coord)

## slope > 0 --> ALT more active
## slope < 0 --> REF more active 
liftover_coord$start_b38 <-  liftover_coord$start_b38+1
eQTLs_mpra <- lapply(subset, function(tissue) as.data.frame(merge(liftover_coord, eQTLs[[tissue]][,c("chrom", "start","ref","alt","slope","gene_id")], 
                                                                  by.x=c('chr_b38', 'start_b38'), 
                                                                  by.y=c("chrom", "start"))))
names(eQTLs_mpra) <- subset
eQTLs_mpra_S <- lapply(subset, function(tissue) as.data.frame(eQTLs_mpra[[tissue]][eQTLs_mpra[[tissue]]$variant_id %in% eQTLs_n[[tissue]],]))
names(eQTLs_mpra_S) <- subset


########################
### GTEx eQTLs mapping to MPRA SNPs #####

#### MPRA data ####
# CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)
# VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)

CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/all_snps_info_risk_CM.txt', sep='\t', header=T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/all_snps_info_risk_VSMC.txt', sep='\t', header=T)

VSMC_results$info_snp <- VSMC_results$change
CM_results$info_snp <- CM_results$change

VSMC_results <- merge(VSMC_results, liftover_coord, all.x = TRUE)
CM_results <- merge(CM_results, liftover_coord, all.x = TRUE)

### merge with eqtl info ####
VSMC_results$variant_id <- paste0(VSMC_results$chr_b38,'_',as.numeric(VSMC_results$start_b38)+1,'_',gsub(':.*','',VSMC_results$E),'_',gsub('.*:','',VSMC_results$E),'_b38')
CM_results$variant_id <- paste0(CM_results$chr_b38,'_',as.numeric(CM_results$start_b38)+1,'_',gsub(':.*','',CM_results$E),'_',gsub('.*:','',CM_results$E),'_b38')


eQTLs_mpra_CM <- lapply(subset, function(tissue) as.data.frame(CM_results[CM_results$info_snp %in% eQTLs_mpra_S[[tissue]]$info_snp,]))
names(eQTLs_mpra_CM) <- subset

eQTLs_mpra_VSMC <- lapply(subset, function(tissue) as.data.frame(VSMC_results[VSMC_results$info_snp %in% eQTLs_mpra_S[[tissue]]$info_snp,]))
names(eQTLs_mpra_VSMC) <- subset

eQTLs_mpra_CM_df <- do.call(rbind,eQTLs_mpra_CM)
eQTLs_mpra_VSMC_df <- do.call(rbind,eQTLs_mpra_VSMC)

# colnames(eQTLs_mpra_CM_df) <- c('variant_id')
# colnames(eQTLs_mpra_VSMC_df) <- c('variant_id')

eQTLs_mpra_CM_df$tissue <- gsub('\\..*','',rownames(eQTLs_mpra_CM_df))
eQTLs_mpra_VSMC_df$tissue <- gsub('\\..*','',rownames(eQTLs_mpra_VSMC_df))

library(RColorBrewer)
n <- 24
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors = sample(col_vector,n)

n_eqtls <- lapply(subset, function(tissue) as.data.frame(nrow(eQTLs[[tissue]])))
names(n_eqtls) <- subset
n_eqtls_df <- do.call(rbind,n_eqtls)
colnames(n_eqtls_df) <- c('n_eqtls')
n_eqtls_df$tissue <- rownames(n_eqtls_df)

eQTLs_mpra_CM_df_funct <- eQTLs_mpra_CM_df[eQTLs_mpra_CM_df$Functional == 'Functional Variant',]
eQTLs_mpra_VSMC_df_funct <- eQTLs_mpra_VSMC_df[eQTLs_mpra_VSMC_df$Functional == 'Functional Variant',]

df2 <- eQTLs_mpra_VSMC_df_funct %>% group_by(tissue) %>%   mutate(count_name_occurr = n())
df2 <- merge(df2, n_eqtls_df)

ggplot(df2, aes(x=reorder(tissue,-n_eqtls), fill= tissue)) +
  geom_bar(stat = 'count', fill = colors, colour='black') + 
  xlab('')+
  coord_flip()+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') +
  ylab('Number of Functional SNPs')

df2$perc <-  df2$count_name_occurr/391
perc <- df2[,c('count_name_occurr','n_eqtls','perc','tissue')] %>% distinct()
## colors CMs <- (#DBAFD7, #673597) VSMCs <- (#FECFA7, #F15A3F)
ggplot(perc, aes(x=reorder(tissue,-n_eqtls), fill= perc, y=perc)) +
  geom_bar(stat = 'identity', colour='black') + 
  scale_fill_gradient(low="#FECFA7",high="#F15A3F")+
  xlab('')+
  coord_flip()+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') +
  ylab('Perc of regulatory variants \n as eQTLs')

ggplot(perc, aes(x=reorder(tissue,-n_eqtls), y=n_eqtls)) +
  geom_bar(stat = 'identity', colour='black', fill='dark grey') + xlab('')+
  coord_flip() +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
                     legend.position = 'none') +
  ylab('Nº of eQTLs')

