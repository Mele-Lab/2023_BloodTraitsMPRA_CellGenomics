#### analyze results from MPRA comp mode
CM_new <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.CM.new_back.005.txt', sep='\t', header = T)
VSMC_new <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.VSMC.new_back.005.txt', sep='\t', header = T)

### get diff and active snps ####
vals_significance <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/CM_vals.significance.5.txt', sep = '\t', header = T)
head(vals_significance)
vals_significance <- vals_significance %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance)

vals_significance[vals_significance$type == 'CONTROL_SNP_INDIV','group'] <- 'Control ALT'
vals_significance[vals_significance$type == 'CONTROL_BUT_HAS_SNP','group'] <- 'Control REF'
vals_significance[vals_significance$type %in% c('WILDTYPE_SNP_INDIV','WILDTYPE_BUT_HAS_SNP') & vals_significance$sig == 'sig', 'group'] <- 'Active element'
vals_significance[vals_significance$type %in% c('WILDTYPE_SNP_INDIV','WILDTYPE_BUT_HAS_SNP') & vals_significance$sig == 'not sig', 'group'] <- 'Not Active element'
vals_significance[vals_significance$type == 'negative control','group'] <- 'Random Sequence'

active_snps <- unique(vals_significance[vals_significance$group == 'Active element','snp_info'])
active_snps <- unique(vals_significance[vals_significance$group == 'Active element','dupe_info'])

snps <- vals_significance[vals_significance$tile_type %in% c('WILDTYPE_BUT_HAS_SNP','WILDTYPE_SNP_INDIV'),c('pos','snp_info','chrom','dupe_info')] %>% distinct()

head(CM_new)
CM_new_active_diff <- CM_new[CM_new$fdr_comp < 0.05 & CM_new$dupe_info %in% active_snps,]
CM_new_active_diff <- merge(CM_new_active_diff, snps, by='dupe_info')
CM_new_active_diff <- CM_new_active_diff[order(abs(CM_new_active_diff$logFC_comp), decreasing = T),]


head(VSMC_new)
VSMC_new_active_diff <- VSMC_new[VSMC_new$fdr_comp < 0.05 & VSMC_new$dupe_info %in% active_snps,]
VSMC_new_active_diff <- merge(VSMC_new_active_diff, snps, by='dupe_info')
VSMC_new_active_diff <- VSMC_new_active_diff[order(abs(VSMC_new_active_diff$logFC_comp), decreasing = T),]

overlap <- CM_new_diff[CM_new_diff$dupe_info %in% VSMC_new_diff$dupe_info,]
overlap <- merge(overlap, snps, by='dupe_info')
overlap <- overlap[order(abs(overlap$logFC_comp), decreasing = T),]

write.table(CM_new_active_diff, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_active_snps_final.txt', quote=F, row.names=F, sep='\t')
write.table(VSMC_new_active_diff, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_active_snps_final.txt', quote=F, row.names=F, sep='\t')
write.table(overlap, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/Overlap_diff_active_snps_final.txt', quote=F, row.names=F, sep='\t')

#### plot distribution of FC ####
CM_new$group <- 'NonFunctionalVariant'
CM_new$group[CM_new$fdr_comp < 0.05] <- 'FunctionalVariant'
CM_new$group[CM_new$fdr_comp < 0.05 & CM_new$dupe_info %in% active_snps] <- 'Functional+ActiveVariant'

CM_new_test <- CM_new[CM_new$tile_type %in% c('WILDTYPE_SNP_INDIV','WILDTYPE_BUT_HAS_SNP'),]

VSMC_new$group <- 'NonFunctionalVariant'
VSMC_new$group[VSMC_new$fdr_comp < 0.05] <- 'FunctionalVariant'
VSMC_new$group[VSMC_new$fdr_comp < 0.05 & VSMC_new$dupe_info %in% active_snps] <- 'Functional+ActiveVariant'

VSMC_new_test <- VSMC_new[VSMC_new$tile_type %in% c('WILDTYPE_SNP_INDIV','WILDTYPE_BUT_HAS_SNP'),]

ggplot(CM_new_test, aes(x=group,y=abs(logFC_comp), color=group)) + geom_violin(fill='white') +
  geom_jitter(aes(fill=group))+
  ylab('abs(log2FC)')+
  xlab('')+
  ggtitle('')+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_brewer(palette ='Set2') +
  scale_color_brewer(palette ='Set2')

# sample size
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
sample_size = CM_new_test %>% group_by(group) %>% dplyr::summarize(num=n())

# Plot
CM_new_test %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(group, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=abs(logFC_comp), fill=group)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("CM FoldChange") +
  xlab("")+
  ylim(c(0,5))

### number of functional variants per sentinel snp ####


# sample size
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
sample_size = VSMC_new %>% group_by(group) %>% dplyr::summarize(num=n())

# Plot
VSMC_new_test %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(group, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=abs(logFC_comp), fill=group)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("VSMC FoldChange") +
  xlab("")
#  ylim(c(0,5))

### Get ranking for all differential SNPs ####
CM_new_diff <- CM_new[CM_new$fdr_comp < 0.05,]
CM_new_diff <- CM_new_diff[order(abs(CM_new_diff$logFC_comp), decreasing = T),]
row.names(CM_new_diff) <- NULL
CM_new_diff$rankingFC <- rownames(CM_new_diff)

CM_new_diff <- CM_new_diff[order(abs(CM_new_diff$fdr_comp), decreasing = F),]
row.names(CM_new_diff) <- NULL
CM_new_diff$rankingFDR <- rownames(CM_new_diff)

VSMC_new_diff <- VSMC_new[VSMC_new$fdr_comp < 0.05,]
VSMC_new_diff <- VSMC_new_diff[order(abs(VSMC_new_diff$logFC_comp), decreasing = T),]
row.names(VSMC_new_diff) <- NULL
VSMC_new_diff$rankingFC <- rownames(VSMC_new_diff)

VSMC_new_diff <- VSMC_new_diff[order(abs(VSMC_new_diff$fdr_comp), decreasing = F),]
row.names(VSMC_new_diff) <- NULL
VSMC_new_diff$rankingFDR <- rownames(VSMC_new_diff)

CM_new_diff$active <- 'Not Active'
CM_new_diff$active[CM_new_diff$dupe_info %in% CM_new_active_diff$dupe_info] <- 'Active'
CM_new_diff$overlap <- 'Not Overlap'
CM_new_diff$overlap[CM_new_diff$dupe_info %in% overlap$dupe_info] <- 'Overlap'
CM_new_diff <- merge(CM_new_diff, snps, by='dupe_info')
CM_new_diff <- CM_new_diff[order(abs(CM_new_diff$fdr_comp), decreasing = F),]

VSMC_new_diff$active <- 'Not Active'
VSMC_new_diff$active[VSMC_new_diff$dupe_info %in% VSMC_new_active_diff$dupe_info] <- 'Active'
VSMC_new_diff$overlap <- 'Not Overlap'
VSMC_new_diff$overlap[VSMC_new_diff$dupe_info %in% overlap$dupe_info] <- 'Overlap'
VSMC_new_diff <- merge(VSMC_new_diff, snps, by='dupe_info')
VSMC_new_diff <- VSMC_new_diff[order(abs(VSMC_new_diff$fdr_comp), decreasing = F),]

write.table(CM_new_diff[CM_new_diff$active == 'Active',], '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_active_snps_final.ranking.txt', quote=F, row.names=F, sep='\t')
write.table(VSMC_new_diff[VSMC_new_diff$active == 'Active',], '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_active_snps_final.ranking.txt', quote=F, row.names=F, sep='\t')

### Get closest gene ####
closest_gene <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/clostest_gene_snps.txt', sep='\t', header=F)
head(CM_new_diff)
snps <- vals_significance[vals_significance$tile_type %in% c('WILDTYPE_BUT_HAS_SNP','WILDTYPE_SNP_INDIV'),c('pos','snp_info','chrom','dupe_info','info')] %>% distinct()
snps <- snps %>%
  separate(info, c("chr", "start","ref","alt"), ":")

CM_new_diff <- merge(CM_new_diff, snps)
CM_new_diff$start <- as.numeric(CM_new_diff$start)
CM_new_diff$end <- CM_new_diff$start + nchar(CM_new_diff$ref) - 1
CM_new_diff$pos <- NA
CM_new_diff$pos <- paste0(CM_new_diff$chrom,':',CM_new_diff$start,':',CM_new_diff$end)

closest_gene$pos <- paste0(closest_gene$V1,':',closest_gene$V2,':',closest_gene$V3)

CM_new_diff_genes <- merge(CM_new_diff, closest_gene[,c('V7','V8','pos')])
CM_new_diff_genes$dist <- gsub('[-+]\\|','',CM_new_diff_genes$V8)
CM_new_diff_genes$ensemblID <- gsub('\\..*','',CM_new_diff_genes$V7)

VSMC_new_diff <- merge(VSMC_new_diff, snps)
VSMC_new_diff$start <- as.numeric(VSMC_new_diff$start)
VSMC_new_diff$end <- VSMC_new_diff$start + nchar(VSMC_new_diff$ref) - 1
VSMC_new_diff$pos <- NA
VSMC_new_diff$pos <- paste0(VSMC_new_diff$chrom,':',VSMC_new_diff$start,':',VSMC_new_diff$end)

VSMC_new_diff_genes <- merge(VSMC_new_diff, closest_gene[,c('V7','V8','pos')])
VSMC_new_diff_genes$dist <- gsub('[-+]\\|','',VSMC_new_diff_genes$V8)
VSMC_new_diff_genes$ensemblID <- gsub('\\..*','',VSMC_new_diff_genes$V7)

## get gene name ####
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

genes_CM <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=(CM_new_diff_genes$ensemblID),mart= mart)
genes_VSMC <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=(VSMC_new_diff_genes$ensemblID),mart= mart)

VSMC_new_diff_genes <- merge(VSMC_new_diff_genes, genes_VSMC, by.x='ensemblID', by.y='ensembl_gene_id', all.x = TRUE)
CM_new_diff_genes <- merge(CM_new_diff_genes, genes_CM, by.x='ensemblID', by.y='ensembl_gene_id', all.x = TRUE)

write.table(CM_new_diff_genes[CM_new_diff_genes$active == 'Active',], '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_active_snps_final.ranking.txt', quote=F, row.names=F, sep='\t')
write.table(VSMC_new_diff_genes[VSMC_new_diff_genes$active == 'Active',], '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_active_snps_final.ranking.txt', quote=F, row.names=F, sep='\t')

## get eQTL info ####
### eQTL data downloaded from GTEx portal ######
gtex_path <- '~/marenostrum/Projects/GTEx_v8/Manuscript/'
TissueInfoFile <- "TissuesInfo.rds" 

TissueInfo <- readRDS(paste0(gtex_path,TissueInfoFile))
tissues <- as.character(TissueInfo$Tissue_id)
subset <- tissues[grep('Arter|Heart',tissues)]

eQTLs_path <- '~/marenostrum/Data/GTEx_v8/cisQTLs/eGenes/GTEx_Analysis_v8_eQTL/'

eQTLs <- lapply(subset, function(tissue) read.table(paste0(eQTLs_path,tissue,".v8.signif_variant_gene_pairs.txt.gz"), sep='\t', header=T) )
names(eQTLs) <- subset

eQTLs_coord <- lapply(subset, function(tissue) eQTLs[[tissue]]$variant_id )
eQTLs_coord <- unlist(eQTLs_coord)
eQTLs_coord <- as.data.frame(eQTLs_coord)
eQTLs_coord <- eQTLs_coord %>% distinct()
eQTLs_coord$chr <- gsub('_.*','',eQTLs_coord$eQTLs_coord)
eQTLs_coord <- eQTLs_coord %>%
  separate(eQTLs_coord, c("chr", "start","ref","alt","genome"), "_")

eQTLs_coord_list <- lapply(subset, function(tissue) eQTLs[[tissue]]$variant_id )
eQTLs_coord_list <- unlist(eQTLs_coord_list)
eQTLs_coord_list <- as.data.frame(eQTLs_coord_list)
eQTLs_coord_list <- eQTLs_coord_list %>% distinct()
eQTLs_coord$info <- eQTLs_coord_list$eQTLs_coord_list

eQTLs_coord$start <- as.numeric(eQTLs_coord$start) - 1 
eQTLs_coord$end <- eQTLs_coord$start + nchar(eQTLs_coord$ref) - 1

write.table(eQTLs_coord[,c('chrom','start','end','info')], '~/marenostrum/Data/Hypertension/GTEx_eQTLs.hg38.bed', quote = F, row.names = F, col.names = F, sep='\t')

#### Get FIMO data ####
# if N_diff_motifs < 0 it means more motifs on ALT than in REF sequence
fimo_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/03_fimo/diff_TFBS.csv', sep='\t', header = T)

VSMC_new_diff_genes_fimo <- merge(VSMC_new_diff_genes, fimo_results[c('SNP','TF_WT','TF_SNP','N_diff_motifs')], by.x='snp_info', by.y='SNP', all.x = TRUE)
CM_new_diff_genes_fimo <- merge(CM_new_diff_genes, fimo_results[c('SNP','TF_WT','TF_SNP','N_diff_motifs')], by.x='snp_info', by.y='SNP', all.x = TRUE)

VSMC_new_diff_genes_fimo <- VSMC_new_diff_genes_fimo[order((VSMC_new_diff_genes_fimo$fdr_comp), decreasing = F),]
CM_new_diff_genes_fimo <- CM_new_diff_genes_fimo[order((CM_new_diff_genes_fimo$fdr_comp), decreasing = F),]

## read liftover coord ####
liftover_coord <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/snp_coord_master.hg38.simple.chr.sorted.bed', sep='\t', header=F)
colnames(liftover_coord) <- c('chr_b38', 'start_b38', 'end_b38','info_snp')

VSMC_new_diff_genes_fimo$info_snp <- paste0(VSMC_new_diff_genes_fimo$chr,':',VSMC_new_diff_genes_fimo$start,':',VSMC_new_diff_genes_fimo$ref,':',VSMC_new_diff_genes_fimo$alt)
CM_new_diff_genes_fimo$info_snp <- paste0(CM_new_diff_genes_fimo$chr,':',CM_new_diff_genes_fimo$start,':',CM_new_diff_genes_fimo$ref,':',CM_new_diff_genes_fimo$alt)

VSMC_new_diff_genes_fimo <- merge(VSMC_new_diff_genes_fimo, liftover_coord, all.x = TRUE)
CM_new_diff_genes_fimo <- merge(CM_new_diff_genes_fimo, liftover_coord, all.x = TRUE)

### merge with eqtl info ####
head(eQTLs$Artery_Coronary)
VSMC_new_diff_genes_fimo$variant_id <- paste0(VSMC_new_diff_genes_fimo$chr_b38,'_',as.numeric(VSMC_new_diff_genes_fimo$start_b38)+1,'_',VSMC_new_diff_genes_fimo$ref,'_',VSMC_new_diff_genes_fimo$alt,'_b38')
VSMC_new_diff_genes_fimo$eQTL_Artery_Coronary[VSMC_new_diff_genes_fimo$variant_id %in% eQTLs$Artery_Coronary$variant_id] <- 'YES'
VSMC_new_diff_genes_fimo$eQTL_Heart_Atrial_Appendage[VSMC_new_diff_genes_fimo$variant_id %in% eQTLs$Heart_Atrial_Appendage$variant_id] <- 'YES'
VSMC_new_diff_genes_fimo$eQTL_Heart_Left_Ventricle[VSMC_new_diff_genes_fimo$variant_id %in% eQTLs$Heart_Left_Ventricle$variant_id] <- 'YES'
VSMC_new_diff_genes_fimo$eQTL_Artery_Aorta[VSMC_new_diff_genes_fimo$variant_id %in% eQTLs$Artery_Aorta$variant_id] <- 'YES'
VSMC_new_diff_genes_fimo$eQTL_Artery_Tibial[VSMC_new_diff_genes_fimo$variant_id %in% eQTLs$Artery_Tibial$variant_id] <- 'YES'

CM_new_diff_genes_fimo$variant_id <- paste0(CM_new_diff_genes_fimo$chr_b38,'_',as.numeric(CM_new_diff_genes_fimo$start_b38)+1,'_',CM_new_diff_genes_fimo$ref,'_',CM_new_diff_genes_fimo$alt,'_b38')
CM_new_diff_genes_fimo$eQTL_Artery_Coronary[CM_new_diff_genes_fimo$variant_id %in% eQTLs$Artery_Coronary$variant_id] <- 'YES'
CM_new_diff_genes_fimo$eQTL_Heart_Atrial_Appendage[CM_new_diff_genes_fimo$variant_id %in% eQTLs$Heart_Atrial_Appendage$variant_id] <- 'YES'
CM_new_diff_genes_fimo$eQTL_Heart_Left_Ventricle[CM_new_diff_genes_fimo$variant_id %in% eQTLs$Heart_Left_Ventricle$variant_id] <- 'YES'
CM_new_diff_genes_fimo$eQTL_Artery_Aorta[CM_new_diff_genes_fimo$variant_id %in% eQTLs$Artery_Aorta$variant_id] <- 'YES'
CM_new_diff_genes_fimo$eQTL_Artery_Tibial[CM_new_diff_genes_fimo$variant_id %in% eQTLs$Artery_Tibial$variant_id] <- 'YES'

VSMC_new_diff_genes_fimo <- VSMC_new_diff_genes_fimo[order((VSMC_new_diff_genes_fimo$fdr_comp), decreasing = F),]
CM_new_diff_genes_fimo <- CM_new_diff_genes_fimo[order((CM_new_diff_genes_fimo$fdr_comp), decreasing = F),]

write.table(CM_new_diff_genes_fimo, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_snps_final.ranking.withinfo.txt', quote=F, row.names=F, sep='\t')
write.table(VSMC_new_diff_genes_fimo, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_snps_final.ranking.withinfo.txt', quote=F, row.names=F, sep='\t')

eqtls <- merge(CM_new_diff_genes_fimo, eQTLs$Artery_Coronary[,c('variant_id','gene_id')], all.x=TRUE)
names_cols <- colnames(eqtls)[0:(ncol(eqtls)-1)]
colnames(eqtls) <- c(names_cols,'gene_eqtl_Artery_Coronary')
eqtls <- eqtls %>% 
  group_by_at(setdiff(names(eqtls), "gene_eqtl_Artery_Coronary")) %>%
  mutate(gene_eqtl_Artery_Coronary = paste0(gene_eqtl_Artery_Coronary, collapse = ",")) %>%
  ungroup() %>%
  distinct()

eqtls <- merge(eqtls, eQTLs$Heart_Atrial_Appendage[,c('variant_id','gene_id')], all.x=TRUE)
names_cols <- colnames(eqtls)[0:(ncol(eqtls)-1)]
colnames(eqtls) <- c(names_cols,'gene_eqtl_Heart_Atrial_Appendage')
eqtls <- eqtls %>% 
  group_by_at(setdiff(names(eqtls), "gene_eqtl_Heart_Atrial_Appendage")) %>%
  mutate(gene_eqtl_Heart_Atrial_Appendage = paste0(gene_eqtl_Heart_Atrial_Appendage, collapse = ",")) %>%
  ungroup() %>%
  distinct()

eqtls <- merge(eqtls, eQTLs$Heart_Left_Ventricle[,c('variant_id','gene_id')], all.x=TRUE)
names_cols <- colnames(eqtls)[0:(ncol(eqtls)-1)]
colnames(eqtls) <- c(names_cols,'gene_eqtl_Heart_Left_Ventricle')
eqtls <- eqtls %>% 
  group_by_at(setdiff(names(eqtls), "gene_eqtl_Heart_Left_Ventricle")) %>%
  mutate(gene_eqtl_Heart_Left_Ventricle = paste0(gene_eqtl_Heart_Left_Ventricle, collapse = ",")) %>%
  ungroup() %>%
  distinct()

eqtls <- merge(eqtls, eQTLs$Artery_Aorta[,c('variant_id','gene_id')], all.x=TRUE)
names_cols <- colnames(eqtls)[0:(ncol(eqtls)-1)]
colnames(eqtls) <- c(names_cols,'gene_eqtl_Artery_Aorta')
eqtls <- eqtls %>% 
  group_by_at(setdiff(names(eqtls), "gene_eqtl_Artery_Aorta")) %>%
  mutate(gene_eqtl_Artery_Aorta = paste0(gene_eqtl_Artery_Aorta, collapse = ",")) %>%
  ungroup() %>%
  distinct()

eqtls <- merge(eqtls, eQTLs$Artery_Tibial[,c('variant_id','gene_id')], all.x=TRUE)
names_cols <- colnames(eqtls)[0:(ncol(eqtls)-1)]
colnames(eqtls) <- c(names_cols,'gene_eqtl_Artery_Tibial')
eqtls <- eqtls %>% 
  group_by_at(setdiff(names(eqtls), "gene_eqtl_Artery_Tibial")) %>%
  mutate(gene_eqtl_Artery_Tibial = paste0(gene_eqtl_Artery_Tibial, collapse = ",")) %>%
  ungroup() %>%
  distinct()

eqtls_vsmc <- merge(VSMC_new_diff_genes_fimo, eQTLs$Artery_Coronary[,c('variant_id','gene_id')], all.x=TRUE)
names_cols <- colnames(eqtls_vsmc)[0:(ncol(eqtls_vsmc)-1)]
colnames(eqtls_vsmc) <- c(names_cols,'gene_eqtl_Artery_Coronary')
eqtls_vsmc <- eqtls_vsmc %>% 
  group_by_at(setdiff(names(eqtls_vsmc), "gene_eqtl_Artery_Coronary")) %>%
  mutate(gene_eqtl_Artery_Coronary = paste0(gene_eqtl_Artery_Coronary, collapse = ",")) %>%
  ungroup() %>%
  distinct()


eqtls_vsmc <- merge(eqtls_vsmc, eQTLs$Heart_Atrial_Appendage[,c('variant_id','gene_id')], all.x=TRUE)
names_cols <- colnames(eqtls_vsmc)[0:(ncol(eqtls_vsmc)-1)]
colnames(eqtls_vsmc) <- c(names_cols,'gene_eqtl_Heart_Atrial_Appendage')
eqtls_vsmc <- eqtls_vsmc %>% 
  group_by_at(setdiff(names(eqtls_vsmc), "gene_eqtl_Heart_Atrial_Appendage")) %>%
  mutate(gene_eqtl_Heart_Atrial_Appendage = paste0(gene_eqtl_Heart_Atrial_Appendage, collapse = ",")) %>%
  ungroup() %>%
  distinct()

eqtls_vsmc <- merge(eqtls_vsmc, eQTLs$Heart_Left_Ventricle[,c('variant_id','gene_id')], all.x=TRUE)
names_cols <- colnames(eqtls_vsmc)[0:(ncol(eqtls_vsmc)-1)]
colnames(eqtls_vsmc) <- c(names_cols,'gene_eqtl_Heart_Left_Ventricle')
eqtls_vsmc <- eqtls_vsmc %>% 
  group_by_at(setdiff(names(eqtls_vsmc), "gene_eqtl_Heart_Left_Ventricle")) %>%
  mutate(gene_eqtl_Heart_Left_Ventricle = paste0(gene_eqtl_Heart_Left_Ventricle, collapse = ",")) %>%
  ungroup() %>%
  distinct()

eqtls_vsmc <- merge(eqtls_vsmc, eQTLs$Artery_Aorta[,c('variant_id','gene_id')], all.x=TRUE)
names_cols <- colnames(eqtls_vsmc)[0:(ncol(eqtls_vsmc)-1)]
colnames(eqtls_vsmc) <- c(names_cols,'gene_eqtl_Artery_Aorta')
eqtls_vsmc <- eqtls_vsmc %>% 
  group_by_at(setdiff(names(eqtls_vsmc), "gene_eqtl_Artery_Aorta")) %>%
  mutate(gene_eqtl_Artery_Aorta = paste0(gene_eqtl_Artery_Aorta, collapse = ",")) %>%
  ungroup() %>%
  distinct()

eqtls_vsmc <- merge(eqtls_vsmc, eQTLs$Artery_Tibial[,c('variant_id','gene_id')], all.x=TRUE)
names_cols <- colnames(eqtls_vsmc)[0:(ncol(eqtls_vsmc)-1)]
colnames(eqtls_vsmc) <- c(names_cols,'gene_eqtl_Artery_Tibial')
eqtls_vsmc <- eqtls_vsmc %>% 
  group_by_at(setdiff(names(eqtls_vsmc), "gene_eqtl_Artery_Tibial")) %>%
  mutate(gene_eqtl_Artery_Tibial = paste0(gene_eqtl_Artery_Tibial, collapse = ",")) %>%
  ungroup() %>%
  distinct()

### new overlap ####
eqtls$overlap <- 'Not Overlap'
eqtls$overlap[eqtls$dupe_info %in% overlap$dupe_info] <- 'Differential Overlap'
eqtls$overlap[eqtls$dupe_info %in% VSMC_new_active_diff$dupe_info] <- 'Differential + Active Overlap'

eqtls_vsmc$overlap <- 'Not Overlap'
eqtls_vsmc$overlap[eqtls_vsmc$dupe_info %in% overlap$dupe_info] <- 'Differential Overlap'
eqtls_vsmc$overlap[eqtls_vsmc$dupe_info %in% CM_new_active_diff$dupe_info] <- 'Differential + Active Overlap'

write.table(eqtls, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_snps_final.ranking.withinfo.eqtl.txt', quote=F, row.names=F, sep='\t')
write.table(eqtls_vsmc, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_snps_final.ranking.withinfo.eqtl.txt', quote=F, row.names=F, sep='\t')

### RepeatMasker ####
repeat_masker <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/snps_repeat_masker.txt', sep='\t', header=T)
oligo_info <- vals_significance[,c('dupe_info','pos')]
colnames(oligo_info) <- c('dupe_info','query')

eqtls <- merge(eqtls, oligo_info, all.x=TRUE)
eqtls_vsmc <- merge(eqtls_vsmc, oligo_info, all.x=TRUE)

eqtls$repeatMasker <- 'NoRepeat'
eqtls_vsmc$repeatMasker <- 'NoRepeat'

eqtls$repeatMasker[eqtls$query %in% repeat_masker$query] <- 'Repeat'
eqtls_vsmc$repeatMasker[eqtls_vsmc$query %in% repeat_masker$query] <- 'Repeat'

eqtls <- eqtls[eqtls$query != 'RANDOM',]
eqtls <- eqtls %>% distinct()

eqtls_vsmc <- eqtls_vsmc[eqtls_vsmc$query != 'RANDOM',]
eqtls_vsmc <- eqtls_vsmc %>% distinct()

eqtls <- eqtls[order((eqtls$fdr_comp), decreasing = F),]
eqtls_vsmc <- eqtls_vsmc[order((eqtls_vsmc$fdr_comp), decreasing = F),]

write.table(eqtls, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', quote=F, row.names=F, sep='\t')
write.table(eqtls_vsmc, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', quote=F, row.names=F, sep='\t')

### Do plots ####
eqtls <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header = T)
head(eqtls)

eqtls_vsmc <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header = T)
head(eqtls_vsmc)

### number of functional variant per sentinel ####
sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)

sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]
head(eqtls)

active_sentinel <- merge(eqtls, sentinel_snps, by.x='snp_info', by.y='snp', all.y=TRUE)
active_sentinel$overlap[is.na(active_sentinel$overlap)] <- 'NoFunctionalVar'
active_sentinel$diff[active_sentinel$overlap == 'NoFunctionalVar'] <- 'NoFunctionalVar'
active_sentinel$diff[active_sentinel$overlap != 'NoFunctionalVar'] <- 'FunctionalVar'
active_sentinel$diff <- NULL
active_sentinel$diff[active_sentinel$overlap != 'NoFunctionalVar' & active_sentinel$active == 'Active'] <- 'FunctionalVar'

active_sentinel_table <- as.data.frame(table(active_sentinel$sentinel, active_sentinel$diff))

ggplot(active_sentinel_table, aes(x=Freq)) + 
  geom_density() + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 33)+
  geom_density(alpha=.2, fill="#FF6666") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') +
  xlab('N.FunctionalVariants/sentinel')

### number eqtls ####
eqtls_vsmc$eQTL_Artery_Aorta[is.na(eqtls_vsmc$eQTL_Artery_Aorta)] <- 'NO'
eqtls_vsmc$eQTL_Heart_Atrial_Appendage[is.na(eqtls_vsmc$eQTL_Heart_Atrial_Appendage)] <- 'NO'
eqtls_vsmc$eQTL_Heart_Left_Ventricle[is.na(eqtls_vsmc$eQTL_Heart_Left_Ventricle)] <- 'NO'
eqtls_vsmc$eQTL_Artery_Coronary[is.na(eqtls_vsmc$eQTL_Artery_Coronary)] <- 'NO'
eqtls_vsmc$eQTL_Artery_Tibial[is.na(eqtls_vsmc$eQTL_Artery_Tibial)] <- 'NO'

eqtls$eQTL_Artery_Aorta[is.na(eqtls$eQTL_Artery_Aorta)] <- 'NO'
eqtls$eQTL_Heart_Atrial_Appendage[is.na(eqtls$eQTL_Heart_Atrial_Appendage)] <- 'NO'
eqtls$eQTL_Heart_Left_Ventricle[is.na(eqtls$eQTL_Heart_Left_Ventricle)] <- 'NO'
eqtls$eQTL_Artery_Coronary[is.na(eqtls$eQTL_Artery_Coronary)] <- 'NO'
eqtls$eQTL_Artery_Tibial[is.na(eqtls$eQTL_Artery_Tibial)] <- 'NO'

ggplot(eqtls, aes(x=eQTL_Artery_Tibial, fill= eQTL_Artery_Tibial)) +
  geom_bar(stat = 'count') + 
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') +
  ylab('Number of Functional SNPs')+
  scale_fill_brewer(palette="Set1")

ggplot(eqtls[eqtls$active == 'Active',], aes(x=eQTL_Artery_Tibial, fill= eQTL_Artery_Tibial)) +
  geom_bar(stat = 'count') + 
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') +
  ylab('Number of Functional SNPs')+
  scale_fill_brewer(palette="Set1")

### Distribution of number of different TF ####
N_TF_WT <- strsplit(eqtls$TF_WT, split = ",")
names(N_TF_WT) <- eqtls$snp_info
N_TF_WT <- lapply(N_TF_WT, length )
N_TF_WT <- unlist(N_TF_WT)
eqtls$N_TF_WT <- N_TF_WT 

N_TF_SNP <- strsplit(eqtls$TF_SNP, split = ",")
names(N_TF_SNP) <- eqtls$snp_info
N_TF_SNP <- lapply(N_TF_SNP, length )
N_TF_SNP <- unlist(N_TF_SNP)
eqtls$N_TF_SNP <- N_TF_SNP 

eqtls$N_diff_TF <- eqtls$N_TF_WT + eqtls$N_TF_SNP

ggplot(eqtls, aes(x=active,y=N_diff_motifs, fill=active)) + 
  geom_violin() + 
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
  legend.position = 'none') +
  xlab('Activity') + 
  ylim(c(-30,30))+
  scale_fill_brewer(palette="Set2")

ggplot(eqtls, aes(x=N_diff_TF)) + 
  geom_density() + 
  geom_histogram(aes(y=..density..), colour="black", fill="white", bins = 35)+
  geom_density(alpha=.2, fill="#FF6666") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') +
  xlab('N. diff TF binding')

### logFC of overlapped ####
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
sample_size = eqtls %>% group_by(overlap) %>% dplyr::summarize(num=n())

# Plot
eqtls %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(overlap, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=abs(logFC_comp), fill=overlap)) +
  geom_violin(width=1) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("CM FoldChange") +
  xlab("")

### FIMO motifs between groups ####
fimo_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/03_fimo/diff_TFBS.csv', sep='\t', header = T)
sentinel_fimo <- merge(sentinel_snps, fimo_results, by.x='snp', by.y='SNP', all.x=TRUE)
sentinel_fimo <- sentinel_fimo[sentinel_fimo$Trait %in% c('DBP','SBP','PP'),]

CM_fimo <- merge(sentinel_fimo, eqtls, by.x='snp',by.y='snp_info',all.x=TRUE)
CM_fimo$type <-  'No Functional Variant'
CM_fimo$type[CM_fimo$fdr_comp < 0.05] <-  'Functional Variant'
CM_fimo$type[CM_fimo$active == 'Active'] <-  'Functional Active Variant'

N_TF_WT <- strsplit(CM_fimo$TF_WT.x, split = ",")
names(N_TF_WT) <- CM_fimo$snp_info
N_TF_WT <- lapply(N_TF_WT, length )
N_TF_WT <- unlist(N_TF_WT)
CM_fimo$N_TF_WT <- N_TF_WT 

N_TF_SNP <- strsplit(CM_fimo$TF_SNP.x, split = ",")
names(N_TF_SNP) <- CM_fimo$snp_info
N_TF_SNP <- lapply(N_TF_SNP, length )
N_TF_SNP <- unlist(N_TF_SNP)
CM_fimo$N_TF_SNP <- N_TF_SNP 

CM_fimo$N_diff_TF <- CM_fimo$N_TF_WT + CM_fimo$N_TF_SNP

VSMC_fimo <- merge(sentinel_fimo, eqtls_vsmc, by.x='snp',by.y='snp_info',all.x=TRUE)
VSMC_fimo$type <-  'No Functional Variant'
VSMC_fimo$type[VSMC_fimo$fdr_comp < 0.05] <-  'Functional Variant'
VSMC_fimo$type[VSMC_fimo$active == 'Active'] <-  'Functional Active Variant'

N_TF_WT <- strsplit(VSMC_fimo$TF_WT.x, split = ",")
names(N_TF_WT) <- VSMC_fimo$snp_info
N_TF_WT <- lapply(N_TF_WT, length )
N_TF_WT <- unlist(N_TF_WT)
VSMC_fimo$N_TF_WT <- N_TF_WT 

N_TF_SNP <- strsplit(VSMC_fimo$TF_SNP.x, split = ",")
names(N_TF_SNP) <- VSMC_fimo$snp_info
N_TF_SNP <- lapply(N_TF_SNP, length )
N_TF_SNP <- unlist(N_TF_SNP)
VSMC_fimo$N_TF_SNP <- N_TF_SNP 

VSMC_fimo$N_diff_TF <- VSMC_fimo$N_TF_WT + VSMC_fimo$N_TF_SNP

library(ggpubr)
my_comparisons <- list( c("Functional Active Variant", "Functional Varian"), c("Functional Active Variant", "No Functional Varian"), c("Functional Varian", "No Functional Varian") )

ggplot(CM_fimo, aes(x=type,y=N_diff_TF, fill=type)) + 
  geom_violin() + 
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') +
  stat_compare_means(label = "p.format", method = "t.test",
                     ref.group = "No Functional Variant") +
  #stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  #stat_compare_means(label.y = 50) +
  xlab('Type') + 
  ylim(c(0,30))+
  scale_fill_brewer(palette="Set2")

#### Do some more plots Repeat and TFBS ####
CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)

repeat_masker <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/snps_repeat_masker.info.txt', sep=' ', header=T)
repeat_masker$query <- gsub('_.*','',repeat_masker$query)
repeat_masker <- repeat_masker %>% distinct()

library(RColorBrewer)
n <- 28
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors = sample(col_vector,n)

ggplot(repeat_masker, aes(x=beginR)) + 
  geom_bar(stat = 'count', colour='black', fill=colors) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') +
  coord_flip()+
  xlab('Family of Repeat')

## enrichment of repeat family #####

get_col <- function(fdr){
  if(abs(fdr) < -log10(0.05)){
    col <- "light grey"
  }else{
    col <- col_numeric(rev(brewer.pal(11,"RdBu")),
                       domain = seq(-3,
                                    3,length.out = 11))(fdr)
  }
}

active_snps <- unique(vals_significance_cm$pos[vals_significance_cm$CM_padj<0.05])
active_snps_vsmc <- unique(vals_significance_vsmc$pos[vals_significance_vsmc$VSMC_padj<0.05])

my_fisher <- function(type){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  repeat_fam1 <- repeat_masker$query[repeat_masker$beginR == type]
  repeat_other <- repeat_masker$query[repeat_masker$beginR != type]
   type_counts_original <- length(repeat_fam1[repeat_fam1 %in% CM_results$query])
   other_type_counts_original <- length(repeat_other[repeat_other %in% CM_results$query])
   # type_counts_original <- length(repeat_fam1[repeat_fam1 %in% VSMC_results$query])
   # other_type_counts_original <- length(repeat_other[repeat_other %in% VSMC_results$query])
  type_counts_random <- length(repeat_fam1)-type_counts_original
  other_type_counts_random <- length(repeat_other)-other_type_counts_original
  
  m <- matrix(c(type_counts_original, type_counts_random, other_type_counts_original, other_type_counts_random), 2,2, byrow = T)
  print(m)
  m[is.na(m)] <- 0
  rownames(m) <- c("Intron", "Other")
  colnames(m) <- c("Original","Random")
  f <- fisher.test(m)
  print(f)
  return(list("f" = f,
              "overlap" = x11))
}
# Two-tailed Fisher test
families <- unique(repeat_masker$beginR)
fisher_results <- lapply(families, function(chr) my_fisher(chr))
names(fisher_results) <- families
# Plot
par(mfrow=c(2,2), mar=c(5,7,4,4), mgp=c(4,1,0))
odds_ratio <- list()
adj.P.Val <- list()
# odds ratio
odds_ratio <- sapply(families, function(chr) fisher_results[[chr]][['f']]$estimate )
# adj.P.Val
adj.P.Val <- p.adjust(sapply(families, function(chr) fisher_results[[chr]][['f']]$p.value), method = "BH")
print(paste0(sum(adj.P.Val<0.05)))
print(paste0("Depleted: ", sum(odds_ratio[which(adj.P.Val<0.05)]<1))) #  depleted
print(paste0("Enriched: ", sum(odds_ratio[which(adj.P.Val<0.05)]>1))) # enriched
# fdr
fdr <- -log10(adj.P.Val)
fdr[fdr>3] <- 3 # for nice color scale, min adj.P.Val 0.01
fdr[odds_ratio<1] <- -fdr[odds_ratio<1] # if depleted blue

names(odds_ratio) <- families
odds_ratio <- as.data.frame(odds_ratio)

odds_ratio$fam <- rownames(odds_ratio)
colnames(odds_ratio) <- c('odds_ratio','Family')
odds_ratio$pval <- adj.P.Val

ggplot(odds_ratio, aes(x=odds_ratio, y=Family)) + geom_bar(stat = 'identity')

threshold_OE <- odds_ratio$pval < 0.05  & odds_ratio$odds_ratio > 1
length(which(threshold_OE))
odds_ratio$thr <- threshold_OE 

ggplot(odds_ratio, aes(x=odds_ratio, y=Family)) + 
  geom_bar(aes(fill=thr),stat = 'identity') +
  ggtitle('location enrichment') +
  ylab('')+
  scale_fill_manual(values = c('grey','Dark red'))+
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') + 
  geom_vline(aes(xintercept=1), linetype="dotted")

### overlap of nearest genes ####
library(UpSetR)
nearest_cm <- CM_results$ensemblID
nearest_vsmc <- VSMC_results$ensemblID

listInput <- list(CM = nearest_cm, VSMC = nearest_vsmc)

upset(fromList(listInput), order.by = "freq")

#### plots of sentinel snps with 1 differential ####
CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)

sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

diff_snps <- unique(CM_results[,'snp_info'])

sentinel_snps$active[sentinel_snps$snp %in% diff_snps] <- 'Active'
sentinel_snps$active[!sentinel_snps$snp %in% diff_snps] <- 'Not Active'

sentinel <- unique(sentinel_snps$sentinel)

sentinel_snps_sent <- sentinel_snps[sentinel_snps$snp %in% sentinel,]

grouped <- sentinel_snps[,c('sentinel','Chr','Trait','active')] %>% distinct()

active_sentinel <- unique(grouped$sentinel[grouped$active == 'Active'])

sentinel_snps_sent$active[sentinel_snps_sent$sentinel %in% active_sentinel] <- 'Active'
dist_snps <- as.data.frame(table(sentinel_snps_sent$active))

# Compute the position of labels
dist_snps <- dist_snps %>% 
  #arrange(desc(factor(Var1))) %>%
  mutate(prop = Freq / sum(dist_snps$Freq) *100) %>%
  mutate(ypos = (prop)- 0.1*prop )

ggplot(dist_snps, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, color="white", stat = 'identity') +
  coord_polar("y")+
  geom_text(aes(y = ypos+1, label = paste0(round(Freq),'')), color = "white", size=5) +
  ylab("")+
  xlab("")+
  ggtitle("Sentinel SNPs with at least one LD differential SNP")+
  theme(legend.position="bottom") +
  theme_void() + 
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_brewer(palette="Set2")

