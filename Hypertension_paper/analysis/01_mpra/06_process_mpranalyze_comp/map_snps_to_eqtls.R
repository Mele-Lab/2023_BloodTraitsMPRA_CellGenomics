### GTEx eQTLs mapping to MPRA SNPs #####

#### MPRA data ####
CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)

CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/all_snps_info_risk_CM.txt', sep='\t', header=T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/all_snps_info_risk_VSMC.txt', sep='\t', header=T)

### variant id -- liftover ####
## read liftover coord ####
liftover_coord <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/snp_coord_master.hg38.simple.chr.sorted.bed', sep='\t', header=F)
colnames(liftover_coord) <- c('chr_b38', 'start_b38', 'end_b38','info_snp')

# VSMC_results$info_snp <- paste0(VSMC_results$chr,':',VSMC_results$start,':',VSMC_results$ref,':',VSMC_results$alt)
# CM_results$info_snp <- paste0(CM_results$chr,':',CM_results$start,':',CM_results$ref,':',CM_results$alt)

VSMC_results$info_snp <- VSMC_results$change
CM_results$info_snp <- CM_results$change

VSMC_results <- merge(VSMC_results, liftover_coord, all.x = TRUE)
CM_results <- merge(CM_results, liftover_coord, all.x = TRUE)

### merge with eqtl info ####
VSMC_results$variant_id <- paste0(VSMC_results$chr_b38,'_',as.numeric(VSMC_results$start_b38)+1,'_',gsub(':.*','',VSMC_results$E),'_',gsub('.*:','',VSMC_results$E),'_b38')
CM_results$variant_id <- paste0(CM_results$chr_b38,'_',as.numeric(CM_results$start_b38)+1,'_',gsub(':.*','',CM_results$E),'_',gsub('.*:','',CM_results$E),'_b38')

### GTEx eQTLs ####
gtex_path <- '~/marenostrum/Projects/GTEx_v8/Manuscript/'
TissueInfoFile <- "TissuesInfo.rds" 

TissueInfo <- readRDS(paste0(gtex_path,TissueInfoFile))
tissues <- as.character(TissueInfo$Tissue_id)
subset <- tissues[grep('Arter|Heart',tissues)]

eQTLs_path <- '~/marenostrum/Data/GTEx_v8/cisQTLs/eGenes/GTEx_Analysis_v8_eQTL/'

eQTLs <- lapply(subset, function(tissue) read.table(paste0(eQTLs_path,tissue,".v8.signif_variant_gene_pairs.txt.gz"), sep='\t', header=T) )
names(eQTLs) <- subset

eQTLs_mpra_CM <- lapply(tissues, function(tissue) as.data.frame(CM_results$variant_id[CM_results$variant_id %in% eQTLs[[tissue]]$variant_id]))
names(eQTLs_mpra_CM) <- tissues

eQTLs_mpra_VSMC <- lapply(tissues, function(tissue) as.data.frame(VSMC_results$variant_id[VSMC_results$variant_id %in% eQTLs[[tissue]]$variant_id]))
names(eQTLs_mpra_VSMC) <- tissues

eQTLs_mpra_CM_df <- do.call(rbind,eQTLs_mpra_CM)
eQTLs_mpra_VSMC_df <- do.call(rbind,eQTLs_mpra_VSMC)

colnames(eQTLs_mpra_CM_df) <- c('variant_id')
colnames(eQTLs_mpra_VSMC_df) <- c('variant_id')

eQTLs_mpra_CM_df$tissue <- gsub('\\..*','',rownames(eQTLs_mpra_CM_df))
eQTLs_mpra_VSMC_df$tissue <- gsub('\\..*','',rownames(eQTLs_mpra_VSMC_df))

library(RColorBrewer)
n <- 46
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors = sample(col_vector,n)

n_eqtls <- lapply(tissues, function(tissue) as.data.frame(nrow(eQTLs[[tissue]])))
names(n_eqtls) <- tissues
n_eqtls_df <- do.call(rbind,n_eqtls)
colnames(n_eqtls_df) <- c('n_eqtls')
n_eqtls_df$tissue <- rownames(n_eqtls_df)

df2 <- eQTLs_mpra_CM_df %>% group_by(tissue) %>%   mutate(count_name_occurr = n())
df2 <- merge(df2, n_eqtls_df)

ggplot(df2, aes(x=reorder(tissue,-n_eqtls), fill= tissue)) +
  geom_bar(stat = 'count', fill = colors, colour='black') + 
  xlab('')+
  coord_flip()+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') +
  ylab('Number of Functional SNPs')


#### enrichment of eQTLs ####

#eQTLs_numbers <- as.data.frame(table(eQTLs_mpra_CM_df$tissue))

eQTLs_mpra_CM_funct <- lapply(tissues, function(tissue) as.data.frame(CM_results$variant_id[CM_results$variant_id %in% eQTLs[[tissue]]$variant_id & CM_results$fdr_comp < 0.05]))
names(eQTLs_mpra_CM_funct) <- tissues

eQTLs_mpra_VSMC_funct <- lapply(tissues, function(tissue) as.data.frame(VSMC_results$variant_id[VSMC_results$variant_id %in% eQTLs[[tissue]]$variant_id & VSMC_results$fdr_comp < 0.05]))
names(eQTLs_mpra_VSMC_funct) <- tissues

eQTLs_mpra_CM_funct_df <- do.call(rbind,eQTLs_mpra_CM_funct)
eQTLs_mpra_VSMC_funct_df <- do.call(rbind,eQTLs_mpra_VSMC_funct)

colnames(eQTLs_mpra_CM_funct_df) <- c('variant_id')
colnames(eQTLs_mpra_VSMC_funct_df) <- c('variant_id')

eQTLs_mpra_CM_funct_df$tissue <- gsub('\\..*','',rownames(eQTLs_mpra_CM_funct_df))
eQTLs_mpra_VSMC_funct_df$tissue <- gsub('\\..*','',rownames(eQTLs_mpra_VSMC_funct_df))

eQTLs_mpra_CM_funct_df <- eQTLs_mpra_CM_funct_df[!is.na(eQTLs_mpra_CM_funct_df$variant_id),]
eQTLs_mpra_VSMC_funct_df <- eQTLs_mpra_VSMC_funct_df[!is.na(eQTLs_mpra_VSMC_funct_df$variant_id),]
eQTLs_mpra_CM_df <- eQTLs_mpra_CM_df[!is.na(eQTLs_mpra_CM_df$variant_id),]

eQTLs_numbers <- as.data.frame(table(eQTLs_mpra_CM_df$tissue))
eQTLs_mpra_CM_funct_numbers <- as.data.frame(table(eQTLs_mpra_CM_funct_df$tissue))
eQTLs_mpra_VSMC_funct_numbers <- as.data.frame(table(eQTLs_mpra_VSMC_funct_df$tissue))

get_col <- function(fdr){
  if(abs(fdr) < -log10(0.05)){
    col <- "light grey"
  }else{
    col <- col_numeric(rev(brewer.pal(11,"RdBu")),
                       domain = seq(-3,
                                    3,length.out = 11))(fdr)
  }
}

my_fisher <- function(type){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  type_counts_original <- eQTLs_numbers$Freq[eQTLs_numbers$Var1 == type] - eQTLs_mpra_CM_funct_numbers$Freq[eQTLs_mpra_CM_funct_numbers$Var1 == type]
  other_type_counts_original <- sum(eQTLs_numbers$Freq[eQTLs_numbers$Var1 != type]) - sum(eQTLs_mpra_CM_funct_numbers$Freq[eQTLs_mpra_CM_funct_numbers$Var1 != type])
  type_counts_random <- eQTLs_mpra_CM_funct_numbers$Freq[eQTLs_mpra_CM_funct_numbers$Var1 == type]
  other_type_counts_random <- sum(eQTLs_mpra_CM_funct_numbers$Freq[eQTLs_mpra_CM_funct_numbers$Var1 != type])
  
  m <- matrix(c(type_counts_random, type_counts_original, other_type_counts_random, other_type_counts_original), 2,2, byrow = T)
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
chrom <- unique(eQTLs_numbers$Var1)
fisher_results <- lapply(chrom, function(chr) my_fisher(chr))
names(fisher_results) <- chrom
# Plot
par(mfrow=c(2,2), mar=c(5,7,4,4), mgp=c(4,1,0))
odds_ratio <- list()
adj.P.Val <- list()
# odds ratio
odds_ratio <- sapply(chrom, function(chr) fisher_results[[chr]][['f']]$estimate )
# adj.P.Val
adj.P.Val <- p.adjust(sapply(chrom, function(chr) fisher_results[[chr]][['f']]$p.value), method = "BH")
print(paste0(sum(adj.P.Val<0.05)))
print(paste0("Depleted: ", sum(odds_ratio[which(adj.P.Val<0.05)]<1))) #  depleted
print(paste0("Enriched: ", sum(odds_ratio[which(adj.P.Val<0.05)]>1))) # enriched
# fdr
fdr <- -log10(adj.P.Val)
fdr[fdr>3] <- 3 # for nice color scale, min adj.P.Val 0.01
fdr[odds_ratio<1] <- -fdr[odds_ratio<1] # if depleted blue

names(odds_ratio) <- chrom
odds_ratio <- as.data.frame(odds_ratio)

odds_ratio$chrom <- rownames(odds_ratio)
colnames(odds_ratio) <- c('odds_ratio','chrom')
odds_ratio$pval <- adj.P.Val

ggplot(odds_ratio, aes(x=odds_ratio, y=chrom)) + geom_bar(stat = 'identity')



threshold_OE <- odds_ratio$pval < 0.05  & odds_ratio$odds_ratio > 1
length(which(threshold_OE))
odds_ratio$thr <- threshold_OE 

ggplot(odds_ratio, aes(x=odds_ratio, y=chrom)) + 
  geom_bar(aes(fill=thr),stat = 'identity') +
  ggtitle('location enrichment') +
  ylab('')+
  scale_fill_manual(values = c('grey','Dark red'))+
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') + 
  geom_vline(aes(xintercept=1), linetype="dotted")


#### do plots of eqtls ####
#df2 <- eQTLs_mpra_CM_df %>% group_by(tissue) %>%   mutate(count_name_occurr = n())
#df2 <- merge(df2, n_eqtls_df)

eQTLs_mpra_CM_funct <- lapply(tissues, function(tissue) as.data.frame(CM_results$variant_id[CM_results$variant_id %in% eQTLs[[tissue]]$variant_id & CM_results$fdr_comp < 0.05]))
names(eQTLs_mpra_CM_funct) <- tissues

eQTLs_mpra_VSMC_funct <- lapply(tissues, function(tissue) as.data.frame(VSMC_results$variant_id[VSMC_results$variant_id %in% eQTLs[[tissue]]$variant_id & VSMC_results$fdr_comp < 0.05]))
names(eQTLs_mpra_VSMC_funct) <- tissues

eQTLs_mpra_CM_funct_df <- do.call(rbind,eQTLs_mpra_CM_funct)
eQTLs_mpra_VSMC_funct_df <- do.call(rbind,eQTLs_mpra_VSMC_funct)

colnames(eQTLs_mpra_CM_funct_df) <- c('variant_id')
colnames(eQTLs_mpra_VSMC_funct_df) <- c('variant_id')

eQTLs_mpra_CM_funct_df$tissue <- gsub('\\..*','',rownames(eQTLs_mpra_CM_funct_df))
eQTLs_mpra_VSMC_funct_df$tissue <- gsub('\\..*','',rownames(eQTLs_mpra_VSMC_funct_df))

eQTLs_mpra_CM_funct_df <- eQTLs_mpra_CM_funct_df[!is.na(eQTLs_mpra_CM_funct_df$variant_id),]
eQTLs_mpra_VSMC_funct_df <- eQTLs_mpra_VSMC_funct_df[!is.na(eQTLs_mpra_VSMC_funct_df$variant_id),]
eQTLs_mpra_CM_df <- eQTLs_mpra_CM_df[!is.na(eQTLs_mpra_CM_df$variant_id),]

eQTLs_numbers <- as.data.frame(table(eQTLs_mpra_CM_df$tissue))
eQTLs_mpra_CM_funct_numbers <- as.data.frame(table(eQTLs_mpra_CM_funct_df$tissue))
eQTLs_mpra_VSMC_funct_numbers <- as.data.frame(table(eQTLs_mpra_VSMC_funct_df$tissue))

DSC_spec_p_full <- plot_grid(number_events, Perc_events, align = "h", nrow = 1, rel_widths = c(0.70, 0.3))

##### new subset GWAS snps  -----
my_data$variant_id <- my_data$info_snp
eqtls <- merge(my_data, eQTLs$Artery_Coronary[,c('variant_id','gene_id')], all.x=TRUE)
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

write.table(eqtls, 'Downloads/2022_01_04_sentinelSNP_eQTL_GTEx_info.txt',
            quote = F, row.names = F, col.names = T, sep = '\t')
