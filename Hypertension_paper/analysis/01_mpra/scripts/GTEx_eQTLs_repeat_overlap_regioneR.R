###### overlap eQTL GTEx SNPs with repeats -------
### repeat sequences -----
### RepeatMasker ####
repeat_masker <- read.delim('~/marenostrum/Data/RepeatMAsker/hg38_repeatmasker_clean.txt', sep='',skip = 3,header=F)
head(repeat_masker)

###gtex snps ----
### GTEx eQTLs ####
gtex_path <- '~/marenostrum/Projects/GTEx_v8/Manuscript/'
TissueInfoFile <- "TissuesInfo.rds" 

TissueInfo <- readRDS(paste0(gtex_path,TissueInfoFile))
tissues <- as.character(TissueInfo$Tissue_id)
subset <- tissues[grep('Arter|Heart',tissues)]

eQTLs_path <- '~/marenostrum_scratch/GTEx/v8/cis_QTLs/cis_eQTLs/GTEx_Analysis_v8_eQTL/'

eQTLs <- lapply(subset, function(tissue) read.delim(paste0(eQTLs_path,tissue,".v8.signif_variant_gene_pairs.txt.gz"), sep='\t', header=T) )
names(eQTLs) <- subset

eQTLs <- lapply(subset, function(tissue) eQTLs[[tissue]] %>%
                  separate(variant_id, c("chrom", "start","ref","alt",'genome'), "_"))
names(eQTLs) <- subset

## plot n detected eqtls #### 
eqtls_number <- as.data.frame(sapply(subset, function(tissue) nrow(eQTLs[[tissue]])))
colnames(eqtls_number) <- c('eQTLs')
eqtls_number$tissue <- rownames(eqtls_number)
eqtls_number$new_name <- c('Coronary artery','Atrial append.','Left ventricle','Aorta','Tibial artery')
ggplot(eqtls_number, aes(y=reorder(new_name, -eQTLs), x = eQTLs)) + 
  geom_bar(stat = 'identity') +
  theme_bw() + ylab('') + xlab('Nº of eQTLs')

liftover_coord <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/snp_coord_master.hg38.simple.chr.sorted.bed', sep='\t', header=F)
colnames(liftover_coord) <- c('chr_b38', 'start_b38', 'end_b38','info_snp')

sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]
liftover_coord <- merge(liftover_coord, sentinel_snps, by.x='info_snp', by.y='change')
liftover_coord <- liftover_coord %>% distinct()

liftover_coord$variant_id <- paste0(liftover_coord$chr_b38,'_',as.numeric(liftover_coord$start_b38)+1,'_',gsub(':.*','',liftover_coord$E),'_',gsub('.*:','',liftover_coord$E),'_b38')
head(liftover_coord)

eQTLs_mpra <- lapply(subset, function(tissue) as.data.frame(merge(liftover_coord, eQTLs[[tissue]][,c("chrom", "start")], 
                                                                  by.x=c('chr_b38', 'start_b38'), 
                                                                  by.y=c("chrom", "start"))))
names(eQTLs_mpra) <- subset

### we need to get info in regions -----
head(eQTLs)
all.genes = liftover_coord %>%
  dplyr::select(chrom=chr_b38, start=start_b38, end=end_b38) %>%
  makeGRangesFromDataFrame

repeat_gr = repeat_masker[,c('V5','V6','V7')] %>%
  as.data.frame %>%
  dplyr::select(chrom=V5, start=V6, end=V7) %>%
  makeGRangesFromDataFrame

eqtls_gr <- liftover_coord[liftover_coord$snp %in% eQTLs_mpra$Artery_Tibial$snp,] %>%
  dplyr::select(chrom=chr_b38, start=start_b38, end=end_b38) %>%
  makeGRangesFromDataFrame

##### make test -------
library(regioneR)
pt <- permTest(A=eqtls_gr, B=repeat_gr, randomize.function=resampleRegions, universe=all.genes,
               evaluate.function=numOverlaps, ntimes=5000,verbose=FALSE)
summary(pt)
pt
plot(pt)


#### all eqtls assayed ----
eqtls_tested <- read.table('~/marenostrum/Data/GTEx_v8/cisQTLs/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt', 
                           sep = '\t', header = T)
