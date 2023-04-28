##### take a look at coloc and susie results -----

### first merge coloc results ----
files <- list.files('../../analysis/06_colocalization/colocQuial_res/', 
                    pattern = 'spread')

all_coloc_res_s <- lapply(c(1:3), function(x) read.table(paste0('../../analysis/06_colocalization/colocQuial_res/',
                                                              files[x]), header = T, sep = '\t'))

all_coloc_res_s_df <- do.call('rbind.data.frame', all_coloc_res_s)
nrow(all_coloc_res_s_df[all_coloc_res_s_df$cond.PP.H4.abf>0.8,])
length(unique(all_coloc_res_s_df$SNP[all_coloc_res_s_df$cond.PP.H4.abf>0.8]))
length(unique(all_coloc_res_s_df$Gene[all_coloc_res_s_df$cond.PP.H4.abf>0.8]))

###trait-related tissues ####3
## info from GTEx data portal
gtex_path <- '~/marenostrum/Projects/GTEx_v8/Manuscript/'
TissueInfoFile <- "TissuesInfo.rds" 

TissueInfo <- readRDS(paste0(gtex_path,TissueInfoFile))
tissues <- as.character(TissueInfo$Tissue_id)
subset <- tissues[grep('Arter|Heart',tissues)]
subset <- tissues[grep('Brain|Pitu|Gland|Thyr|Panc|Nerve|Kidn|Uret|Blad',tissues)]

all_coloc_res_s_df$tissue <- gsub('ENSG[0-9]*\\.[0-9]*_','', all_coloc_res_s_df$GeneID.Tissue)
tissue_blood <- all_coloc_res_s_df[all_coloc_res_s_df$tissue %in% subset,]
nrow(tissue_blood[tissue_blood$cond.PP.H4.abf>0.8,])
length(unique(tissue_blood$SNP[tissue_blood$cond.PP.H4.abf>0.8]))
length(unique(tissue_blood$Gene[tissue_blood$cond.PP.H4.abf>0.8]))

write.table(unique(tissue_blood$Gene[tissue_blood$cond.PP.H4.abf>0.8]), 
            '../../analysis/06_colocalization/genes_coloc_blood_tissues.txt',
            quote = F, row.names = F, col.names = F, sep = '\n')

### first merge susie results ----
files <- list.files('../../analysis/06_colocalization/colocQuial_res/', 
                    pattern = 'susie')

all_coloc_res <- lapply(c(1:3), function(x) read.table(paste0('../../analysis/06_colocalization/colocQuial_res/',
                                                              files[x]), header = F, sep = '\t', skip = 1))

all_coloc_res_df <- do.call('rbind.data.frame', all_coloc_res)
cols <- c('SNP','Gene','GeneID-Tissue','Trait','number','nsnps','hit1','hit2','PP.H0.abf','PP.H1.abf',
          'PP.H2.abf','PP.H3.abf','PP.H4.abf','idx1','idx2')
colnames(all_coloc_res_df) <- cols
all_coloc_res_df$cond.PP.H4.abf <- all_coloc_res_df$PP.H4.abf/(all_coloc_res_df$PP.H3.abf+all_coloc_res_df$PP.H4.abf)

nrow(all_coloc_res_df[all_coloc_res_df$cond.PP.H4.abf>0.8,])
length(unique(all_coloc_res_df$SNP[all_coloc_res_df$cond.PP.H4.abf>0.8]))
length(unique(all_coloc_res_df$Gene[all_coloc_res_df$cond.PP.H4.abf>0.8]))

all_coloc_res_df$tissue <- gsub('ENSG[0-9]*\\.[0-9]*_','', all_coloc_res_df$`GeneID-Tissue`)
tissue_blood_s <- all_coloc_res_df[all_coloc_res_df$tissue %in% subset,]
nrow(tissue_blood_s[tissue_blood_s$cond.PP.H4.abf>0.8,])
length(unique(tissue_blood_s$SNP[tissue_blood_s$cond.PP.H4.abf>0.8]))
length(unique(tissue_blood_s$Gene[tissue_blood_s$cond.PP.H4.abf>0.8]))


#### plots -----
## read differential results 
CM_results_s <- read.table('../../analysis/CM_eqtls_rep_cons_dhs_tbs_scores.txt', sep = '\t',
                           header = T)

VSMC_results_s <- read.table('../../analysis/VSMC_eqtls_rep_cons_dhs_tbs_scores.txt', sep = '\t',
                             header = T)

### plot general overview results  #####
nrow(tissue_blood[tissue_blood$cond.PP.H4.abf>0.8,])
length(unique(tissue_blood$SNP[tissue_blood$cond.PP.H4.abf>0.8]))
length(unique(tissue_blood$Gene[tissue_blood$cond.PP.H4.abf>0.8]))

data <- data.frame(
  name = c("coloc","reg-coloc"),
  eGenes = c(length(unique(tissue_blood$Gene[tissue_blood$cond.PP.H4.abf>0.8])),
             length(unique(tissue_blood$Gene[tissue_blood$cond.PP.H4.abf>0.8 & tissue_blood$SNP %in% CM_results_s$sentinel])))
)

# Increase bottom margin
par(mar=c(6,4,4,4))


# Basic Barplot
my_bar <- barplot(data$eGenes , border=F , names.arg=data$name , 
                  las=2 , 
                  col=c(rgb(0.3,0.1,0.4,0.6) , rgb(0.3,0.5,0.4,0.6)) , 
                  ylim=c(0,200) , 
                  main="" , ylab = 'number of eGenes')

# Add the text 
text(my_bar, data$eGenes-10 , paste( data$eGenes, sep="") ,cex=1) 

#Legende
legend("topleft", legend = c("coloc","regulatory-coloc") , 
       col = c(rgb(0.3,0.1,0.4,0.6) , rgb(0.3,0.5,0.4,0.6)) , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.05, 0.05))

data <- data.frame(
  name = c("coloc","reg-coloc"),
  eGenes = c(length(unique(tissue_blood$SNP[tissue_blood$cond.PP.H4.abf>0.8])),
             length(unique(tissue_blood$SNP[tissue_blood$cond.PP.H4.abf>0.8 & tissue_blood$SNP %in% CM_results_s$sentinel])))
)

# Increase bottom margin
par(mar=c(4,6,4,4))


# Basic Barplot
my_bar <- barplot(data$eGenes , border=F , names.arg=data$name , 
                  las=1 , 
                  col=c(rgb(0.3,0.1,0.4,0.6) , rgb(0.3,0.5,0.4,0.6)),
                  main="" , horiz = T, xlim = c(0,80), xlab = 'number of sentinel variants')

# Add the text 
text(my_bar, c(1,2) , paste( data$eGenes, sep="") ,cex=1) 

#Legende
legend("topleft", legend = c("coloc","regulatory-coloc") , 
       col = c(rgb(0.3,0.1,0.4,0.6) , rgb(0.3,0.5,0.4,0.6)) , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.05, 0.05))


### rs11191548 SFXN2
#### rs3774372 ULK
#### rs59476975 MAP4
#provide the header to the eqtl data
eQTL_all_header = c("chrom_b38", "chromStart_b38", "chromEnd_b38", "eGeneID", "A1_eqtl", "A2_eqtl", "build", "tss_distance", "ma_samples", "ma_count", "maf", "pvalue_eQTL", "slope", "slope_se")
#name of column in header representing geneID
eQTL_all_geneID = "eGeneID" 
#name of column in header representing chromosome
eQTL_all_chrom = "chrom_b38"
#name of column in header representing end coordinate 
eQTL_all_chromEnd = "chromEnd_b38"
#name of column in header representing p-value
eQTL_all_pvalue = "pvalue_eQTL"

leadSNP_DF_Res = read.delim('../../analysis/06_colocalization/colocQuial_res/rs3774372/ULK4_ENSG00000168038.10_Artery_Aorta_PP_coloc_results_full.txt')
leadSNP_DF = read.delim('../../analysis/06_colocalization/colocQuial_res/rs3774372/ULK4_ENSG00000168038.10_Artery_Aorta_PP_coloc_input_data.txt')
leadSNP_DF[["chrom_b38"]] = as.integer(gsub('[a-zA-Z]', '', leadSNP_DF[["chrom_b38"]])) 
#leadSNP_DF = leadSNP_DF %>% dplyr::select(SNP, all_of(eQTL_all_chrom), all_of(trait_BPcol), all_of(eQTL_all_pvalue), all_of(trait_Pcol))

# library(locuscomparer)
# 
# gwas <- as.data.frame(leadSNP_DF[,c('SNP','P')])
# colnames(gwas) <- c('rsid','pval')
# eqtl <- as.data.frame(leadSNP_DF[,c('SNP','pvalue_eQTL')])
# colnames(eqtl) <- c('rsid','pval')
# 
# out_prefix <- '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/06_colocalization/colocQuial_res/rs11191548/NT5C2_ENSG00000076685.18_Artery_Tibial_SBP'
# 
# locuscomp_gwas_str = paste(out_prefix,"gwas_locuscomp.txt",sep="_")
# locuscomp_eqtl_str = paste(out_prefix,"eqtl_locuscomp.txt",sep="_")
# 
# write.table(gwas, file=locuscomp_gwas_str, sep="\t", row.names=F, quote=F, col.names = T)
# write.table(eqtl, file=locuscomp_eqtl_str, sep="\t", row.names=F, quote=F, col.names = T)
# 
# RA_plot <- locuscompare(in_fn1 = locuscomp_gwas_str, 
#                         in_fn2 = locuscomp_eqtl_str, 
#                         title = 'GWAS', title2 = 'eQTL')

head(leadSNP_DF)
leadSNP_DF$mpra_reg <- 'not reg variant'
leadSNP_DF$mpra_reg[leadSNP_DF$SNP %in% CM_results_s$snp_info] <- 'reg in CM'
leadSNP_DF$mpra_reg[leadSNP_DF$SNP %in% VSMC_results_s$snp_info] <- 'reg in VSMC'
leadSNP_DF$mpra_reg[leadSNP_DF$SNP %in% CM_results_s$snp_info & leadSNP_DF$SNP %in% VSMC_results_s$snp_info] <- 'reg in both'
table(leadSNP_DF$mpra_reg)

ggplot(leadSNP_DF, aes(x=-log10(P), y = -log10(pvalue_eQTL)))+
  geom_point(alpha=0.3)+
  geom_point(aes(color=mpra_reg))+
  ggrepel::geom_text_repel( 
    data=leadSNP_DF %>% filter(mpra_reg == 'reg in both'), # Filter data first
    aes(label=SNP),
    box.padding = 0.5,
    point.padding = 0.2,
    nudge_x = .15,
    nudge_y = .5,
    segment.curvature = -0.1,
    arrow = arrow(length = unit(0.015, "npc"))
  ) + 
  theme_classic()


#### plot genomic region ####
#### plot coloc results ------
library(Gviz)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

snp_coord <- list(paste0('1:',15464479,':',15464479+1), paste0('10:',102854593,':',102854593+1),
                  paste0('8:',140049552,':',140049552+1), paste0('2:',164186950,':',164186950+1))

names(snp_coord) <- levels(metadata$V3)[2:5]
snp_coord_df = lapply(snp_coord, function (x) {res=strsplit(x, ':')}) %>%
  unlist %>%
  as.numeric %>%
  matrix(ncol=3, byrow=T) %>%
  as.data.frame %>%
  dplyr::select(chrom=V1, start=V2, end=V3) %>%
  mutate(chrom=paste0('chr', chrom))

rownames(snp_coord_df) <- levels(metadata$V3)[2:5]

for (snp in levels(metadata$V3)[2:5]) {
  mycoords.gr <- makeGRangesFromDataFrame(coord_snps[snp,])
  genes <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), mycoords.gr)
  entrez_DE <- mapIds(org.Hs.eg.db,keys=genes$gene_id,column="ENSEMBL", keytype="ENTREZID", multiVals="first")
  symbol_DE <- mapIds(org.Hs.eg.db,keys=genes$gene_id,column="SYMBOL", keytype="ENTREZID", multiVals="first")
  
  genes$ensembl <- entrez_DE
  genes$symbol <- symbol_DE
  
  snp_coord_gr <- makeGRangesFromDataFrame(snp_coord_df[snp,])
  atrack <- AnnotationTrack(snp_coord_gr, name = 'SNP', id = ' ')
  
  chr <- as.character(unique(seqnames(mycoords.gr)))
  gen <- genome(genes)
  genes <- genes[genes$ensembl %in% rownames(list_results[[snp]]),]
  #atrack <- AnnotationTrack(genes, name = "genes", id = genes$ensembl)
  
  grtrack <- GeneRegionTrack(genes, genome = gen,
                             chromosome = chr, name = "Gene Model", 
                             transcriptAnnotation = "symbol",
                             background.panel = "#FFFEDB",
                             background.title = "darkblue")
  
  
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = "hg38", chromosome = chr)
  
  list_results[[snp]]$ensembl <- rownames(list_results[[snp]])
  merged_res <- merge(list_results[[snp]], genes)
  merged_res$sig <- 'Not sig'
  merged_res$sig[merged_res$perm_pval<0.05] <-'sig'
  dtrack <- DataTrack(data = merged_res[,c('rho')], start = merged_res$start,
                      end = merged_res$end, chromosome = merged_res$seqnames, genome = "hg38", 
                      name = "Rho")
  
  plotTracks(list(itrack, gtrack, atrack, grtrack, dtrack),featureAnnotation = "id", 
             fontcolor.feature = "darkblue", just.feature = 'above') 
  
  #plotTracks(dtrack)
  
}
