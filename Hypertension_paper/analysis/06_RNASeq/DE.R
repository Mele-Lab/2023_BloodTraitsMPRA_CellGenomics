library("DESeq2")
library(dplyr)
library(tidyr)

### read data ----
counts_PE <- read.table('../../analysis/RNASeq/PE.counts.all.txt', sep = '\t', header = T)
head(counts_PE)
counts_PE <- counts_PE[,colnames(counts_PE)[!colnames(counts_PE) %in% c( "rs3824754_73_1", "rs3824754_74_2")]]

### read info -----
metadata <- read.table('../../analysis/RNASeq/metadata.PE.txt', sep = '\t', header = F)
head(metadata)
metadata <- metadata[metadata$V4 !=11,]
rownames(metadata) <- metadata$V2

### merge counts from same gene
counts_PE$ensembl <- gsub('\\..*','', counts_PE$Geneid)
counts_PE_sum <- counts_PE %>% 
 group_by(ensembl) %>%
 summarise_if(is.numeric, sum, na.rm = TRUE)

counts_PE_sum <- as.data.frame(counts_PE_sum)
rownames(counts_PE_sum) <- counts_PE_sum$ensembl


#### DE all together ------
suppressMessages(library(edgeR))
suppressMessages(library(limma))
options(connectionObserver = NULL)
library(org.Hs.eg.db)

### filter lowly expressed genes ####
# Obtain CPMs
myCPM <- cpm(counts_PE_sum[,c(2:21)])
# Have a look at the output
head(myCPM)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)
# Summary of how many TRUEs there are in each row
# There are 10975 genes that have TRUEs in all 14 samples.
table(rowSums(thresh))
# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
summary(keep)
# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
plot(myCPM[,1],counts_PE_sum[,2])
# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],counts_PE_sum[,2],ylim=c(0,50),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.8)
abline(h=15)
### filter genes
geneExpr = DGEList( counts_PE_sum[,c(2:21)] )
geneExpr$samples$snp <- metadata$V3
geneExpr$samples$cond <- metadata$V5
geneExpr = calcNormFactors( geneExpr )
y <- geneExpr[keep, keep.lib.sizes=FALSE]


metadata$V4 <- as.factor(metadata$V4)
metadata$V5 <- as.factor(metadata$V5)
metadata$V3 <- as.factor(metadata$V3)


#### test correlations between allelic frequency and expression ------
#### window of +- 500kb
### df of gene / distance / rho

#### first calculate correlation between expression and allelic frequency
metadata <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/RNASeq/metadata.PE.txt', sep = '\t', header = F)
head(metadata)
rownames(metadata) <- metadata$V2
metadata$V5 <- as.numeric(metadata$V5)
metadata$V5[metadata$V3 == 'WT'] <- 0
metadata$V3 <- as.factor(metadata$V3)
summary(metadata)

list_results <- list()
for (snp in levels(metadata$V3)[1:4]) {
  metadata_filt <- metadata[metadata$V3 %in% c('WT',snp),]
  results <- list()
  for (gene in rownames(y$counts)) {
    expresion <- y$counts[gene,metadata_filt$V2]
    varfreq <- metadata_filt$V5
    res2 <-cor.test(expresion,varfreq,  method = "spearman")
    results[[gene]] <- list("rho" = res2$estimate[[1]],
                            "pvalue" = res2$p.value)
  }
  list_results[[snp]] <- do.call('rbind.data.frame', results)
}

list_results <- list()
library(DESeq2)
for (snp in levels(metadata$V3)[1:4]) {
  metadata_filt <- metadata[metadata$V3 %in% c('WT',snp),]
  gene_expr <- counts_PE_sum[,rownames(metadata_filt)]
  
  myCPM <- cpm(gene_expr)
  thresh <- myCPM > 0.5

  keep <- rowSums(thresh) >= 2
  summary(keep)
  
  ### filter genes
  geneExpr = DGEList( gene_expr )
  geneExpr$samples$snp <- metadata_filt$V3
  geneExpr$samples$cond <- metadata_filt$V5
  geneExpr = calcNormFactors( geneExpr )
  y <- geneExpr[keep, keep.lib.sizes=FALSE]
  
  m1 <- model.matrix(~ V5 , metadata_filt)
  dds <- DESeqDataSetFromMatrix(countData = y$counts,
                                colData = metadata_filt,
                                design = m1)
  
  dds <- estimateSizeFactors(dds)
  human_counts <- counts(dds, normalized=TRUE)
  results <- list()
  for (gene in rownames(human_counts)) {
    expresion <- human_counts[gene,metadata_filt$V2]
    varfreq <- metadata_filt$V5
    res2 <-cor.test(expresion,varfreq,  method = "pearson")
    results[[gene]] <- list("rho" = res2$estimate[[1]],
                            "pvalue" = res2$p.value)
  }
  list_results[[snp]] <- do.call('rbind.data.frame', results)
}


for (snp in levels(metadata$V3)[1:4]) {

  adj.P.Val <- p.adjust(list_results[[snp]]$pvalue, method = "BH")
  list_results[[snp]]$padj <- adj.P.Val
  
}

###### permutation test 1 -----
#### Correlate pearson distance ~ rho 
#### 50000 or 100000 permuted rho values vs distance (or independent_test)
## histogram pearson correlations
BiocManager::install('Homo.sapiens')
library(Homo.sapiens)
library(biomaRt)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

filterlist <- list(paste0('1:',15464479-500000,':',15464479+500000), paste0('10:',102854593-500000,':',102854593+500000),
                   paste0('8:',140049552-500000,':',140049552+500000), paste0('2:',164186950-500000,':',164186950+500000))
filterlist <- list(paste0('1:',15464479-5000000,':',15464479+5000000), paste0('10:',102854593-5000000,':',102854593+5000000),
                   paste0('8:',140049552-5000000,':',140049552+5000000), paste0('2:',164186950-5000000,':',164186950+5000000))

names(filterlist) <- levels(metadata$V3)[1:4]

coord_snps = lapply(filterlist, function (x) {res=strsplit(x, ':')}) %>%
  unlist %>%
  as.numeric %>%
  matrix(ncol=3, byrow=T) %>%
  as.data.frame %>%
  dplyr::select(chrom=V1, start=V2, end=V3) %>%
  mutate(chrom=paste0('chr', chrom))
rownames(coord_snps) <- levels(metadata$V3)[1:4]

genes_chrom <- list()
library(biomaRt)
mart.hs <- useMart("ensembl", "hsapiens_gene_ensembl")

for (snp in levels(metadata$V3)[1:4]) {
  genes <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position","hgnc_symbol"),
        filters    = c("chromosome_name"), 
        values     = coord_snps[snp,'chrom'], 
        mart       = mart.hs)
  genes_chrom[[snp]] <- genes$ensembl_gene_id
}

# genes <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), mycoords.gr)
# entrez_DE <- mapIds(org.Hs.eg.db,keys=genes$gene_id,column="ENSEMBL", keytype="ENTREZID", multiVals="first")
# genes$ensembl <- entrez_DE

##### permutation test 2 -----
### Compute rho genome-wide == as "permuted"
### plot density of genome wide rho and observed specific rho of nearby gene 
permutation_results <- list()
for (snp in levels(metadata$V3)[1:4]) {
  mycoords.gr <- makeGRangesFromDataFrame(coord_snps[snp,])
  genes <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), mycoords.gr)
  entrez_DE <- mapIds(org.Hs.eg.db,keys=genes$gene_id,column="ENSEMBL", keytype="ENTREZID", multiVals="first")
  genes$ensembl <- entrez_DE
  
  list_results[[snp]]$perm_pval <- 1
  list_results[[snp]][,'perm_pval_padj'] <- 1
  for (gene in unique(genes$ensembl[!is.na(genes$ensembl)])) {
    if (gene %in% rownames(list_results[[snp]])) {
      rho <- list_results[[snp]][gene,'rho']
      higher <- 0 + nrow(list_results[[snp]][abs(list_results[[snp]]$rho )> abs(rho),])
      permuted_pvlue <- higher/nrow(list_results[[snp]])
      list_results[[snp]][gene,'perm_pval'] <- permuted_pvlue
      
      
      hist((list_results[[snp]]$rho), breaks=50, col='grey', main=paste0(snp,'-',gene,';',permuted_pvlue), las=1, xlab='rho')
      abline(v=(rho), lwd=3, col="red") 
    }
    
  }
  list_results[[snp]][unique(genes$ensembl[!is.na(genes$ensembl)]),'perm_pval_padj'] <- p.adjust(list_results[[snp]][unique(genes$ensembl[!is.na(genes$ensembl)]),'perm_pval'], method = "BH")
  

}

#### pval for all genes!
for (snp in levels(metadata$V3)[2:5]) {
  
  list_results[[snp]]$perm_pval_all <- 1
  for (gene in unique(rownames(list_results[[snp]]))) {
    rho <- list_results[[snp]][gene,'rho']
    higher <- 0 + nrow(list_results[[snp]][abs(list_results[[snp]]$rho )> abs(rho),])
    permuted_pvlue <- higher/nrow(list_results[[snp]])
    list_results[[snp]][gene,'perm_pval_all'] <- permuted_pvlue
    
  }
  
  
}

for (snp in levels(metadata$V3)[2:5]) {
  
  adj.P.Val <- p.adjust(list_results[[snp]]$perm_pval_all, method = "BH")
  list_results[[snp]]$perm_pval_all_adj <- adj.P.Val
  
}

lapply(levels(metadata$V3)[1:4], function(snp) print(sum(list_results[[snp]]$perm_pval_padj<0.05), na.omit=T))
lapply(levels(metadata$V3)[1:4], function(snp) print(sum(list_results[[snp]]$perm_pval_all_adj<0.05)))



### plot distribution of rho values of nearby genes ----
library(ggrepel)
library(ggridges)
for (snp in levels(metadata$V3)[2:5]) {
  filt <- list_results[[snp]]
  filt <- filt[!is.na(filt$rho),]
  list_results[[snp]] <- filt
  # print(ggplot(filt, aes(x = rho, y = gene_interact, fill = gene_interact)) +
  #   geom_density_ridges() +
  #   theme_ridges() + 
  #   ggtitle(snp) + 
  #   theme(legend.position = "none"))
  
  p2 <- ggplot(data=filt, aes(x=abs(rho), group=gene_interact, fill=gene_interact)) +
    geom_density(adjust=1.5, alpha=.4) +
    ggtitle(snp) + 
    theme_minimal()
  print(p2)
  
}

library(ggpubr)
for (snp in levels(metadata$V3)[2:5]) {
  filt <- list_results[[snp]]
  filt <- filt[!is.na(filt$rho),]
  list_results[[snp]] <- filt
  
  filt <- filt[filt$gene_chrom == 'YES',]
  # print(ggplot(filt, aes(x = rho, y = gene_interact, fill = gene_interact)) +
  #   geom_density_ridges() +
  #   theme_ridges() + 
  #   ggtitle(snp) + 
  #   theme(legend.position = "none"))
  
  sample_size = filt %>% group_by(gene_interact) %>% summarize(num=n())
  
  # Plot
  filt <- filt %>%
    left_join(sample_size) %>%
    mutate(myaxis = paste0(gene_interact, "\n", "n=", num))
  
  p2 <- ggplot(data = filt, aes(x=myaxis, y=abs(rho), fill=gene_interact)) +
    geom_violin(width=1.) +
    geom_boxplot(width=0.1, color="black", alpha=0.5, fill='light grey') +
    #scale_fill_viridis(discrete = TRUE) +
    ggtitle(snp) + 
    theme_minimal() +
    xlab("") + stat_compare_means()
  
  # p2 <- ggplot(data=filt, aes(y=abs(rho), x=gene_interact, fill=gene_interact)) +
  #   geom_violin() + #geom_jitter()+
  #   geom_boxplot(width=0.05, color='black', fill='light grey')+
  #   ggtitle(snp) + 
  #   theme_minimal()
  print(p2)
  
  p3 <- ggplot(data = filt, aes(x=abs(rho), colour=gene_interact)) +
    stat_ecdf() +
    theme_minimal() 
  print(p3)
  # 
}

lapply(levels(metadata$V3)[2:5], function(snp) print(nrow(list_results[[snp]][list_results[[snp]]$perm_pval_all<0.05 & list_results[[snp]]$gene_interact == 'YES',])))

genes <- lapply(levels(metadata$V3)[2:5], function(snp) print(rownames(list_results[[snp]][list_results[[snp]]$perm_pval_all<0.05 & list_results[[snp]]$gene_interact == 'YES',])))
names(genes) <- levels(metadata$V3)[2:5]
saveRDS(genes, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/RNASeq/genes_in_contact_significant_rho.rds')

#### plot rho from 1MB window ------
library(Gviz)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

coord_snps = lapply(filterlist, function (x) {res=strsplit(x, ':')}) %>%
  unlist %>%
  as.numeric %>%
  matrix(ncol=3, byrow=T) %>%
  as.data.frame %>%
  dplyr::select(chrom=V1, start=V2, end=V3) %>%
  mutate(chrom=paste0('chr', chrom))
rownames(coord_snps) <- levels(metadata$V3)[2:5]

for (snp in levels(metadata$V3)[2:5]) {
  mycoords.gr <- makeGRangesFromDataFrame(coord_snps[snp,])
  genes <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), mycoords.gr)
  entrez_DE <- mapIds(org.Hs.eg.db,keys=genes$gene_id,column="ENSEMBL", keytype="ENTREZID", multiVals="first")
  genes$ensembl <- entrez_DE
  
  chr <- as.character(unique(seqnames(mycoords.gr)))
  #gen <- genome(genes)
  genes <- genes[genes$ensembl %in% rownames(list_results[[snp]]),]
  atrack <- AnnotationTrack(genes, name = "genes", id = genes$ensembl)
  
  #plotTracks(itrack)
  
  gtrack <- GenomeAxisTrack()
  #plotTracks(list(gtrack, atrack))
  itrack <- IdeogramTrack(genome = "hg38", chromosome = chr)
  #plotTracks(list(itrack, gtrack,atrack))
  
  list_results[[snp]]$ensembl <- rownames(list_results[[snp]])
  merged_res <- merge(list_results[[snp]], genes)
  merged_res$sig <- 'Not sig'
  merged_res$sig[merged_res$perm_pval_padjl<0.05] <-'sig'
  dtrack <- DataTrack(data = merged_res[,c('rho')], start = merged_res$start,
                      end = merged_res$end, chromosome = merged_res$seqnames, genome = "hg38", 
                      name = "Rho")
  data(geneModels)

  plotTracks(list(itrack, gtrack, atrack, dtrack),featureAnnotation = "id", 
             fontcolor.feature = "darkblue", just.feature = 'above') 
  
  #plotTracks(dtrack)
  
}


