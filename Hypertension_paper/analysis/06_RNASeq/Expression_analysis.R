library(ggplot2)
library(ggpubr)
library(dplyr)
library(hrbrthemes)
library(viridis)

### Parse RNAseq data for TF expression ###

### input files #####
VSMC_rep_1 <- read.table('../../analysis/RNASeq/feature_counts/V1_21d/V1_21d.counts.txt', sep='\t', header=TRUE, skip = 1)

VSMC_all <- read.table('../../analysis/RNASeq/VSMC_all.counts.voom.txt', header=T)
HEK_all <- read.table('../../analysis/RNASeq/HEK_all.counts.voom.txt',header=T)

#### publicly available data for CMs
CM_1 <- read.table('marenostrum/Data/Hypertension/CM/RNASeq/GSM3262978_RZY637_RNA_D80_1-chrM.rpkm.gz', sep = '\t', header = T)
CM_2 <- read.table('marenostrum/Data/Hypertension/CM/RNASeq/GSM3262979_RZY643_RNA_D80_2-chrM.rpkm.gz', sep = '\t', header = T)

rownames(CM_1) <- CM_1$Geneid
CM_1 <- CM_1[,c('Length','bam.RZY637_RNA_D80_1.nodup.bam')]
colnames(CM_1) <- c('Length','RPKM')
CM_1$Geneid <- rownames(CM_1)

rownames(CM_2) <- CM_2$Geneid
CM_2 <- CM_2[,c('Length','bam.RZY643_RNA_D80_2.nodup.bam')]
colnames(CM_2) <- c('Length','RPKM')
CM_2$Geneid <- rownames(CM_2)

all_counts_CM <- merge(CM_1, CM_2, by=c('Geneid','Length'))
rownames(all_counts_CM) <- all_counts_CM$Geneid
all_counts_CM <- all_counts_CM[,c('Length','RPKM.x','RPKM.y')]
colnames(all_counts_CM) <- c('Length','RPKM_CM1','RPKM_CM2')


### rpkm to tpm ####
TPM_CM = apply( subset(all_counts_CM, select = c(-Length)), 2, 
       function(x) ((x)/sum(as.numeric(x)))*10^6)


VSMC_all <- merge(VSMC_all, VSMC_rep_1[,c('Geneid','Length')], by.x='gene_id',by.y='Geneid')
rownames(VSMC_all) <- VSMC_all$gene_id
VSMC_all <- VSMC_all[,c('Length','rep1','rep2')]

### calculate TPM ####
rpk <- apply( subset(VSMC_all, select = c(-Length)), 2, 
              function(x) x/(VSMC_all$Length/1000))
#normalize by the sample size using rpk values
tpm_VSMC <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)


HEK_all <- merge(HEK_all, VSMC_rep_1[,c('Geneid','Length')], by.x='gene_id',by.y='Geneid')
rownames(HEK_all) <- HEK_all$gene_id
HEK_all <- HEK_all[,c('Length','rep1','rep2')]

### calculate TPM ####
rpk <- apply( subset(HEK_all, select = c(-Length)), 2, 
              function(x) x/(HEK_all$Length/1000))
#normalize by the sample size using rpk values
tpm_HEK <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)


all_counts <- merge(VSMC_all, HEK_all, by = c('gene_id','Length'))
#all_counts_gene_length <- merge(all_counts, VSMC_rep_1[,c('Geneid','Length')], by.x='gene_id',by.y='Geneid')
rownames(all_counts) <- all_counts$gene_id
all_counts_gene_length <- all_counts[,c('VSMC','HEK','Length')]

### calculate TPM ####
rpk <- apply( subset(all_counts_gene_length, select = c(-Length)), 2, 
              function(x) x/(all_counts_gene_length$Length/1000))
#normalize by the sample size using rpk values
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)

# y <- DGEList(counts=PDX_cts, genes = lengths)
# y <- calcNormFactors(y)
# RPKM <- rpkm(y)


### TF to look for ####
tf <- unique(c('ZNF879','FOXC2','FOXA1','ZNF35','FOXC1','FOXL1','FOXF1','FOXJ2','ARID5A','ZNF177',
        'ZNF324','ZNF174','TOPORS','SOX3','LCOR','SIX1','ZNF224','ZNF276','ZNF135','PRDM14',
        'NR2F1','NR2C1','ZNF774','MITF','ZNF586','ZNF177','ZNF707','SIX1','ZNF224','RORA','ZNF445',
        'ZNF329','ZNF283','TP63','ZNF524'))

tpm <- as.data.frame(tpm)
tpm$gene_id <- gsub('\\..*','',rownames(tpm))

### get symbol ####
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

symbol <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=tpm$gene_id,mart= mart)
tpm = merge(tpm,symbol,by.x="gene_id",by.y="ensembl_gene_id", all.x=TRUE)


for (snp in tf) {
  print(ggplot(subs, aes(x=tile_type,y=median_rep, color=tile_type)) + geom_boxplot(fill='white') +
          geom_jitter(aes(fill=tile_type))+
          ylab('median counts/barcode')+
          xlab('')+
          ggtitle(snp)+
          theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
          scale_fill_brewer(palette ='Set2') +
          scale_color_brewer(palette ='Set2'))
}

ggplot(df_subs, aes(x=tile_type,y=median_rep, color=tile_type)) + geom_boxplot(fill='white') +
  geom_jitter(aes(fill=tile_type))+
  ylab('median counts/barcode')+
  xlab('')+
  facet_wrap(~dupe, scales = "free")+
  ggtitle('')+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_brewer(palette ='Set2') +
  scale_color_brewer(palette ='Set2')

library(RColorBrewer)
n <- 32
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colors = sample(col_vector,n)

ggplot(tpm[tpm$hgnc_symbol %in% tf,], aes(y=hgnc_symbol, x=VSMC)) + 
  geom_bar(stat='identity', fill=colors) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  ylab('')+
  xlab('VSMC TPM')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))


#### TF for Phil -----
library(AnnotationDbi)
library(org.Hs.eg.db)
genes_vsmc <- mapIds(org.Hs.eg.db,keys=gsub('\\..*','',rownames(tpm_VSMC)),column="SYMBOL", keytype="ENSEMBL", multiVals="first")
tpm_VSMC <- as.data.frame(tpm_VSMC)
tpm_VSMC$gene_id <- genes_vsmc
  
all_counts <- merge(tpm_VSMC, TPM_CM, by.x = c('gene_id'), by.y=0)



### TF to look for ####
tf <- unique(c('SMARCC1','ULK4','MAP4','ESR1','PDE5A','INSR','CFDP1','CPEB4','CYP17A1','CNNM2',
               'NT5C2','RPEL1','BORCS7'))

all_counts$median_VSMC <- rowMeans(all_counts[,c('rep1','rep2')])
all_counts$median_cM <- rowMeans(all_counts[,c('RPKM_CM1','RPKM_CM2')])

library(reshape)
counts_res <- melt(all_counts, id=c("gene_id","rep1",'rep2','RPKM_CM1','RPKM_CM2'))

ggplot(counts_res[counts_res$gene_id %in% tf,], aes(y=gene_id, x=(value), fill=variable)) + 
  geom_bar(stat='identity',position=position_dodge()) +
  theme_minimal() +
  ylab('')+
  xlab('TPM')+
  scale_fill_manual(values=c('#EC7847','#7C56A1'))

  
#### Check expression of TF ####
human_TF <- read.table('../../data/03_fimo/Kaia_FIMO/curated_motif_map.txt', sep = '\t', header = T)
TPM_CM <- as.data.frame(TPM_CM)
TPM_CM_tf <- TPM_CM[human_TF$gene_name,]
TPM_CM_tf_filt <- TPM_CM_tf[-grep('NA\\.',rownames(TPM_CM_tf)),]

tpm_HEK <- as.data.frame(tpm_HEK)
tpm_HEK_tf <- tpm_HEK[human_TF$gene_id,]
tpm_HEK_tf_filt <- tpm_HEK_tf[-grep('NA\\.',rownames(tpm_HEK_tf)),]

tpm_VSMC <- as.data.frame(tpm_VSMC)
tpm_VSMC_tf <- tpm_VSMC[human_TF$gene_id,]
tpm_VSMC_tf_filt <- tpm_VSMC_tf[-grep('NA\\.',rownames(tpm_VSMC_tf)),]

#### How many TF are expressed on each cell line? ####
View(tpm_VSMC_tf_filt[tpm_VSMC_tf_filt>1,])

nsamples = ncol(tpm_VSMC_tf_filt)
keep     = (rowSums(is.na(tpm_VSMC_tf_filt),na.rm=F) + rowSums(tpm_VSMC_tf_filt<1,na.rm=T) ) <= nsamples/2
tpm_VSMC_tf_filt     = tpm_VSMC_tf_filt[keep,]
tpm_VSMC_tf_filt$geneid <- gsub('\\..*','',rownames(tpm_VSMC_tf_filt))

nsamples = ncol(tpm_HEK_tf_filt)
keep     = (rowSums(is.na(tpm_HEK_tf_filt),na.rm=F) + rowSums(tpm_HEK_tf_filt<1,na.rm=T) ) <= nsamples/2
tpm_HEK_tf_filt     = tpm_HEK_tf_filt[keep,]
tpm_HEK_tf_filt$geneid <- gsub('\\..*','',rownames(tpm_HEK_tf_filt))

nsamples = ncol(TPM_CM_tf_filt)
keep     = (rowSums(is.na(TPM_CM_tf_filt),na.rm=F) + rowSums(TPM_CM_tf_filt<1,na.rm=T) ) <= nsamples/2
TPM_CM_tf_filt     = TPM_CM_tf_filt[keep,]
TPM_CM_tf_filt$geneid <- gsub('\\..*','',rownames(TPM_CM_tf_filt))

tpm_VSMC_tf_filt <- tpm_VSMC_tf_filt %>% distinct()
tpm_HEK_tf_filt <- tpm_HEK_tf_filt %>% distinct()
TPM_CM_tf_filt <- TPM_CM_tf_filt %>% distinct()

names_tf_cm <- unique(human_TF$gene_name[human_TF$gene_name %in% TPM_CM_tf_filt$geneid])
names_tf_hek <- unique(human_TF$gene_name[human_TF$gene_id %in% tpm_HEK_tf_filt$geneid])
names_tf_vsmc <- unique(human_TF$gene_name[human_TF$gene_id %in% tpm_VSMC_tf_filt$geneid])

listInput <- list(CM = names_tf_cm, HEK = names_tf_hek, 
                  VSMC = names_tf_vsmc)

library(UpSetR)
upset(fromList(listInput), order.by = "freq")

# library
library(VennDiagram)

#Make the plot
venn.diagram(
  x = list(names_tf_cm , 
    names_tf_hek , 
    names_tf_vsmc),
  category.names = c("CM" , "HEK" , "VSMC"),
  filename = '../../plots/venn.svg',
  output = TRUE ,
  imagetype="svg",
  # Output features
  height = 100 , 
  width = 100 , 
  resolution = 300,
  compression = "lzw")

if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")

library(ggvenn)
ggvenn(
  listInput, 
  fill_color = c("#999999", "#E69F00", "#56B4E9"),
  stroke_size = 0.6, set_name_size = 4
)


