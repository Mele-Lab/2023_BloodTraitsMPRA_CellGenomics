library(ggplot2)
library(reshape)
library(dplyr)
library(ggrepel)

### Read GENCODE annotation ####
AnnotationFile <- "gencode.v26.GRCh38.genes.bed"
gene_annotation_path <- "~/marenostrum/Data/GTEx_v8/GeneAnnotation/"

Annotation <- read.delim(paste0(gene_annotation_path,AnnotationFile),header = F)
colnames(Annotation) <- c("chr","start","end","strand","feature","Ensembl_ID","gene_name",
                          "Biotype","source")

filtered_annotations <- Annotation[Annotation$Biotype %in% c("protein_coding","lincRNA"),]

gencode_v26 <- table(filtered_annotations$chr, filtered_annotations$Biotype)
write.table(gencode_v26, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/gencode_v26_pc_lnc.txt',
            sep = '\t')

order = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11',
          'chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21',
          'chr22','chrX','chrY','chrM')
filtered_annotations$chr <- factor(filtered_annotations$chr,levels = order)

ggplot(filtered_annotations, aes(x=factor(chr),fill=Biotype)) +
  geom_bar(width=1, color="white") +
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  scale_fill_brewer(palette="Set1")


### READ GWAS info ####
sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)

sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

only_sentine_traitl <- read.table('~/Desktop/Escritorio - MacBook Pro de Winona - 1/PhD/sentinel_snp.txt', sep = '\t', header = T)
only_sentine_traitl <- only_sentine_traitl[,c('Chr','rsID','Trait')]
only_sentine_traitl <- only_sentine_traitl[!is.na(only_sentine_traitl$Chr),]
only_sentine_traitl <- only_sentine_traitl[!only_sentine_traitl$rsID %in% sentinel_snps$sentinel,]

only_sentinel <- sentinel_snps[,c('sentinel','Chr','Trait')] %>% distinct()

#### How many SNPs per sentinel ####
number_snps <- as.data.frame(table(sentinel_snps$sentinel))
ggplot(number_snps, aes(x=Freq)) + 
  geom_density(adjust = 1/3, color='darkgrey')+
  xlab('N of LD SNPs')+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))

hist(number_snps$Freq, col = "steelblue", frame = FALSE,
     breaks = 60) + abline(v=median(number_snps$Freq), col="red")

# ggplot(number_snps, aes(x=Freq)) + 
#   geom_density(adjust = 1/3, color='darkgrey')+
#   xlab('N of LD SNPs')+
#   theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))

# Basic piechart
order = c('1','2','3','4','5','6','7','8','9','10','11','12','14','15','16','17','18','19','20','22','X')
sentinel_snps$Chr <- factor(sentinel_snps$Chr,levels = order)

ggplot(only_sentinel, aes(x=factor(Chr),fill=Trait)) +
  geom_bar(width=1, color="white") +
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme(legend.position="bottom") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  scale_fill_brewer(palette="Set2")

only_sentinel.os = only_sentinel %>% 
  group_by(Trait) %>% 
  count() %>% 
  ungroup()%>% 
  arrange(desc(Trait)) %>%
  mutate(percentage = round(freq/sum(freq),4)*100,
         lab.pos = cumsum(percentage)-.5*percentage)

traits_df <- as.data.frame(table(only_sentinel$Trait))

# Compute the position of labels
traits_df <- traits_df %>% 
  arrange(desc(Var1)) %>%
  mutate(prop = Freq / sum(traits_df$Freq) *100) %>%
  mutate(ypos = cumsum(prop)- 0.6*prop )

ggplot(traits_df, aes(x="", y=prop, fill=Var1)) +
  geom_bar(width=1, color="white", stat = 'identity') +
  coord_polar("y")+
  geom_text(aes(y = ypos, label = Freq), color = "white", size=6) +
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme(legend.position="bottom") +
  theme_void() + 
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_brewer(palette="Set2")


chr_1_dbp <- nrow(sentinel_snps[sentinel_snps$Chr==1 & sentinel_snps$Trait =='DBP',])
chr_1_other <- nrow(sentinel_snps[sentinel_snps$Chr==1 & sentinel_snps$Trait !='DBP',])
other_dpb <- nrow(sentinel_snps[sentinel_snps$Chr!=1 & sentinel_snps$Trait =='DBP',])
other_other <- nrow(sentinel_snps[sentinel_snps$Chr!=1 & sentinel_snps$Trait !='DBP',])
fisher.test(matrix(c(chr_1_dbp,other_dpb,chr_1_other,other_other),nrow=2,ncol=2))

## enrichment ####

get_col <- function(fdr){
  if(abs(fdr) < -log10(0.05)){
    col <- "light grey"
  }else{
    col <- col_numeric(rev(brewer.pal(11,"RdBu")),
                       domain = seq(-3,
                                    3,length.out = 11))(fdr)
  }
}

# enrichment of DSE in particular types of AS events
my_fisher <- function(chr, trait){
  #             DS      Not DS
  # type
  # other_types
  print(chr)
  print(trait)
  chr_counts <- table(only_sentinel[only_sentinel$Chr == chr, "Trait"])
  not_chr_counts <- table(only_sentinel[only_sentinel$Chr != chr, "Trait"])
  identical(sum(chr_counts)+sum(not_chr_counts),nrow(only_sentinel))
  x11 <- chr_counts[trait]
  x12 <- sum(chr_counts[names(chr_counts)!=trait])
  x21 <- not_chr_counts[trait]
  x22 <- sum(not_chr_counts[names(not_chr_counts)!=trait])
  m <- matrix(c(x11, x12, x21, x22), 2,2, byrow = T)
  print(m)
  m[is.na(m)] <- 0
  rownames(m) <- c("trait", "other_traits")
  colnames(m) <- c("chr","other_chr")
  f <- fisher.test(m)
  print(f)
  return(list("f" = f,
              "overlap" = x11))
}
# Two-tailed Fisher test
types <- c('DBP','PP','SBP')
chrom <- unique(only_sentinel$Chr)
fisher_results <- lapply(c(types), function(e)
  lapply(chrom, function(chr) my_fisher(chr,e)))
names(fisher_results) <- types
for(e in types){names(fisher_results[[e]]) <- chrom}
# Plot
par(mfrow=c(2,2), mar=c(5,7,4,4), mgp=c(4,1,0))
odds_ratio <- list()
adj.P.Val <- list()
for(i in 1:3){
  e <- types[i]
  print(e)
  # odds ratio
  odds_ratio[[e]] <- sapply(chrom, function(chr) fisher_results[[e]][[chr]][['f']]$estimate )
  # adj.P.Val
  adj.P.Val[[e]] <- p.adjust(sapply(chrom, function(chr) fisher_results[[e]][[chr]][['f']]$p.value), method = "BH")
  print(paste0(e,": ", sum(adj.P.Val[[e]]<0.05)))
  print(paste0("Depleted: ", sum(odds_ratio[[e]][which(adj.P.Val[[e]]<0.05)]<1))) #  depleted
  print(paste0("Enriched: ", sum(odds_ratio[[e]][which(adj.P.Val[[e]]<0.05)]>1))) # enriched
  # fdr
  fdr <- -log10(adj.P.Val[[e]])
  fdr[fdr>3] <- 3 # for nice color scale, min adj.P.Val 0.01
  fdr[odds_ratio[[e]]<1] <- -fdr[odds_ratio[[e]]<1] # if depleted blue
  # odds_ratio[[e]][odds_ratio[[e]]<1] <- 1/odds_ratio[[e]][odds_ratio[[e]]<1] # so point size is equal whhether enriched or depleted
  # odds_ratio[[e]][odds_ratio[[e]] == Inf] <- 50
  # # Plot
  # points(1:length(chrom),
  #        rep(4-i,length(chrom)),
  #        cex = sqrt(odds_ratio[[e]]/pi)*2,
  #        col = sapply(fdr, function(i) get_col(i)),
  #        pch =16
  # )
  # axis(1, at=axTicks(1), labels=prettyNum(axTicks(1), big.mark = ","))
  # # round(mean(d[,1]))
  # # round(median(d[,1]))
  # # abline(v=round(mean(d[,1])), lty = 2)
  names(odds_ratio[[e]]) <- chrom
  odds_ratio[[e]] <- as.data.frame(odds_ratio[[e]])
  
  order = c('1','2','3','4','5','6','7','8','9','10','11','12','14','15','16','17','18','19','20','22','X')
  odds_ratio[[e]]$chrom <- rownames(odds_ratio[[e]])
  odds_ratio[[e]]$chrom <- factor(odds_ratio[[e]]$chrom,levels = order)
  colnames(odds_ratio[[e]]) <- c('odds_ratio','chrom')
  odds_ratio[[e]]$pval <- adj.P.Val[[e]]
  
  ggplot(odds_ratio[[e]], aes(x=odds_ratio, y=chrom)) + geom_bar(stat = 'identity')
  
}

threshold_OE <- odds_ratio$DBP$pval < 0.05  & odds_ratio$DBP$odds_ratio > 1
length(which(threshold_OE))
odds_ratio$DBP$thr <- threshold_OE 

ggplot(odds_ratio$DBP, aes(x=odds_ratio, y=chrom)) + 
  geom_bar(aes(fill=thr),stat = 'identity') +
  ggtitle('DBP enrichment') + 
  scale_fill_manual(values = c('grey','Dark red'))+
  coord_flip() + 
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') + 
  geom_vline(aes(xintercept=1), linetype="dotted")


threshold_OE <- odds_ratio$PP$pval < 0.05  & odds_ratio$PP$odds_ratio > 1
length(which(threshold_OE))
odds_ratio$PP$thr <- threshold_OE
ggplot(odds_ratio$PP, aes(x=odds_ratio, y=chrom)) +
  geom_bar(aes(fill=thr),stat = 'identity') +
  ggtitle('PP enrichment') + 
  scale_fill_manual(values = c('grey','Dark red'))+
  coord_flip() + 
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') + 
  geom_vline(aes(xintercept=1), linetype="dotted")

threshold_OE <- odds_ratio$SBP$pval < 0.05  & odds_ratio$SBP$odds_ratio > 1
length(which(threshold_OE))
odds_ratio$SBP$thr <- threshold_OE
ggplot(odds_ratio$SBP, aes(x=odds_ratio, y=chrom)) + 
  geom_bar(aes(fill=thr),stat = 'identity') +
  ggtitle('SBP enrichment') + 
  scale_fill_manual(values = c('grey','Dark red'))+
  coord_flip() + 
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'), 
        legend.position = 'none') + 
  geom_vline(aes(xintercept=1), linetype="dotted")

### SNPs in Genomic locations ####
vep_results <- read.table("~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/vep_snp_list_predictor.txt", sep = '\t')
head(vep_results)

vep_basic <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/vep.hs38.canonicalInfo.txt', sep = '\t')

vep_filtered <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/vep.hs38.canonical.oneConsPerVariant.txt', sep = '\t')

## 4608, 4555 mapped
### get unique results ####
vep_results <- vep_results[,c('V1','V2','V4','V5','V6','V7')]
vep_basic <- vep_basic[,c('X.Uploaded_variation','Location','Consequence','IMPACT','SYMBOL','Gene')]
vep_filtered <- vep_filtered[,c('X.Uploaded_variation','Location','Consequence','IMPACT','SYMBOL','Gene')]

vep_results <- vep_results %>% distinct() 
vep_results$effect <- gsub(',.*','',vep_results$V4)

vep_basic <- vep_basic %>% distinct() 
vep_basic$effect <- gsub(',.*','',vep_basic$Consequence)

vep_filtered <- vep_filtered %>% distinct() 
vep_filtered$effect <- gsub(',.*','',vep_filtered$Consequence)

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

vep_results.os = vep_filtered %>% 
  group_by(effect) %>% 
  count() %>% 
  ungroup()%>% 
  arrange(desc(effect)) %>%
  mutate(percentage = round(freq/sum(freq),4)*100,
         lab.pos = cumsum(percentage)-.5*percentage)

ggplot(vep_results.os, aes(x = "", y = percentage,fill=effect)) +
  geom_bar(width=1, color="white",stat = "identity") +
  coord_polar("y", start = 0)+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme(legend.position="left") +
  #geom_label_repel(aes(label = paste(percentage,"%", sep = "")), size=3, show.legend = F) +
  guides(fill = guide_legend(title = "Effect"))+
  #geom_text(aes(y = lab.pos, label = paste(percentage,"%", sep = "")), col = "white") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_manual(values=mycolors)

ggplot(vep_results.os, aes(x = freq,fill=effect)) +
  geom_bar(width=1, color="white",stat = "count") +
  coord_polar("y", start = 0)+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme(legend.position="left") +
  guides(fill = guide_legend(title = "Effect"))+
  #geom_text(aes(y = lab.pos, label = paste(percentage,"%", sep = "")), col = "white") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())+
  scale_fill_manual(values=mycolors)

ggplot(vep_results.os, aes(y = effect,fill=effect)) +
  geom_bar(width=1, color="white",stat = "count") +
  geom_text(aes(label=..count..),stat="count", hjust=0)+
  ylab("")+
  xlab("")+
  ggtitle("")+
  xlim(c(0,4000))+
  theme(legend.position="none") +
  guides(fill = guide_legend(title = "Effect"))+
  #geom_text(aes(y = lab.pos, label = paste(percentage,"%", sep = "")), col = "white") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_manual(values=mycolors)


### Gencode + fantom annotation ####

gencode_fantom <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/all_results_SNPS_overlap_genome.bed', 
                             sep = '\t', header = F,col.names = paste0("V",seq_len(17)), fill = TRUE)

gencode_fantom_filt <- gencode_fantom[,c('V4','V5')]
gencode_fantom_filt <- gencode_fantom_filt %>% distinct() 

gencode_genes = gencode_fantom_filt[gencode_fantom_filt$V5 %in% c('Gene','Intergenic'),] %>% 
  group_by(V5) %>% 
  count() %>% 
  ungroup()%>% 
  arrange(desc(V5)) %>%
  mutate(percentage = round(n/sum(n),4)*100,
         lab.pos = cumsum(percentage)-.5*percentage)

gencode_all = gencode_fantom_filt %>% 
  filter(V5 != "NA") %>% 
  group_by(V5) %>% 
  count() %>% 
  ungroup()%>% 
  arrange(desc(V5)) %>%
  mutate(percentage = round(n/sum(n),4)*100,
         lab.pos = cumsum(percentage)-.5*percentage)

ggplot(gencode_genes, aes(x = "", y = percentage,fill=V5)) +
  geom_bar(width=1, color="white",stat = "identity") +
  coord_polar("y", start = 0)+
  theme_minimal()+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme(legend.position="left") +
  #geom_label_repel(aes(label = paste(percentage,"%", sep = "")), size=3, show.legend = F) +
  #guides(fill = guide_legend(title = "Effect"))+
  geom_text(aes(y = lab.pos, label = paste(percentage,"%", sep = "")), col = "white") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_manual(values=mycolors)

gencode_sp = gencode_fantom_filt[!gencode_fantom_filt$V5 %in% c('Gene','Intergenic'),] %>% 
  group_by(V5) %>% 
  count() %>% 
  ungroup()%>% 
  arrange(desc(V5)) %>%
  mutate(percentage = round(n/sum(n),4)*100,
         lab.pos = cumsum(percentage)-.5*percentage)

ggplot(gencode_sp, aes(x = "", y = percentage,fill=V5)) +
  geom_bar(width=1, color="white",stat = "identity") +
  coord_polar("y", start = 0)+
  theme_minimal()+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme(legend.position="left") +
  #geom_label_repel(aes(label = paste(percentage,"%", sep = "")), size=3, show.legend = F, max.overlaps = 10) +
  #guides(fill = guide_legend(title = "Effect"))+
  #geom_text(aes(y = lab.pos, label = paste(percentage,"%", sep = "")), col = "white") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_brewer(palette ='Set3')


vals <- gencode_sp$n
val_names <- sprintf("%s (%s)", c("TSS", "Intronic", "Exon","Enhancer","CAGE_peak"), gencode_sp$percentage)
names(vals) <- val_names

waffle::waffle(vals/10)

gencode_sp$ymax <- cumsum(gencode_sp$percentage)
gencode_sp$ymin <- c(0, head(gencode_sp$ymax, n=-1))
gencode_sp$labelPosition <- (gencode_sp$ymax + gencode_sp$ymin) / 2 
gencode_sp$label <- paste0(gencode_sp$V5, "\n value: ", gencode_sp$n)

gencode_sp$labelPosition[gencode_sp$V5 == 'TSS'] <- 5
gencode_sp$labelPosition[gencode_sp$V5 == 'Exon'] <- 75
gencode_sp$labelPosition[gencode_sp$V5 == 'Enhancer'] <- 85
gencode_sp$labelPosition[gencode_sp$V5 == 'CAGE_peak'] <- 95

ggplot(gencode_sp, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=V5)) +
  geom_rect() +
  geom_text( x=2, aes(y=labelPosition, label=label), size=3, color='black') + # x here controls label position (inner / outer)
  scale_fill_brewer(palette='Set3') +
  scale_color_brewer(palette='Set3') +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void()


### random SNPs ####

snpmart = useEnsembl(biomart = "snp", dataset="hsapiens_snp")
getBM(attributes = c('refsnp_id','chr_name','allele','chrom_start','chrom_strand'), 
      mart = snpmart)


### Genome average og introns, TSS... ####
#${data}genes_sorted.hg19.bed 
#${data}exon_sorted.hg19.chr.bed 
#${data}intergenic_sorted.hg19.bed 
#${data}gencode.v34.introns.hg19.chr.sorted.bed 
#${fantom}human_permissive_enhancers_phase_1_and_2.sorted.bed 
#${fantom}TSS_human.sorted.bed 
#${fantom}hg19.cage_peak_phase1and2combined_coord.sorted.bed

## genes
genes <- read.table('~/marenostrum/Data/gene_annotation/gencode/release_34/genes_sorted.hg19.bed', sep = '\t')
genes$type <- 'genes'

## exons 
exons <- read.table('~/marenostrum/Data/gene_annotation/gencode/release_34/exon_sorted.hg19.chr.bed', sep = '\t')
exons$type <- 'exons'

## intergenic
intergenic <- read.table('~/marenostrum/Data/gene_annotation/gencode/release_34/intergenic_sorted.hg19.bed', sep = '\t')
intergenic$type <- 'intergenic'

## introns
introns <- read.table('~/marenostrum/Data/gene_annotation/gencode/release_34/gencode.v34.introns.hg19.chr.sorted.bed', sep = '\t')
introns$type <- 'introns'

## enhancers
enhancers <- read.table('~/marenostrum/Data/FANTOM/human_permissive_enhancers_phase_1_and_2.bed.gz', sep = '\t')
enhancers$type <- 'enhancers'

## TSS
TSS <- read.delim('~/marenostrum/Data/FANTOM/TSS_human.bed.gz', sep = '\t', header = F)
TSS$type <- 'TSS'

## CAGE_peak
CAGE_peak <- read.table('~/marenostrum/Data/FANTOM/hg19.cage_peak_phase1and2combined_coord.bed.gz', sep = '\t')
CAGE_peak$type <- 'CAGE_peak'

### combine all data
genes[,c("V4","V5","V6","V7","V8","V9","V10","V11","V12")] <- NA
exons[,c("V4","V5","V6","V7","V8","V9","V10","V11","V12")] <- NA
introns[,c("V4","V5","V6","V7","V8","V9","V10","V11","V12")] <- NA
intergenic[,c("V4","V5","V6","V7","V8","V9","V10","V11","V12")] <- NA
TSS[,c("V4","V5","V6","V7","V8","V9","V10","V11","V12")] <- NA
CAGE_peak[,c("V4","V5","V6","V7","V8","V9","V10","V11","V12")] <- NA

all <- rbind(genes, exons, introns, intergenic, enhancers, TSS, CAGE_peak)

## plot results

gencode_all_genome_filt <- all[all$type %in% c('TSS','introns','exons','enhancers','CAGE_peak'),]
gencode_all_genome = gencode_all_genome_filt %>% 
  filter(type != "NA") %>% 
  group_by(type) %>% 
  count() %>% 
  ungroup()%>% 
  arrange(desc(type)) %>%
  mutate(percentage = round(n/sum(n),4)*100,
         lab.pos = cumsum(percentage)-.5*percentage)


ggplot(gencode_all_genome, aes(x = "", y = percentage,fill=type)) +
  geom_bar(width=1, color="white",stat = "identity") +
  coord_polar("y", start = 0)+
  theme_minimal()+
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme(legend.position="left") +
  #geom_label_repel(aes(label = paste(percentage,"%", sep = "")), size=3, show.legend = F, max.overlaps = 10) +
  #guides(fill = guide_legend(title = "Effect"))+
  #geom_text(aes(y = lab.pos, label = paste(percentage,"%", sep = "")), col = "white") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_brewer(palette ='Set3')


vals <- gencode_all_genome$n
val_names <- sprintf("%s (%s)", c("TSS", "Intronic", "Exon","Enhancer","CAGE_peak"), gencode_all_genome$percentage)
names(vals) <- val_names

waffle::waffle(vals/1000)

gencode_sp$ymax <- cumsum(gencode_sp$percentage)
gencode_sp$ymin <- c(0, head(gencode_sp$ymax, n=-1))
gencode_sp$labelPosition <- (gencode_sp$ymax + gencode_sp$ymin) / 2 
gencode_sp$label <- paste0(gencode_sp$V5, "\n value: ", gencode_sp$n)

gencode_sp$labelPosition[gencode_sp$V5 == 'TSS'] <- 5
gencode_sp$labelPosition[gencode_sp$V5 == 'Exon'] <- 75
gencode_sp$labelPosition[gencode_sp$V5 == 'Enhancer'] <- 85
gencode_sp$labelPosition[gencode_sp$V5 == 'CAGE_peak'] <- 95

ggplot(gencode_sp, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=V5)) +
  geom_rect() +
  geom_text( x=2, aes(y=labelPosition, label=label), size=3, color='black') + # x here controls label position (inner / outer)
  scale_fill_brewer(palette='Set3') +
  scale_color_brewer(palette='Set3') +
  coord_polar(theta="y") +
  xlim(c(-1, 4)) +
  theme_void()

### Plot random sets of SNPs distributions ####
## set1
vep_1 <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/SNPsnap_test_run_hypertension/vep_set1/vep_set1_results.txt', sep = '\t')

## set2
vep_2 <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/SNPsnap_test_run_hypertension/vep_set2/vep_set2_results.txt', sep = '\t')

## set3
vep_3 <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/SNPsnap_test_run_hypertension/vep_snp3/vep_snp3_results.txt', sep = '\t')

## set4
vep_4 <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/SNPsnap_test_run_hypertension/vep_set4/vep_set4_results.txt', sep = '\t')

## set5
vep_5 <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/SNPsnap_test_run_hypertension/vep_set5/vep_set5_results.txt', sep = '\t')

### RANDOM SET OF SNPS
random <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/vep_random_dbsnp.txt', sep = '\t')

## 4608, 4555 mapped
vep_filtered$set <- 'Original'
vep_filtered <- vep_filtered[,c('X.Uploaded_variation','Location','Consequence','IMPACT','SYMBOL','Gene','set')]
vep_1$set <- 'Set1'
vep_1 <- vep_1[,c('X.Uploaded_variation','Location','Consequence','IMPACT','SYMBOL','Gene','set')]
vep_2$set <- 'Set2'
vep_2 <- vep_2[,c('X.Uploaded_variation','Location','Consequence','IMPACT','SYMBOL','Gene','set')]
vep_3$set <- 'Set3'
vep_3 <- vep_3[,c('X.Uploaded_variation','Location','Consequence','IMPACT','SYMBOL','Gene','set')]
vep_4$set <- 'Set4'
vep_4 <- vep_4[,c('X.Uploaded_variation','Location','Consequence','IMPACT','SYMBOL','Gene','set')]
vep_5$set <- 'Set5'
vep_5 <- vep_5[,c('X.Uploaded_variation','Location','Consequence','IMPACT','SYMBOL','Gene','set')]
random$set <- 'Random'
random <- random[,c('X.Uploaded_variation','Location','Consequence','IMPACT','SYMBOL','Gene','set')]

vep_results <- rbind(vep_filtered, vep_1, vep_2, vep_3, vep_4, vep_5)
vep_results <- rbind(vep_filtered, random)

### get unique results ####
vep_results$effect <- gsub(',.*','',vep_results$Consequence)

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

library(reshape2)
table_vep <- as.data.frame(table(vep_results$effect, vep_results$set))

ggplot(table_vep, aes(x=Var1, y = Freq,fill=Var2)) +
  geom_bar(width=1, color="white",stat = "identity",position="dodge") +
  #geom_text(aes(label=Freq),stat="identity", hjust=0, position = position_dodge(width= 1))+
  #facet_wrap(~Var2, scale = "fixed")+
  #stat_compare_means(aes(group = Var2), label = "p.format", label.y = 3600) +
  geom_signif(tip_length = 0) +
  ylab("")+
  xlab("")+
  ggtitle("")+
  ylim(c(0,4000))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  guides(fill = guide_legend(title = "Set"))+
  #geom_text(aes(y = lab.pos, label = paste(percentage,"%", sep = "")), col = "white") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_manual(values=mycolors)


m <- matrix(c(original_intron, random_intron, original_other, random_other), 2,2, byrow = T)
print(m)
m[is.na(m)] <- 0
rownames(m) <- c("Intron", "Other")
colnames(m) <- c("Original","Random")
f <- fisher.test(m)
print(f)

#### Enrichment of distribution ####

get_col <- function(fdr){
  if(abs(fdr) < -log10(0.05)){
    col <- "light grey"
  }else{
    col <- col_numeric(rev(brewer.pal(11,"RdBu")),
                       domain = seq(-3,
                                    3,length.out = 11))(fdr)
  }
}

# enrichment of DSE in particular types of AS events
my_fisher <- function(type){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  type_counts_original <- table_vep[table_vep$Var1 == type & table_vep$Var2 == 'Original',]$Freq
  other_type_counts_original <- sum(table_vep[table_vep$Var1 != type & table_vep$Var2 == 'Original',]$Freq)
  type_counts_random <- table_vep[table_vep$Var1 == type & table_vep$Var2 == 'Random',]$Freq
  other_type_counts_random <- sum(table_vep[table_vep$Var1 != type & table_vep$Var2 == 'Random',]$Freq)

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
chrom <- unique(table_vep$Var1)
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

### LDlink ####
### token = ac7862288f62

LDproxy_batch(snp = unique(sentinel_snps$sentinel), 
              token = "ac7862288f62"
)

snps_details <- lapply(unique(sentinel_snps$sentinel)[1:79], function(i) {
  print(i)
  read.table(paste0(i,'.txt'), sep='\t', header = T)
  })
names(snps_details) <- unique(sentinel_snps$sentinel)[1:79]

snps_details_2 <- lapply(unique(sentinel_snps$sentinel)[81:135], function(i) {
  print(i)
  read.table(paste0(i,'.txt'), sep='\t', header = T)
})
names(snps_details_2) <- unique(sentinel_snps$sentinel)[81:135]

snps_details_df <- do.call('rbind.data.frame', snps_details)
snps_details_df_2 <- do.call('rbind.data.frame', snps_details_2)

all_snps <- rbind(snps_details_df, snps_details_df_2)
all_snps <- all_snps[all_snps$R2>0.8,]
all_snps$sentinel <- gsub('\\..*','',rownames(all_snps))

all_snps$RS_Number[all_snps$RS_Number %in% sentinel_snps$snp]

snps_sentinel_ours <- all_snps[all_snps$RS_Number %in% sentinel_snps$snp,]
snps_sentinel_ours$sentinel <- gsub('\\..*','',rownames(snps_sentinel_ours))

snps_sentinel_ours <- merge(snps_sentinel_ours, sentinel_snps, by.x = c('RS_Number','sentinel'), by.y = c('snp','sentinel'))
snps_sentinel_ours <- snps_sentinel_ours %>% distinct()
sentinel_snps$snp[! sentinel_snps$snp %in% snps_sentinel_ours$RS_Number]

write.table(snps_sentinel_ours, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/PDLink.results.134Sentinel.txt', 
            sep = '\t', quote = F, row.names = F, col.names = T)

### Get information per sentinel SNP ####
sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)

sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

snps_sentinel_ours <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/PDLink.results.134Sentinel.txt', sep='\t', header = T)

CM_new <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.CM.new_back.005.txt', sep='\t', header = T)
VSMC_new <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.VSMC.new_back.005.txt', sep='\t', header = T)

library(dplyr)
library(tidyr)
vals_significance_CM <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/CM_vals.significance.5.txt', sep = '\t', header = T)
head(vals_significance_CM)
vals_significance_CM <- vals_significance_CM %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance_CM)

vals_significance <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/VSMC_vals.significance.4.txt', sep = '\t', header = T)
head(vals_significance)
vals_significance <- vals_significance %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance)

active_snps_vsmc <- unique(vals_significance$snp_info[vals_significance$sig == 'sig' & vals_significance$tile_type != 'RANDOM'])
active_snps_cm <- unique(vals_significance_CM$snp_info[vals_significance_CM$sig == 'sig' & vals_significance_CM$tile_type != 'RANDOM'])

snps <- vals_significance[vals_significance$tile_type %in% c('WILDTYPE_BUT_HAS_SNP','WILDTYPE_SNP_INDIV'),c('pos','snp_info','chrom','dupe_info','info')] %>% distinct()
snps <- snps %>%
  separate(info, c("chr", "start","ref","alt"), ":")

CM_new <- CM_new[CM_new$tile_type %in% c('WILDTYPE_SNP_INDIV','WILDTYPE_BUT_HAS_SNP'),]
VSMC_new <- VSMC_new[VSMC_new$tile_type %in% c('WILDTYPE_SNP_INDIV','WILDTYPE_BUT_HAS_SNP'),]

CM_new <- merge(CM_new, snps[,c('snp_info','dupe_info')])
VSMC_new <- merge(VSMC_new, snps[,c('snp_info','dupe_info')])

CM_new <- merge(CM_new, sentinel_snps, by.x = 'snp_info', by.y = 'snp', all.y = T)
VSMC_new <- merge(VSMC_new, sentinel_snps, by.x = 'snp_info', by.y = 'snp', all.y = T)

CM_new <- CM_new %>% distinct()
VSMC_new <- VSMC_new %>% distinct()

CM_new$Activity <- 'Not Active'
CM_new$Activity[CM_new$snp_info %in% active_snps_cm] <- 'Active'
VSMC_new$Activity <- 'Not Active'
VSMC_new$Activity[VSMC_new$snp_info %in% active_snps_vsmc] <- 'Active'

CM_new$Functional <- 'Not Functional Variant'
CM_new$Functional[CM_new$fdr_comp < 0.05] <- 'Functional Variant'
VSMC_new$Functional <-  'Not Functional Variant'
VSMC_new$Functional[VSMC_new$fdr_comp < 0.05] <- 'Functional Variant'

write.table(CM_new, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/all_snps_act_diff_CM.txt', 
            quote = F, row.names = F, col.names = T, sep = '\t')
write.table(VSMC_new, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/all_snps_act_diff_VSMC.txt', 
            quote = F, row.names = F, col.names = T, sep = '\t')

### plot some counts ####
counts_dir = "~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/01_counts"
mpranalyze_dir = paste0(counts_dir,"/mpranalyze_files")

CM_counts = read.table(paste0(counts_dir,"/CM__all_counts_final.5.txt"), sep="\t", header = T)
head(CM_counts)

### VSMC
VSMC_counts = read.table(paste0(counts_dir,"/VSMC__all_counts_final.4.txt"), sep="\t", header = T)
head(VSMC_counts)

### index 
index <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/Hypertension__pooled.index.txt', sep = '\t', header = T)
head(index)

colnames(VSMC_counts) = c("barcode", "dna_1", "VSMC_rep1", "VSMC_rep2", "VSMC_rep3","VSMC_rep6")
colnames(CM_counts)  = c("barcode", "dna_1", "CM_rep1", "CM_rep2", "CM_rep3", "CM_rep4","CM_rep5")

df = merge(VSMC_counts,index[,c("barcode", "element", "tile_type", "tile_id","snp","name")], by="barcode")
df <- df %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(df)

df_CM = merge(CM_counts,index[,c("barcode", "element", "tile_type", "tile_id","snp","name")], by="barcode")
df_CM <- df_CM %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(df_CM)

#df["barc_id"] = df.apply(get_barc_id, axis=1).astype(int)

df['median_rep'] = rowMeans(df[,c("VSMC_rep1", "VSMC_rep2", "VSMC_rep3","VSMC_rep6")], na.rm = T)
df['median_rep'] = rowMeans(df[,c("CM_rep1", "CM_rep2", "CM_rep3", "CM_rep4")], na.rm = T)

df$dupe <- gsub('\\..*','',df$tile_id)

for (snp in c('2075')) {
  subs <- df[df$dupe == snp,]
  snp <- unique(subs$snp[subs$snp != 'none'])
  print(ggplot(subs, aes(x=tile_type,y=median_rep, color=tile_type)) + geom_boxplot(fill='white') +
          geom_jitter(aes(fill=tile_type))+
          ylab('median counts/barcode')+
          xlab('')+
          ggtitle(snp)+
          theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
          scale_fill_brewer(palette ='Set2') +
          scale_color_brewer(palette ='Set2'))
}

active_sentinel <- as.data.frame(table(CM_new$sentinel, CM_new$Activity))
active_sentinel <- active_sentinel[active_sentinel$Var2 == 'Active',]
colnames(active_sentinel) <- c('sentinel','Active_info','ActiveSNPs')

functional_sentinel <- as.data.frame(table(CM_new$sentinel, CM_new$Functional))
functional_sentinel <- functional_sentinel[functional_sentinel$Var2 == 'Functional Variant',]
colnames(functional_sentinel) <- c('sentinel','Functional_info','FunctionalSNPs')

number_sentinel <- as.data.frame(table(CM_new$sentinel))
colnames(number_sentinel) <- c('sentinel','NumberSNPs')

CM_info <- merge(number_sentinel, merge(active_sentinel[,c('sentinel','ActiveSNPs')],
                                        functional_sentinel[,c('sentinel','FunctionalSNPs')]))

active_sentinel_vsmc <- as.data.frame(table(VSMC_new$sentinel, VSMC_new$Activity))
active_sentinel_vsmc <- active_sentinel_vsmc[active_sentinel_vsmc$Var2 == 'Active',]
colnames(active_sentinel_vsmc) <- c('sentinel','Active_info','ActiveSNPs')

functional_sentinel_vsmc <- as.data.frame(table(VSMC_new$sentinel, VSMC_new$Functional))
functional_sentinel_vsmc <- functional_sentinel_vsmc[functional_sentinel_vsmc$Var2 == 'Functional Variant',]
colnames(functional_sentinel_vsmc) <- c('sentinel','Functional_info','FunctionalSNPs')

VSMC_info <- merge(number_sentinel, merge(active_sentinel_vsmc[,c('sentinel','ActiveSNPs')],
                                        functional_sentinel_vsmc[,c('sentinel','FunctionalSNPs')]))

write.table(CM_info, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/sentinel_summary_Act_funct_CM.txt', 
            quote = F, row.names = F, col.names = T, sep = '\t')
write.table(VSMC_info, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/sentinel_summary_Act_funct_VSMC.txt', 
            quote = F, row.names = F, col.names = T, sep = '\t')

### correlate N SNPs active and N SNPs per locus ####
library(ggpubr)
sp <- ggscatter(CM_info, x = "NumberSNPs", y = "FunctionalSNPs",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 3, label.y = 150)


### correlate FC with N of SNPs ####
head(CM_new)
head(CM_info)

CM_new_info <- merge(CM_new, CM_info)
VSMC_new_info <- merge(VSMC_new, VSMC_info)
CM_new_info$abs_FC <- abs(CM_new_info$logFC_comp)
VSMC_new_info$abs_FC <- abs(VSMC_new_info$logFC_comp)

sp <- ggscatter(VSMC_new_info, x = "NumberSNPs", y = "abs_FC",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 3, label.y = 5)

#### Correlate logFC with N. TF disrupted ####
CM_Tf <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header = T)
VSMC_Tf <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header = T)

N_TF_WT <- strsplit(CM_Tf$TF_WT, split = ",")
names(N_TF_WT) <- CM_Tf$snp_info
N_TF_WT <- lapply(N_TF_WT, length )
N_TF_WT <- unlist(N_TF_WT)
CM_Tf$N_TF_WT <- N_TF_WT 

N_TF_SNP <- strsplit(CM_Tf$TF_SNP, split = ",")
names(N_TF_SNP) <- CM_Tf$snp_info
N_TF_SNP <- lapply(N_TF_SNP, length )
N_TF_SNP <- unlist(N_TF_SNP)
CM_Tf$N_TF_SNP <- N_TF_SNP 

CM_Tf$N_diff_TF <- CM_Tf$N_TF_WT + CM_Tf$N_TF_SNP

N_TF_WT <- strsplit(VSMC_Tf$TF_WT, split = ",")
names(N_TF_WT) <- VSMC_Tf$snp_info
N_TF_WT <- lapply(N_TF_WT, length )
N_TF_WT <- unlist(N_TF_WT)
VSMC_Tf$N_TF_WT <- N_TF_WT 

N_TF_SNP <- strsplit(VSMC_Tf$TF_SNP, split = ",")
names(N_TF_SNP) <- VSMC_Tf$snp_info
N_TF_SNP <- lapply(N_TF_SNP, length )
N_TF_SNP <- unlist(N_TF_SNP)
VSMC_Tf$N_TF_SNP <- N_TF_SNP 

VSMC_Tf$N_diff_TF <- VSMC_Tf$N_TF_WT + VSMC_Tf$N_TF_SNP

CM_Tf$abs_FC <- abs(CM_Tf$logFC_comp)
VSMC_Tf$abs_FC <- abs(VSMC_Tf$logFC_comp)

sp <- ggscatter(CM_Tf, x = "N_diff_TF", y = "abs_FC",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 3, label.y = 5)

#### Correlate Risk with FC #####
risk_known <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/Nature_Gen_SNPs/152_prevreported_snps.txt', sep='\t', header = T)
risk_new <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/Nature_Gen_SNPs/115_validNovel.txt', sep='\t', header = T)

head(risk_known)
head(risk_new)

risk_all <- rbind(risk_known, risk_new)

LDlink <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/PDLink.results.134Sentinel.txt', sep='\t', header = T)
VSMC_new_info_corr <- merge(VSMC_new_info, LDlink, all.x = T)
CM_new_info_corr <- merge(CM_new_info, LDlink, all.x = T)

head(VSMC_new_info)
head(risk_known)
VSMC_new_info_corr_risk <- merge(VSMC_new_info_corr, risk_all, by.x = 'sentinel', by.y = 'rsID', all.x = T)
CM_new_info_corr_risk <- merge(CM_new_info_corr, risk_all, by.x = 'sentinel', by.y = 'rsID', all.x = T)

### SBP 
filter <- CM_new_info_corr_risk[CM_new_info_corr_risk$Trait == 'SBP',]
head(filter)
filter$sbp_beta <- gsub(',','\\.',filter$sbp_beta)
filter$sbp_beta <- as.numeric(filter$sbp_beta)
filter$abs_beta <- abs(filter$sbp_beta)

sp <- ggscatter(filter[filter$fdr_comp < 0.05,], x = "abs_beta", y = "abs_FC",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 0.2, label.y = 4)

### DBP 
filter <- CM_new_info_corr_risk[CM_new_info_corr_risk$Trait == 'DBP',]
head(filter)
filter$dbp_beta <- gsub(',','\\.',filter$dbp_beta)
filter$dbp_beta <- as.numeric(filter$dbp_beta)
filter$abs_beta <- abs(filter$dbp_beta)

sp <- ggscatter(filter[filter$fdr_comp < 0.05,], x = "abs_beta", y = "abs_FC",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 0.2, label.y = 4)

### PP
filter <- VSMC_new_info_corr_risk[VSMC_new_info_corr_risk$Trait == 'PP',]
head(filter)
filter$pp_beta <- gsub(',','\\.',filter$pp_beta)
filter$pp_beta <- as.numeric(filter$pp_beta)
filter$abs_beta <- abs(filter$pp_beta)

sp <- ggscatter(filter[filter$fdr_comp < 0.05,], x = "abs_beta", y = "abs_FC",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 0.2, label.y = 4)

### only with sentinel ####
VSMC_new_info_corr_risk <- merge(VSMC_new_info_corr, risk_all, by.x = 'snp_info', by.y = 'rsID', all.x = T)
CM_new_info_corr_risk <- merge(CM_new_info_corr, risk_all, by.x = 'snp_info', by.y = 'rsID', all.x = T)

### SBP 
filter <- VSMC_new_info_corr_risk[VSMC_new_info_corr_risk$Trait == 'DBP',]
head(filter)
filter$dbp_beta <- gsub(',','\\.',filter$dbp_beta)
filter$dbp_beta <- as.numeric(filter$dbp_beta)
filter$abs_beta <- abs(filter$dbp_beta)

sp <- ggscatter(filter, x = "abs_beta", y = "abs_FC",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 0.2, label.y = 4)


write.table(CM_new_info_corr_risk, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/all_snps_info_risk_CM.txt', 
            quote = F, row.names = F, col.names = T, sep = '\t')
write.table(VSMC_new_info_corr_risk, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/all_snps_info_risk_VSMC.txt', 
            quote = F, row.names = F, col.names = T, sep = '\t')

### correlate Risk with n of SNPs ####
filter <- VSMC_new_info_corr_risk[VSMC_new_info_corr_risk$Trait == 'DBP',]
head(filter)
filter$dbp_beta <- gsub(',','\\.',filter$dbp_beta)
filter$dbp_beta <- as.numeric(filter$dbp_beta)
filter$abs_beta <- abs(filter$dbp_beta)

sp <- ggscatter(filter, x = "abs_beta", y = "NumberSNPs",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 0.2, label.y = 150)

#### SBP
filter <- VSMC_new_info_corr_risk[VSMC_new_info_corr_risk$Trait == 'SBP',]
head(filter)
filter$sbp_beta <- gsub(',','\\.',filter$sbp_beta)
filter$sbp_beta <- as.numeric(filter$sbp_beta)
filter$abs_beta <- abs(filter$sbp_beta)

sp <- ggscatter(filter, x = "abs_beta", y = "NumberSNPs",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 0.2, label.y = 150)

### PP 
filter <- CM_new_info_corr_risk[CM_new_info_corr_risk$Trait == 'DBP',]
head(filter)
filter$dbp_beta <- gsub(',','\\.',filter$dbp_beta)
filter$dbp_beta <- as.numeric(filter$dbp_beta)
filter$abs_beta <- abs(filter$dbp_beta)

sp <- ggscatter(filter, x = "abs_beta", y = "NumberSNPs",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
sp + stat_cor(method = "pearson", label.x = 0.2, label.y = 150)

#### enrichment of genomic location ChiPSeeker #####
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)

bed_oligos <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/bed_oligos_testtiles.bed')
bed_snps <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/snp_coord_master.hg19.simple.chr.sortedplus1.bed')
colnames(bed_snps) <- c('chr','start','end','snp')
bed_snps <- GRanges(bed_snps)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnnoList=lapply(list(bed_snps), annotatePeak, TxDb=txdb)
plotAnnoPie(peakAnnoList[[1]])

random <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/vep_random_dbsnp.txt', sep = '\t')
random$set <- 'Random'
random <- random[,c('X.Uploaded_variation','Location','Consequence','IMPACT','SYMBOL','Gene','set')]

bed_random <- random[,c('Location','X.Uploaded_variation')]
bed_random <- bed_random %>%
  separate(Location, c("chr", "loc"), ":")
bed_random <- bed_random %>%
  separate(loc, c("start", "end"), "-")
bed_random$chr <- paste0('chr',bed_random$chr)
bed_random$start <- as.numeric(bed_random$start) - 1
bed_random <- GRanges(bed_random)

peakAnnoList_random=lapply(list(bed_random), annotatePeak, TxDb=txdb)
plotAnnoPie(peakAnnoList_random[[1]])

tile_snps <- as.data.frame(peakAnnoList[[1]])
random_snps <- as.data.frame(peakAnnoList_random[[1]])

#### Get enrichment ####

get_col <- function(fdr){
  if(abs(fdr) < -log10(0.05)){
    col <- "light grey"
  }else{
    col <- col_numeric(rev(brewer.pal(11,"RdBu")),
                       domain = seq(-3,
                                    3,length.out = 11))(fdr)
  }
}

# enrichment of DSE in particular types of AS events
my_fisher <- function(type){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  type_counts_original <- nrow(tile_snps[grep(type,tile_snps$annotation),])
  other_type_counts_original <- nrow(tile_snps) - type_counts_original
  type_counts_random <- nrow(random_snps[grep(type,random_snps$annotation),])
  other_type_counts_random <- nrow(random_snps) - type_counts_random
  
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
chrom <- c('Promoter','UTR','Intron',"3' UTR", "5' UTR","Distal Intergenic","Downstream",
           "Exon")
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



threshold_OE <- odds_ratio$pval < 0.05 
# & odds_ratio$odds_ratio > 1
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

#### compare enrichment of genomic locations differential vs not differential #####
CM_new <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.CM.new_back.005.txt', sep='\t', header = T)
VSMC_new <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.VSMC.new_back.005.txt', sep='\t', header = T)

diff_CM <- CM_results$info_snp
diff_VSMC <- VSMC_results$info_snp

vals_significance_cm <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/CM_vals.significance.5.txt', sep = '\t', header = T)
head(vals_significance_cm)
vals_significance_cm <- vals_significance_cm %>%
  separate(name, c("pos", "snp_info","info"), "__")
active_snps_cm <- unique(vals_significance_cm$info[vals_significance_cm$CM_padj<0.05])

vals_significance_vsmc <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/VSMC_vals.significance.4.txt', sep = '\t', header = T)
head(vals_significance_vsmc)
vals_significance_vsmc <- vals_significance_vsmc %>%
  separate(name, c("pos", "snp_info","info"), "__")
active_snps_vsmc <- unique(vals_significance_vsmc$info[vals_significance_vsmc$VSMC_padj<0.05])


#### enrichment #####
my_fisher <- function(type){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  type_data <- tile_snps[grep(type,tile_snps$annotation),]
  other <- tile_snps[!tile_snps$snp %in% type_data$snp,]
  type_counts_diff <-  nrow(type_data[type_data$snp %in% active_snps_cm,])
  other_type_counts_diff <- nrow(other[other$snp %in% active_snps_cm,])
  type_counts_notdiff<- nrow(type_data[!type_data$snp %in% active_snps_cm,])
  other_type_counts_notdiff <- nrow(other[!other$snp %in% active_snps_cm,])
  
  m <- matrix(c(type_counts_diff, type_counts_notdiff, other_type_counts_diff, other_type_counts_notdiff), 2,2, byrow = T)
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
chrom <- c('Promoter','UTR','Intron',"3' UTR", "5' UTR","Distal Intergenic","Downstream",
           "Exon")
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



threshold_OE <- odds_ratio$pval < 0.05 
# & odds_ratio$odds_ratio > 1
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
