##### phastcons ------
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
BiocManager::install("GenomicScores")
library(GenomicScores)

#### vcf file -----
sentinel_snps = read.table('../../data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)

sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

sentinel_snps <- sentinel_snps %>%
  separate(change, c('#CHROM','POS','REF','ALT'), ":")
head(sentinel_snps)

vcf <- sentinel_snps
vcf$FILTER <- 'PASS'
vcf <- vcf[,c('#CHROM','POS','snp','REF','ALT','FILTER')]
colnames(vcf) <- c('#CHROM','POS','ID','REF','ALT','FILTER')
write.table(vcf, '../../data/design/snps.vcf',
            sep = '\t', col.names = T, row.names = F, quote = F)

file <- '../../data/design/snps.vcf'
vcf <- readVcf(file)##, "hg19")
seqlevelsStyle(vcf)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevelsStyle(txdb)

#### change nomenclature ----
seqlevelsStyle(vcf) <- seqlevelsStyle(txdb)

loc <- locateVariants(vcf, txdb, AllVariants())
loc[1:3]

table(loc$LOCATION)

#### score -----
library(phastCons100way.UCSC.hg19)
phast <- phastCons100way.UCSC.hg19
class(phast)

loc$PHASTCONS <- score(phast, loc)
loc[1:3]

### plots ----
x <- split(loc$PHASTCONS, loc$LOCATION)
mask <- elementNROWS(x) > 0
boxplot(x[mask], ylab="phastCons score", las=1, cex.axis=1.2, cex.lab=1.5, col="gray")
points(1:length(x[mask])+0.25, sapply(x[mask], mean, na.rm=TRUE), pch=23, bg="black")

#### merge with other data -----
sentinel_snps = read.table('../../data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

vals_significance_cm <- read.table('../../data/02_activs/CM_vals.significance.5.txt', sep = '\t', header = T)
head(vals_significance_cm)
vals_significance_cm <- vals_significance_cm %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance_cm)

vals_significance_vsmc <- read.table('../../data/02_activs/VSMC_vals.significance.4.txt', sep = '\t', header = T)
head(vals_significance_vsmc)
vals_significance_vsmc <- vals_significance_vsmc %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance_vsmc)

mital_names <- read.table('../../data/design/name_Mital_SNPS.txt', sep='\t', header=T)
mital_details <- index$parse_details[index$name %in% mital_names$name]
mital_name <- unique(index$dupe_info[index$name %in% mital_names$name])

active_snps_cm <- unique(vals_significance_cm$snp_info[vals_significance_cm$CM_padj<0.05 ])
active_snps_vsmc <- unique(vals_significance_vsmc$snp_info[vals_significance_vsmc$VSMC_padj<0.05])

CM_results <- read.table('../../analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_ordered.CM.005.new_back.txt', sep='\t', header = T)
VSMC_results <- read.table('../../analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_ordered.VSMC.005.new_back.txt', sep='\t', header = T)

#### parse data ----
head(CM_results)
diff_snps <- unique(CM_results$snp[CM_results$fdr_comp<0.05])
diff_snps_vsmc <- unique(VSMC_results$snp[VSMC_results$fdr_comp<0.05])

sentinel_snps$activeDiff <- 'Not ActiveDiff'
sentinel_snps$activeDiff[sentinel_snps$snp %in% diff_snps_vsmc & sentinel_snps$snp %in% active_snps_vsmc] <- 'ActiveDiff'
sentinel_snps$activeDiffCM <- 'Not ActiveDiff'
sentinel_snps$activeDiffCM[sentinel_snps$snp %in% diff_snps & sentinel_snps$snp %in% active_snps_cm] <- 'ActiveDiff'
sentinel_snps$activeCM <- 'Not Active'
sentinel_snps$activeVSMC <- 'Not Active'
sentinel_snps$activeCM[sentinel_snps$snp %in% active_snps_cm] <- 'Active'
sentinel_snps$activeVSMC[sentinel_snps$snp %in% active_snps_vsmc] <- 'Active'
sentinel_snps$DiffCM <- 'Not Diff'
sentinel_snps$DiffVSMC <- 'Not Diff'
sentinel_snps$DiffCM[sentinel_snps$snp %in% diff_snps] <- 'Diff'
sentinel_snps$DiffVSMC[sentinel_snps$snp %in% diff_snps_vsmc] <- 'Diff'

table(sentinel_snps$DiffCM)

#### read length block -----
length <- read.delim('../../data/design/sentinel_snps_length_block.txt',
                     sep = '\t', header = T)
#1,000,000

number_snps <- as.data.frame(table(sentinel_snps$sentinel[sentinel_snps$activeDiffCM == 'ActiveDiff']))
number_snps_all <- as.data.frame(table(sentinel_snps$sentinel))
number_snps_all <- merge(number_snps_all, length, by.x='Var1', by.y='sentinel')
number_snps_all$clust <- number_snps_all$Freq/(number_snps_all$length/1000000)

number_snps <- merge(number_snps, length, by.x='Var1', by.y='sentinel')
number_snps <- merge(number_snps, number_snps_all, by='Var1')

#### all pool, clustering of snps ------
phastcons_scores <- as.data.frame(names(loc))
phastcons_scores$score <- loc$PHASTCONS
phastcons_scores <- phastcons_scores %>% distinct()
colnames(phastcons_scores) <- c('snp','score')

phastcons_scores <- merge(phastcons_scores, sentinel_snps)

phastcons_scores_g <- phastcons_scores %>% 
  group_by(sentinel) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  as.data.frame()  

number_snps_phylop <- merge(phastcons_scores_g, number_snps_all, by.x='sentinel',by.y='Var1')

ggplot(number_snps_phylop[number_snps_phylop$length>300,], aes(x=clust, y = score)) + 
  geom_point(color='darkgrey')+
  xlab('Nº SNPs / Mb')+
  ylab('mean Phastcons score')+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=FALSE) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  theme_minimal()+
  geom_point() + 
  geom_text_repel(
    data = number_snps_phylop[number_snps_phylop$length>300,],
    aes(label = sentinel),
    size = 4)

library("ggpubr")
ggscatter(number_snps_phylop[number_snps_phylop$length>300,], x = "clust", y = "score", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Clustering", ylab = "mean(PhyloP)")

#### phastcons per group -----
head(phastcons_scores)

p <- ggboxplot(phastcons_scores, x = "DiffCM", y = "score",
               color = "DiffCM", palette = "jco")
p + stat_compare_means(aes(group = DiffCM) , label = "p.format")

library(ggridges)
library(viridis)
library(hrbrthemes)

ggplot(phastcons_scores, aes(x = score, y = sig_both, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  labs(title = '') +
  theme_minimal() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )

library(plyr)
phastcons_scores$sig_both <- 'no sig'
phastcons_scores$sig_both[phastcons_scores$DiffCM == 'Diff'] <- 'sig CM'
phastcons_scores$sig_both[phastcons_scores$DiffVSMC == 'Diff'] <- 'sig VSMC'
phastcons_scores$sig_both[phastcons_scores$DiffCM == 'Diff' & phastcons_scores$DiffVSMC == 'Diff'] <- 'sig both'


mu <- ddply(phastcons_scores, "sig_both", summarise, grp.mean=mean(score, na.rm = TRUE))
head(mu)

p<-ggplot(phastcons_scores, aes(x=score, color=sig_both)) +
  geom_density()+
  xlim(c(-0.1,1.3))+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=sig_both),
             linetype="dashed") + theme_minimal()
p


p <- ggboxplot(phastcons_scores, x = "sig_both", y = "score",
               color = "sig_both", palette = "jco")
p +   stat_compare_means(method = "anova", label.y = 1.3)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "no sig")  

#### per clustering ------
number_snps <- as.data.frame(table(sentinel_snps$sentinel[sentinel_snps$DiffVSMC == 'Diff']))

number_snps <- merge(number_snps, length, by.x='Var1', by.y='sentinel')
number_snps <- merge(number_snps, number_snps_all, by='Var1')

number_snps_phylop <- merge(phastcons_scores_g, number_snps_all, by.x='sentinel',by.y='Var1')

ggplot(number_snps_phylop[number_snps_phylop$Freq.x<200,], aes(x=Freq.x, y = score)) + 
  geom_point(color='darkgrey')+
  xlab('Nº SNPs / Mb')+
  ylab('mean Phastcons score')+
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=FALSE) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  theme_minimal()+
  geom_point() + 
  geom_text_repel(
    data = number_snps_phylop[number_snps_phylop$Freq.x<200,],
    aes(label = sentinel),
    size = 4)

library("ggpubr")
ggscatter(number_snps_phylop[number_snps_phylop$Freq<200,], x = "Freq", y = "score", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Clustering", ylab = "mean(PhastCons)")


### number of ultraconserved elements ------
ultraconserved <- phastcons_scores[!is.na(phastcons_scores$score) & phastcons_scores$score > 0.8,]
ggplot(ultraconserved, aes(x=sig_both, fill=as.factor(sig_both) )) + 
  geom_bar( ) +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal()+ xlab('')+
  theme(legend.position="none") 

write.table(ultraconserved, '../../analysis/04_PhyloP/ultraconserved_snps_phastcons.txt',
            quote = F, sep = '\t', row.names = F, col.names = T)


### ultraconserved elements ----
ultra_elements <- read.table('../../data/design/overlap_ultraconserved.bed',
                             header=F, sep = '\t')
head(ultra_elements)
head(sentinel_snps)

sentinel_snps$ultraconserved <- 'No'
sentinel_snps$ultraconserved[sentinel_snps$change %in% ultra_elements$V4] <- 'YES'
table(sentinel_snps$ultraconserved)
table(sentinel_snps$ultraconserved, sentinel_snps$DiffVSMC)

chisq <- fisher.test(table(sentinel_snps$ultraconserved, sentinel_snps$activeDiffCM))
chisq


write.table(sentinel_snps[sentinel_snps$ultraconserved == 'YES' & sentinel_snps$DiffCM=='Diff',], '../../analysis/04_PhyloP/Reg_CM_ultraconserved_elemnts_overlap_snps_phastcons.txt',
            quote = F, sep = '\t', row.names = F, col.names = T)

