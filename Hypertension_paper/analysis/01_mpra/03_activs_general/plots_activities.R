library(ggplot2)
library(reshape)
library(dplyr)
library(ggrepel)
library(ggpubr)
library(tidyr)

vals_significance_vsmc <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/VSMC_vals.significance.4.txt', sep = '\t', header = T)
head(vals_significance_vsmc)
vals_significance_vsmc <- vals_significance_vsmc %>%
  separate(name, c("pos", "snp_info","info"), "__")

vals_significance <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/CM_vals.significance.5.txt', sep = '\t', header = T)
head(vals_significance)
vals_significance <- vals_significance %>%
  separate(name, c("pos", "snp_info","info"), "__")

mital <- read.table("~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/name_Mital_SNPS.txt", header = T)

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

df_VSMC = merge(VSMC_counts,index[,c("barcode", "element", "tile_type", "tile_id","snp","name")], by="barcode")
df_CM = merge(CM_counts,index[,c("barcode", "element", "tile_type", "tile_id","snp","name")], by="barcode")


#### Plot distribution of random activities ####
control_random <- df_VSMC[grep('rs62012628', df_VSMC$name),]
control_random$median_rep <- rowMeans(control_random[,c("VSMC_rep1", "VSMC_rep2", "VSMC_rep3", "VSMC_rep6")], na.rm = T)

control_random <- df_CM[grep('RANDOM', df_CM$tile_type),]
control_random$median_rep <- rowMeans(control_random[,c("CM_rep1", "CM_rep2", "CM_rep3", "CM_rep4","CM_rep5")], na.rm = T)
control_random$ratio <- control_random$median_rep/control_random$dna_1

ggplot(control_random, aes(x = name, y = ratio)) +
  #geom_jitter(aes(color=tile_type),alpha=0.4)+
  geom_point(aes(color=snp),pch = 21, alpha=0.3, position = position_jitterdodge())+
  geom_violin(outlier.shape = NA, alpha=0.5)+
  ylab("Ratio RNA/DNA") +
  xlab("") +
  theme(legend.position="top", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  ylim(0,10)
  #scale_color_brewer(palette="Set1") +

#### activities random ####
ggplot(positive_controls_CM, aes(x = activity, y = ..density..)) +
  geom_density(adjust=1.5, alpha=.5, fill='light blue', color= 'light blue') +
  geom_density(data= control_random,aes(x = ratio, y = -..density..), fill= "grey", color = 'grey') +
  ylab("Density") +
  xlab("Ratio RNA/DNA") +
  xlim(c(0,2))+
  theme(legend.position="right", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(color="", fill = "")+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  scale_color_brewer(palette="Set3") +
  scale_fill_brewer(palette = "Set3")

#### density plot of activities ####
head(vals_significance)

controls <- c('rs6445040','rs9270898','rs60002611','rs117104239','rs4360494','rs62012628')

vals_significance_vsmc$group <- NA
vals_significance_vsmc[vals_significance_vsmc$snp_info %in% controls,'group'] <- 'Control'
vals_significance_vsmc[vals_significance_vsmc$snp_info%in% controls,'group'] <- 'Control'
vals_significance_vsmc[vals_significance_vsmc$type %in% c('WILDTYPE_SNP_INDIV','WILDTYPE_BUT_HAS_SNP') & vals_significance_vsmc$sig == 'sig', 'group'] <- 'Active element'
vals_significance_vsmc[vals_significance_vsmc$type %in% c('WILDTYPE_SNP_INDIV','WILDTYPE_BUT_HAS_SNP') & vals_significance_vsmc$sig == 'not sig', 'group'] <- 'Not Active element'
vals_significance_vsmc[vals_significance_vsmc$type == 'negative control','group'] <- 'Random Sequence'

vals_significance$group <- NA
vals_significance[vals_significance$snp_info %in% controls,'group'] <- 'Control'
vals_significance[vals_significance$snp_info%in% controls,'group'] <- 'Control'
vals_significance[vals_significance$type %in% c('WILDTYPE_SNP_INDIV','WILDTYPE_BUT_HAS_SNP') & vals_significance$sig == 'sig', 'group'] <- 'Active element'
vals_significance[vals_significance$type %in% c('WILDTYPE_SNP_INDIV','WILDTYPE_BUT_HAS_SNP') & vals_significance$sig == 'not sig', 'group'] <- 'Not Active element'
vals_significance[vals_significance$type == 'negative control','group'] <- 'Random Sequence'

### mital seqs ####
mital_names <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/name_Mital_SNPS.txt', sep='\t', header=T)
mital_details <- index$parse_details[index$name %in% mital_names$name]
### filter mital seqs --> other project
vals_significance_vsmc <- vals_significance_vsmc[!vals_significance_vsmc$parse_details %in% mital_details,]
vals_significance <- vals_significance[!vals_significance$parse_details %in% mital_details,]

ggplot(vals_significance_vsmc[!vals_significance_vsmc$group == 'Random Sequence',], aes(x = VSMC_log, y = ..density..,group = group,fill=group, color = group)) +
  geom_density(adjust=1.5, alpha=.5) +
  geom_density(data= vals_significance_vsmc[vals_significance_vsmc$group == 'Random Sequence',],aes(x = VSMC_log, y = -..density..), fill= "grey", color = 'grey') +
  ylab("Density") +
  xlab("MPRA Activity") +
  xlim(c(-2,2))+
  ylim(c(-4,6))+
  theme(legend.position="right", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(color="", fill = "")+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  scale_color_brewer(palette="Set2") +
  scale_fill_brewer(palette = "Set2")

ggplot(vals_significance[!vals_significance$group == 'Random Sequence' ,], aes(x = CM_log, y = ..density..,group = group,fill=group, color = group)) +
  geom_density(adjust=1.5, alpha=.5) +
  geom_density(data= vals_significance[vals_significance$group == 'Random Sequence',],aes(x = CM_log, y = -..density..), fill= "grey", color = 'grey') +
  ylab("Density") +
  xlab("MPRA Activity") +
  xlim(c(-2,2))+
  ylim(c(-4,6))+
  theme(legend.position="right", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  labs(color="", fill = "")+
  theme(panel.background = element_rect(fill = 'white', colour = 'white'))+
  scale_color_brewer(palette="Set2") +
  scale_fill_brewer(palette = "Set2")

### how active sequences look like? #####
active_cm <- vals_significance[vals_significance$group == 'Random Sequence',]
ggplot(active_cm, aes(x = tile_type, y = CM)) +
  #geom_violin(fill='#BADACF')+
  geom_violin(fill='dark grey')+
  geom_boxplot(width=0.25)+
  ylab("MPRA Activity") +
  xlab("") +
  theme(legend.position="top", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  #labs(color="", fill = "")+
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  ylim(c(0,0.6))

ggplot(active_cm[active_cm$element != 'NA',], aes(x = group, y = VSMC_log)) +
  geom_violin(fill='#BADACF')+
  #geom_violin(fill='dark grey')+
  ylab("MPRA Activity") +
  xlab("") +
  theme(legend.position="top", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  #labs(color="", fill = "")+
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  ylim(c(-0.5,3.5))

### How many active sentinel SNPs ####
sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)

sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

dist_activies <- as.data.frame(table(vals_significance$group))
dist_activies <- dist_activies[dist_activies$Var1 %in% c('Active element','Not Active element'),]

# Compute the position of labels
dist_activies <- dist_activies %>%
  arrange(desc(Var1)) %>%
  mutate(prop = Freq / sum(dist_activies$Freq) *100) %>%
  mutate(ypos = cumsum(prop)- 0.6*prop )

ggplot(dist_activies, aes(x="", y=prop, fill=Var1)) +
  geom_bar(width=1, color="white", stat = 'identity') +
  coord_polar("y")+
  geom_text(aes(y = ypos+1, label = paste0(round(prop),'%')), color = "white", size=5) +
  ylab("")+
  xlab("")+
  ggtitle("")+
  theme(legend.position="bottom") +
  theme_void() +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_brewer(palette="Set2")

active_snps <- unique(vals_significance[vals_significance$group == 'Active element','snp_info'])

sentinel_snps$active[sentinel_snps$snp %in% active_snps] <- 'Active'
sentinel_snps$active[!sentinel_snps$snp %in% active_snps] <- 'Not Active'

sentinel <- unique(sentinel_snps$sentinel)

sentinel_snps_sent <- sentinel_snps[sentinel_snps$snp %in% sentinel,]

grouped <- sentinel_snps[,c('sentinel','Chr','Trait','active')] %>% distinct()

active_sentinel <- unique(grouped$sentinel[grouped$active == 'Active'])

sentinel_snps_sent$active[sentinel_snps_sent$sentinel %in% active_sentinel] <- 'Active'
dist_snps <- as.data.frame(table(sentinel_snps_sent$active))

# Compute the position of labels
dist_snps <- dist_snps %>%
  arrange(desc(Var1)) %>%
  mutate(prop = Freq / sum(dist_snps$Freq) *100) %>%
  mutate(ypos = cumsum(prop)- 0.6*prop )

ggplot(dist_snps, aes(x="", y=prop, fill=Var1)) +
  geom_bar(width=1, color="white", stat = 'identity') +
  coord_polar("y")+
  geom_text(aes(y = ypos+1, label = paste0(round(prop),'%')), color = "white", size=5) +
  ylab("")+
  xlab("")+
  ggtitle("Sentinel SNPs with at least one LD active SNP")+
  theme(legend.position="bottom") +
  theme_void() +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_fill_brewer(palette="Set2")

### location of active tiles ####
vep_filtered <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/vep.hs38.canonical.oneConsPerVariant.txt', sep = '\t')
vep_filtered <- vep_filtered[,c('X.Uploaded_variation','Location','Consequence','IMPACT','SYMBOL','Gene')]

vep_filtered <- vep_filtered %>% distinct()
vep_filtered$effect <- gsub(',.*','',vep_filtered$Consequence)

vep_filtered$set <- 'Original'
vep_filtered <- vep_filtered[,c('X.Uploaded_variation','Location','Consequence','IMPACT','SYMBOL','Gene','set')]

vep_filtered$set[vep_filtered$X.Uploaded_variation %in% active_snps] <- 'Active'

### get unique results ####
vep_filtered$effect <- gsub(',.*','',vep_filtered$Consequence)

get_col <- function(fdr){
  if(abs(fdr) < -log10(0.05)){
    col <- "light grey"
  }else{
    col <- col_numeric(rev(brewer.pal(11,"RdBu")),
                       domain = seq(-3,
                                    3,length.out = 11))(fdr)
  }
}

table_vep <- as.data.frame(table(vep_filtered$effect, vep_filtered$set))

# enrichment of locations in active sequences
my_fisher <- function(type){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  type_counts_original <- table_vep[table_vep$Var1 == type & table_vep$Var2 == 'Original',]$Freq
  other_type_counts_original <- sum(table_vep[table_vep$Var1 != type & table_vep$Var2 == 'Original',]$Freq)
  type_counts_random <- table_vep[table_vep$Var1 == type & table_vep$Var2 == 'Active',]$Freq
  other_type_counts_random <- sum(table_vep[table_vep$Var1 != type & table_vep$Var2 == 'Active',]$Freq)

  m <- matrix(c(type_counts_random, type_counts_original, other_type_counts_random, other_type_counts_original), 2,2, byrow = T)
  print(m)
  m[is.na(m)] <- 0
  rownames(m) <- c("Intron", "Other")
  colnames(m) <- c("Active","NonActive")
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

vep_filtered_active$effect <- gsub(',.*','',vep_filtered_active$Consequence)
table_active <- as.data.frame(table(vep_filtered$effect[vep_filtered$set == 'Active']))

dist_activies <- table_active %>%
  arrange(desc(Var1)) %>%
  mutate(prop = Freq / sum(table_active$Freq) *100) %>%
  mutate(ypos = cumsum(prop)- 0.6*prop )

ggplot(dist_activies, aes(x=Var1, y=prop, fill=Var1)) +
  geom_bar(width=1, color="white", stat = 'identity') +
  #coord_polar("y")+
  geom_text(aes(y = prop+5, label = Freq), color = "black", size=5) +
  ylab("%")+
  xlab("")+
  ggtitle("")+
  theme(legend.position="right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  #theme_void() +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  scale_fill_brewer(palette="Set3")

### active tiles across the genome ####

snps <- vals_significance[vals_significance$tile_type %in% c('WILDTYPE_BUT_HAS_SNP','WILDTYPE_SNP_INDIV'),c('pos','snp_info','chrom')] %>% distinct()

my_fisher <- function(type){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  type_counts_active <- nrow(snps[snps$chrom == type & snps$snp_info %in% active_snps,])
  other_type_counts_active <- nrow(snps[snps$chrom != type & snps$snp_info %in% active_snps,])
  type_counts_nonactive <- nrow(snps[snps$chrom == type & !snps$snp_info %in% active_snps,])
  other_type_counts_nonactive <- nrow(snps[snps$chrom != type & !snps$snp_info %in% active_snps,])

  m <- matrix(c(type_counts_active, type_counts_nonactive, other_type_counts_active, other_type_counts_nonactive), 2,2, byrow = T)
  print(m)
  m[is.na(m)] <- 0
  rownames(m) <- c("chr", "Other")
  colnames(m) <- c("Active","Non-active")
  f <- fisher.test(m)
  print(f)
  return(list("f" = f,
              "overlap" = x11))
}
# Two-tailed Fisher test
chrom <- unique(snps$chrom)
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

snps_filt <- snps[snps$snp_info %in% active_snps,]
table_active <- as.data.frame(table(snps_filt$chrom))

dist_activies <- table_active %>%
  arrange(desc(Var1)) %>%
  mutate(prop = Freq / sum(table_active$Freq) *100) %>%
  mutate(ypos = cumsum(prop)- 0.6*prop )

ggplot(dist_activies, aes(x=Var1, y=prop, fill=Var1)) +
  geom_bar(width=1, color="white", stat = 'identity') +
  #coord_polar("y")+
  geom_text(aes(y = prop+5, label = Freq), color = "black", size=5) +
  ylab("%")+
  xlab("")+
  ggtitle("")+
  theme(legend.position="right") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')+
  #theme_void() +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))
  #scale_fill_brewer(palette="Set3")




