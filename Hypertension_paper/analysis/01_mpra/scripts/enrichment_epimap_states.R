library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(reshape)

#### read SNP info ####
sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]
oligos <-  sentinel_snps %>% separate(coord, c("chr", "start","end"), ":")
oligos <- oligos[,c('chr','start','end','snp')]

#### read active tiles ####
### get diff and active snps ####
vals_significance_cm <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/CM_vals.significance.5.txt', sep = '\t', header = T)
head(vals_significance_cm)
vals_significance_cm <- vals_significance_cm %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance_cm)

vals_significance_vsmc <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/VSMC_vals.significance.4.txt', sep = '\t', header = T)
head(vals_significance_vsmc)
vals_significance_vsmc <- vals_significance_vsmc %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance_vsmc)

CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)

## index

### merge results ####
mital_names <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/name_Mital_SNPS.txt', sep='\t', header=T)
mital_details <- mital_names %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(mital_details)

sentinel_snps$ActiveCM <- 'Not Active'
sentinel_snps$ActiveVSMC <- 'Not Active'

active_snps_cm <- unique(vals_significance_cm$snp_info[vals_significance_cm$CM_padj<0.05 & !vals_significance_cm$parse_details %in% mital_details])
active_snps_vsmc <- unique(vals_significance_vsmc$snp_info[vals_significance_vsmc$VSMC_padj<0.05 & !vals_significance_vsmc$parse_details %in% mital_details])

sentinel_snps$ActiveCM[sentinel_snps$snp %in% active_snps_cm] <- 'Active'
sentinel_snps$ActiveVSMC[sentinel_snps$snp %in% active_snps_vsmc] <- 'Active'

diff_snps <- unique(CM_results$snp[CM_results$fdr_comp<=0.05 & CM_results$is_ctrl != 'control'])
diff_snps_vsmc <- unique(VSMC_results$snp[VSMC_results$fdr_comp<=0.05 & VSMC_results$is_ctrl != 'control'])

sentinel_snps$DiffCM <- 'Not Diff'
sentinel_snps$DiffVSMC <- 'Not Diff'
sentinel_snps$DiffCM[sentinel_snps$snp %in% diff_snps] <- 'Diff'
sentinel_snps$DiffVSMC[sentinel_snps$snp %in% diff_snps_vsmc] <- 'Diff'

### read DHS overlap ####
BrainVascularSmooth_dhs <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap/BrainVascularSmooth_SNPS_overlap.bed',
                                      header = F, sep = '\t', stringsAsFactors = FALSE)

Cardiac_muscle_derived_dhs <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap/Cardiac_muscle_derived_SNPS_overlap.bed',
                                         header = F, sep = '\t',stringsAsFactors = FALSE)

CardiacMyocyte_dhs <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap/CardiacMyocyte_SNPS_overlap.bed',
                                 header = F, sep = '\t',stringsAsFactors = FALSE)

CoronaryArtery1_dhs <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap/CoronaryArtery1_SNPS_overlap.bed',
                                  header = F, sep = '\t',stringsAsFactors = FALSE)

CoronaryArtery2_dhs <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap/CoronaryArtery2_SNPS_overlap.bed',
                                  header = F, sep = '\t',stringsAsFactors = FALSE)

SmoothMuscleDeriv_dhs <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap/SmoothMuscleDeriv_SNPS_overlap.bed',
                                    header = F, sep = '\t',stringsAsFactors = FALSE)


#### enrichment -- chi-squared first ####

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
my_fisher <- function(type, df){
  #             DS      Not DS
  # type
  # other_types
  print(type)
  df$New <- df$V4
  df$New[df$New %in% c('EnhG1','EnhG2','EnhA1','EnhA2','EnhWk','EnhBiv')] <- 'Enh'
  df$New[df$New %in% c('TssA','TssFlnk','TssFlnkU','TssFlnkD','TssBiv')] <- 'Prom'
  df$New[df$New %in% c('ReprPC','ReprPCWk','Het')] <- 'Repressed'
  
  type_df <- df[df$New == type,]
  other_type <- df[df$New != type,]
  type_diff <- length(unique(type_df$V13[type_df$V13 %in% sentinel_snps$change[sentinel_snps$DiffCM == 'Diff']]))
  type_notdiff <-length(unique(type_df$V13[type_df$V13 %in% sentinel_snps$change[sentinel_snps$DiffCM != 'Diff']]))
  other_type_diff <- length(unique(other_type$V13[other_type$V13 %in% sentinel_snps$change[sentinel_snps$DiffCM == 'Diff']]))
  other_type_notdiff <- length(unique(other_type$V13[other_type$V13 %in% sentinel_snps$change[sentinel_snps$DiffCM != 'Diff']]))
  
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, other_type_diff, other_type_notdiff), 2,2, byrow = T)
  print(m)
  
  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c(type, "Other")
  colnames(m) <- c("Diff","Not Diff")
  print(m)
  f <- fisher.test(m)
  print(f)
  return(list("f" = f,
              "overlap" = type_df$V13))
}

### try all tissues
families <- c('Enh','Prom','Repressed')
tissues <- list(CardiacMyocyte_dhs, Cardiac_muscle_derived_dhs, BrainVascularSmooth_dhs,
                CoronaryArtery1_dhs,CoronaryArtery2_dhs,SmoothMuscleDeriv_dhs)
names(tissues) <- c('CardiacMyocyte','Cardiac_muscle_derived','BrainVascularSmooth',
                    'CoronaryArtery1','CoronaryArtery2','SmoothMuscleDeriv')
fisher_results <- lapply(tissues, function(tis) lapply(families, function(chr) my_fisher(chr, tis)))
names(fisher_results) <- names(tissues)


#### Bubble plot ####
for (name in names(tissues)) {
  names(fisher_results[[name]]) <- families
}

odds_ratio <- lapply(names(tissues), function(tis) lapply(families, function(chr) (fisher_results[[tis]][[chr]][['f']]$estimate)))
names(odds_ratio) <- names(tissues)
for (name in names(tissues)) {
  names(odds_ratio[[name]]) <- families
}

odds_ratio_df <- do.call('rbind.data.frame',odds_ratio)
odds_ratio_df[is.na(odds_ratio_df)] <- 0
odds_ratio_df$tissue <- rownames(odds_ratio_df)

pval <- lapply(names(tissues), function(tis) p.adjust(lapply(families, function(chr)  fisher_results[[tis]][[chr]][['f']]$p.value), method = "BH"))
#pval <- lapply(names(tissues), function(tis) (lapply(families, function(chr)  fisher_results[[tis]][[chr]][['f']]$p.value)))

names(pval) <- names(tissues)
for (name in names(tissues)) {
  names(pval[[name]]) <- families
}

pval_df <- do.call('rbind.data.frame',pval)
pval_df[is.na(pval_df)] <- 1
rownames(pval_df) <- names(tissues)
colnames(pval_df) <- families
pval_df$tissue <- rownames(pval_df)

pv <- melt(pval_df, id.vars = 'tissue')
colnames(pv) <- c('tissue','type','pval')
odds <- melt(odds_ratio_df, id.vars = 'tissue')
colnames(odds) <- c('tissue','type','oddsRatio')

all_Res <- merge(odds,pv)

ggplot(all_Res, aes(x=type, y=tissue, size=oddsRatio, fill=-log10(pval)*sign(log2(oddsRatio)))) +
  geom_point(alpha=0.8, shape=21, color="black") +
  scale_size(range = c(0.5, 12), name="OddsRatio") +
  scale_fill_gradient2(low="dark blue",mid="white", high="red", name='-log10(FDR)')+
  #theme_ipsum() +
  theme(legend.position="right") +
  ylab("") +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        panel.grid = element_line(colour = 'light grey'))

# cm 662b8fff
# vsmc f05928ff
all_Res$sig <- 'not sig'
all_Res$sig[all_Res$pval<0.05] <- 'sig'
ggplot(all_Res, aes(x=type, y=tissue, size=oddsRatio, fill=-log10(pval)*sign(log2(oddsRatio)))) +
  geom_point(aes(color=sig),alpha=0.8, shape=21) +
  scale_size(range = c(0.5, 12), name="OddsRatio", limits = c(0,3)) +
  scale_fill_gradient2(low="dark blue",mid="white", high="#f05928ff", name='-log10(pval)')+
  scale_color_manual(values = c("not sig" = "white",
                                "sig" = "black"))+
  #theme_ipsum() +
  theme(legend.position="right") +
  ylab("") +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        panel.grid = element_line(colour = 'light grey'))

#### other files ####
files_enh <- list.files("~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap/",
                        pattern = '_enh_SNPS_overlap.bed', full.names = TRUE)
files_prom <- list.files("~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap/",
                         pattern = '_prom_SNPS_overlap.bed', full.names = TRUE)

tissues <- c('CARD.MUSCL_','SM.MUSCLE_','CORONARY.ATY_','CORONARY.ATY2_','CARDIAC.MYOCYT_','BRN.VASC.SMTH.MUSC_')
tissues <- c('CARDIAC.MYOCYT_','CARD.MUSCL_','BRN.VASC.SMTH.MUSC_','CORONARY.ATY_','CORONARY.ATY2_','SM.MUSCLE_')

my_fisher <- function(tissue,type){
  #             DS      Not DS
  # type
  # other_types
  if (type =='enh') {
    file <- files_enh[grep(tissue, files_enh)]
    df <- read.table(file, header = F, sep = '\t')
  } else {
    file <- files_prom[grep(tissue, files_prom)]
    df <- read.table(file, header = F, sep = '\t')
  }
  
  print(nrow(df))
  
  type_df <- nrow(df)
  other_type <- 4610-type_df
  type_diff <- nrow(df[df$V8 %in% sentinel_snps$change[sentinel_snps$DiffCM == 'Diff'],])
  type_notdiff <- nrow(df[df$V8 %in% sentinel_snps$change[sentinel_snps$DiffCM != 'Diff'],])
  #other_type_diff <- nrow(other_type[other_type$V13 %in% sentinel_snps$change[sentinel_snps$ActiveCM == 'Active'] & other_type$V13 %in% sentinel_snps$change[sentinel_snps$DiffCM == 'Diff'],])
  #other_type_notdiff <- nrow(other_type[other_type$V13 %in% sentinel_snps$change[sentinel_snps$ActiveCM != 'Active'] | other_type$V13 %in% sentinel_snps$change[sentinel_snps$DiffCM != 'Diff'],])
  diff <- nrow(sentinel_snps[sentinel_snps$DiffCM=='Diff',])
  
  ### test significance
  f <- phyper(type_diff-1, type_df, nrow(sentinel_snps)-type_df, diff, lower.tail = F, log.p = F)
  print(f)
  ### odds ratio
  ### test significance
  m <- matrix(c(type_diff, type_notdiff, diff-type_diff , 4610-diff-type_notdiff), 2,2, byrow = T)
  #print(m)
  
  m[is.na(m)] <- 0
  #m <- m[c(type,paste0('No ',type)),]
  rownames(m) <- c("ENH", "Other")
  colnames(m) <- c("Diff","Not Diff")
  print(m)
  f_o <- fisher.test(m)
  print(f_o$estimate)
  print(f_o)
  
  return(list("f" = f,
              "overlap" = type_diff,
              "odds_ratio" = f_o$estimate,
              "total" = df$V8,
              "fisher" = f_o))
}

##do fisher
fisher_results <- lapply(tissues, function(tis) lapply(c('enh','prom'), function(chr) my_fisher(tis, chr)))
names(fisher_results) <- tissues

for (name in tissues) {
  names(fisher_results[[name]]) <- c('enh','prom')
}

odds_ratio <- list()
adj.P.Val <- list()
# odds ratio
odds_ratio <- lapply(tissues, function(chr) sapply(c('enh','prom'), function(type) fisher_results[[chr]][[type]][['odds_ratio']] ))
# adj.P.Val
adj.P.Val <- lapply(tissues, function(chr) sapply(c('enh','prom'), function(type) fisher_results[[chr]][[type]][['f']] ))
### elements into each category
enh_prom_elements <- lapply(tissues, function(chr) p.adjust(sapply(c('enh','prom'), function(type) fisher_results[[chr]][[type]][['total']], method = "BH")))

names(enh_prom_elements) <- tissues
saveRDS(enh_prom_elements, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/03_EpiMapOverlap/enh_prom_elements_overlap_per_tissue.rds')

names(odds_ratio) <- tissues
names(adj.P.Val) <- tissues
odds_ratio <- do.call('rbind.data.frame', odds_ratio)
colnames(odds_ratio) <- c('Enh', 'Prom')
rownames(odds_ratio) <- tissues

odds_ratio$fam <- gsub('_*','',rownames(odds_ratio))
odds_ratio <- melt(odds_ratio, id.vars = 'fam')
colnames(odds_ratio) <- c('fam','type','oddsratio')

adj.P.Val <- do.call('rbind.data.frame', adj.P.Val)
colnames(adj.P.Val) <- c('Enh', 'Prom')
rownames(adj.P.Val) <- tissues
adj.P.Val$fam <- gsub('_*','',rownames(adj.P.Val))
adj.P.Val <- melt(adj.P.Val, id.vars = 'fam')
colnames(adj.P.Val) <- c('fam','type','pval')

all <- merge(odds_ratio, adj.P.Val)
all$thres <- FALSE
all$thres[all$pval < 0.05] <- TRUE

## bubble
# cm 662b8fff
# vsmc f05928ff
ggplot(all, aes(x=type, y=fam, size=oddsratio, fill=-log10(pval)*sign(log2(oddsratio)))) +
  geom_point(alpha=0.8, shape=21, color="black") +
  scale_size(range = c(0.5, 12), name="OddsRatio", limits = c(0,2)) +
  scale_fill_gradient2(low="dark blue",mid="white", high="#662b8fff", name='-log10(FDR)')+
  #theme_ipsum() +
  theme(legend.position="right") +
  ylab("") +
  xlab("") +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        panel.grid = element_line(colour = 'light grey'))

