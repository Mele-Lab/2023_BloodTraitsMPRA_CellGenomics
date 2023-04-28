library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(reshape)

#### read SNP info ####
sentinel_snps = read.table('../../../data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]
oligos <-  sentinel_snps %>% separate(coord, c("chr", "start","end"), ":")
oligos <- oligos[,c('chr','start','end','snp')]

#### read active tiles ####
### get diff and active snps ####
vals_significance_cm <- read.table('../../../data/02_activs/CM_vals.significance.5.txt', sep = '\t', header = T)
head(vals_significance_cm)
vals_significance_cm <- vals_significance_cm %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance_cm)

vals_significance_vsmc <- read.table('../../../data/02_activs/VSMC_vals.significance.4.txt', sep = '\t', header = T)
head(vals_significance_vsmc)
vals_significance_vsmc <- vals_significance_vsmc %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance_vsmc)

CM_results <- read.table('../../../analysis/01_mpra/06_process_mpranalyze_comp/CM_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)
VSMC_results <- read.table('../../../analysis/01_mpra/06_process_mpranalyze_comp/VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt', sep='\t', header=T)

## index

### merge results ####
mital_names <- read.table('../../../data/design/name_Mital_SNPS.txt', sep='\t', header=T)
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
BrainVascularSmooth_dhs <- read.table('../../../analysis/03_EpiMapOverlap/Overlap_beds/BrainVascularSmooth_SNPS_overlap.bed',
                     header = F, sep = '\t', stringsAsFactors = FALSE)

Cardiac_muscle_derived_dhs <- read.table('../../../analysis/03_EpiMapOverlap/Overlap_beds/Cardiac_muscle_derived_SNPS_overlap.bed',
                                      header = F, sep = '\t',stringsAsFactors = FALSE)

CardiacMyocyte_dhs <- read.table('../../../analysis/03_EpiMapOverlap/Overlap_beds/CardiacMyocyte_SNPS_overlap.bed',
                                      header = F, sep = '\t',stringsAsFactors = FALSE)

CoronaryArtery1_dhs <- read.table('../../../analysis/03_EpiMapOverlap/Overlap_beds/CoronaryArtery1_SNPS_overlap.bed',
                                      header = F, sep = '\t',stringsAsFactors = FALSE)

CoronaryArtery2_dhs <- read.table('../../../analysis/03_EpiMapOverlap/Overlap_beds/CoronaryArtery2_SNPS_overlap.bed',
                                      header = F, sep = '\t',stringsAsFactors = FALSE)

SmoothMuscleDeriv_dhs <- read.table('../../../analysis/03_EpiMapOverlap/Overlap_beds/SmoothMuscleDeriv_SNPS_overlap.bed',
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
  type_diff <- nrow(type_df[type_df$V13 %in% sentinel_snps$change[sentinel_snps$DiffVSMC == 'Diff'],])
  type_notdiff <- nrow(type_df[type_df$V13 %in% sentinel_snps$change[sentinel_snps$DiffVSMC != 'Diff'],])
  other_type_diff <- nrow(other_type[other_type$V13 %in% sentinel_snps$change[sentinel_snps$DiffVSMC == 'Diff'],])
  other_type_notdiff <- nrow(other_type[other_type$V13 %in% sentinel_snps$change[sentinel_snps$DiffVSMC != 'Diff'],])

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

#adj.P.Val <-sapply(families, function(chr)  p.adjust(fisher_results[[chr]][['f']]), method = "BH"))
pval <- lapply(names(tissues), function(tis) p.adjust(lapply(families, function(chr)  fisher_results[[tis]][[chr]][['f']]$p.value), method = "BH"))
pval <- lapply(names(tissues), function(tis) (lapply(families, function(chr)  fisher_results[[tis]][[chr]][['f']]$p.value)))

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

### proportion of our SNPs in each category ####
get_prop <- function(df){
  df$New <- df$V4
  df$New[df$New %in% c('EnhG1','EnhG2','EnhA1','EnhA2','EnhWk','EnhBiv')] <- 'Enh'
  df$New[df$New %in% c('TssA','TssFlnk','TssFlnkU','TssFlnkD','TssBiv')] <- 'Prom'

  proportions <- as.data.frame((table(df$New)/4610)*100)
  print(proportions)

  return(proportions)
  }

length(Cardiac_muscle_derived_dhs$V13[Cardiac_muscle_derived_dhs$V13 %in% sentinel_snps$change])

#families <- c('Enh','Prom','Tx','TxWk','ZNF/Rpts','Het','ReprPC','ReprPCWk','Quies')
tissues <- list(CardiacMyocyte_dhs, Cardiac_muscle_derived_dhs, BrainVascularSmooth_dhs,
                CoronaryArtery1_dhs,CoronaryArtery2_dhs,SmoothMuscleDeriv_dhs)
names(tissues) <- c('CardiacMyocyte','Cardiac_muscle_derived','BrainVascularSmooth',
                    'CoronaryArtery1','CoronaryArtery2','SmoothMuscleDeriv')
proportions <- lapply(tissues, function(tis) get_prop(tis))
names(proportions) <- names(tissues)

df_prop <- do.call('rbind.data.frame', proportions)
df_prop$tissue <- gsub('\\..*','',rownames(df_prop))

ggplot(df_prop, aes(fill=Var1, y=Freq, x=tissue)) +
  geom_bar(position="fill", stat="identity")+
  coord_flip()+
  scale_fill_brewer(palette="Set3") +
  xlab('')+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        axis.text = element_text(size = 10))


### Do enhancers and promoters at the same time ####
files_enh <- list.files("../../../analysis/03_EpiMapOverlap/Overlap_beds/",
                    pattern = '_enh_SNPS_overlap.bed', full.names = TRUE)
files_prom <- list.files("analysis/03_EpiMapOverlap/Overlap_beds/",
                         pattern = '_prom_SNPS_overlap.bed', full.names = TRUE)

tissues <- c('CARD.MUSCL_','SM.MUSCLE_','CORONARY.ATY_','CORONARY.ATY2_','CARDIAC.MYOCYT_','BRN.VASC.SMTH.MUSC_')
tissues <- c('CARDIAC.MYOCYT_','CARD.MUSCL_','BRN.VASC.SMTH.MUSC_','CORONARY.ATY_','CORONARY.ATY2_','SM.MUSCLE_')

proportion <- lapply(names(tissues), function(tis) lapply(families, function(chr) (fisher_results[[tis]][[chr]][['overlap']])))
names(proportion) <- tissues

for (name in tissues) {
  names(proportion[[name]]) <- families
}
saveRDS(proportion, '../../../analysis/03_EpiMapOverlap/proportions_enh_prom_repr.rds')

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
  type_diff <- nrow(df[df$V8 %in% sentinel_snps$change[sentinel_snps$ActiveCM == 'Active'],])
  type_notdiff <- nrow(df[df$V8 %in% sentinel_snps$change[sentinel_snps$ActiveCM != 'Active'],])
  #other_type_diff <- nrow(other_type[other_type$V13 %in% sentinel_snps$change[sentinel_snps$ActiveCM == 'Active'] & other_type$V13 %in% sentinel_snps$change[sentinel_snps$DiffCM == 'Diff'],])
  #other_type_notdiff <- nrow(other_type[other_type$V13 %in% sentinel_snps$change[sentinel_snps$ActiveCM != 'Active'] | other_type$V13 %in% sentinel_snps$change[sentinel_snps$DiffCM != 'Diff'],])
  diff <- nrow(sentinel_snps[sentinel_snps$ActiveVSMC=='Active',])

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
adj.P.Val <-  lapply(tissues, function(chr) p.adjust(sapply(c('enh','prom'), function(type)  fisher_results[[chr]][[type]][['f']]), method = "BH"))
### elements into each category
enh_prom_elements <- lapply(tissues, function(chr) sapply(c('enh','prom'), function(type) fisher_results[[chr]][[type]][['total']] ))
names(enh_prom_elements) <- tissues
saveRDS(enh_prom_elements, '../../../analysis/03_EpiMapOverlap/enh_prom_elements_overlap_per_tissue.rds')

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

ggplot(all, aes(x=oddsratio, y=fam)) +
  geom_bar(aes(fill=type, alpha = as.factor(thres), color=type),stat = 'identity', position=position_dodge()) +
  ggtitle('location enrichment') +
  ylab('')+
  scale_fill_brewer(palette="Reds")+
  scale_color_brewer(palette="Reds")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50')) +
  geom_vline(aes(xintercept=1), linetype="dotted")

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


#### DHS regions overlap 2.5kb with SNPs #####
DHS <- read.table('../../../analysis/03_EpiMapOverlap/ALL_DHS_Epimap_core_SNPS_overlap.bed', header = F, sep = '\t')

head(DHS)

### N of SNPs per DHS ####
number_snps <- as.data.frame(table(DHS$V4))
ggplot(number_snps, aes(x=Freq)) +
  geom_density(adjust = 1/3, color='darkgrey')+
  xlab('N of LD SNPs')+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))

hist(number_snps$Freq, col = "steelblue", frame = FALSE,
     breaks = 60)

summary(number_snps$Freq)

# ### N of DHS per SNP ####
# number_snps <- as.data.frame(table(DHS$V8))
# ggplot(number_snps, aes(x=Freq)) +
#   geom_density(adjust = 1/3, color='darkgrey')+
#   xlab('N of LD SNPs')+
#   theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))
#
# hist(number_snps$Freq, col = "steelblue",
#      breaks = 60)
#
# summary(number_snps$Freq)
#
# ### prop of DHS overlap per group #####
# df <- sentinel_snps
# df$DHS <- 'No DHS'
# df$DHS[df$change %in% DHS$V8] <- 'DHS'
#
# cm_dhs <- as.data.frame(table(df$DHS, df$ActiveCM))
# cm_dhs$Perc <- ""
# for (element in c('DHS','No DHS')){
#   for (elem2 in c('Active','Not Active')) {
#     cm_dhs$Perc[cm_dhs$Var1 == element & cm_dhs$Var2 == elem2] <- cm_dhs$Freq[cm_dhs$Var1 == element & cm_dhs$Var2 == elem2]/sum(cm_dhs$Freq[cm_dhs$Var2 == elem2])
#
#   }
# }
#
# ggplot(cm_dhs[cm_dhs$Var1 == 'DHS',], aes(y=as.numeric(Perc), x=Var2, fill=Var2)) +
#   geom_bar(stat = 'identity') +
#   ylim(c(0,1))+
#   ylab('Perc')+
#   theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
#         legend.position = 'none', axis.text = element_text(size = 10)) +
#   xlab('')
#
# ### test significance
# tab <- table(df$DHS, df$ActiveCM)
# test <- fisher.test(tab)
# test2 <- chisq.test(tab)


### how many are enhancers? ####
### these data have been dowloaded from the EpiMap repository
enhancers <- read.table('~/marenostrum/Data/Hypertension/EpiMap/ENH_masterlist_indices_0indexed.tsv', header = F)
map_row_id <- read.table('~/marenostrum/Data/Hypertension/EpiMap/masterlist_DHSs_733samples_WM20180608_all_coords_hg19_r25_e100_names.core.srt.tsv', header = T)
old_new_id <- read.table('~/marenostrum/Data/Hypertension/EpiMap/masterlist_DHSs_733samples_WM20180608_all_chunkIDs2indexIDs.txt', header = F)
head(enhancers)
head(map_row_id)
head(DHS)
head(old_new_id)

enhancers_old_id <- map_row_id$name[map_row_id$cls %in% enhancers$V1]
new_id_enhancers <- old_new_id$V1[old_new_id$V2 %in% enhancers_old_id]

DHS$enhancers <- 'No'
DHS$enhancers[DHS$V4 %in% new_id_enhancers] <- 'Enhancer'

length(unique(DHS$V8[DHS$enhancers != 'No']))
#[1] 3149
length(unique(DHS$V4[DHS$enhancers != 'No']))
#[1] 4803

### how many are promoters? ####
### this table has been downloaded from the EpiMap repository
promoters <- read.table('~/marenostrum/Data/Hypertension/EpiMap/PROM_masterlist_indices_0indexed.tsv', header = F)

promoters_old_id <- map_row_id$name[map_row_id$cls %in% promoters$V1]
new_id_promoters <- old_new_id$V1[old_new_id$V2 %in% promoters_old_id]

DHS$promoters <- 'No'
DHS$promoters[DHS$V4 %in% new_id_promoters] <- 'Promoter'

length(unique(DHS$V8[DHS$promoters != 'No']))
#[1] 225
length(unique(DHS$V4[DHS$promoters != 'No']))
#[1] 509

### how many are dyadic? ####
### this table has been downloaded from the EpiMap repository
dyadic <- read.table('~/marenostrum/Data/Hypertension/EpiMap/DYADIC_masterlist_indices_0indexed.tsv', header = F)

dyadic_old_id <- map_row_id$name[map_row_id$cls %in% dyadic$V1]
new_id_dyadic <- old_new_id$V1[old_new_id$V2 %in% dyadic_old_id]

DHS$dyadic <- 'No'
DHS$dyadic[DHS$V4 %in% new_id_dyadic] <- 'Dyadic'

length(unique(DHS$V8[DHS$dyadic != 'No']))
#[1] 347
length(unique(DHS$V4[DHS$dyadic != 'No']))
#[1] 313

#### write table ----
write.table(DHS, '../../../analysis/03_EpiMapOverlap/DHS_overlap_with_info.txt',
            sep = '\t', quote = F, col.names = T, row.names = F)

### make proportion plots ####
sentinel_snps$promoter <- 'No promoter'
sentinel_snps$promoter[sentinel_snps$change %in% DHS$V8[DHS$promoters != 'No']] <- 'promoter'

sentinel_snps$enhancer <- 'No enhancer'
sentinel_snps$enhancer[sentinel_snps$change %in% DHS$V8[DHS$enhancers != 'No']] <- 'enhancer'

cm_dhs <- as.data.frame(table(sentinel_snps$promoter, sentinel_snps$ActiveCM))
cm_dhs$Perc <- ""
for (element in c('promoter','No promoter')){
  for (elem2 in c('Active','Not Active')) {
    cm_dhs$Perc[cm_dhs$Var1 == element & cm_dhs$Var2 == elem2] <- cm_dhs$Freq[cm_dhs$Var1 == element & cm_dhs$Var2 == elem2]/sum(cm_dhs$Freq[cm_dhs$Var2 == elem2])

  }
}

ggplot(cm_dhs[cm_dhs$Var1 == 'promoter',], aes(y=as.numeric(Perc), x=Var2, fill=Var2)) +
  geom_bar(stat = 'identity') +
  ylim(c(0,0.2))+
  ylab('Perc')+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        legend.position = 'none', axis.text = element_text(size = 10)) +
  xlab('')

### test significance
df$promoter <- factor(df$promoter, levels = c('promoter','No promoter'))
df$enhancer <- factor(df$enhancer, levels = c('enhancer','No enhancer'))
tab <- table(df$promoter, df$ActiveCM)
test <- fisher.test(tab)
test2 <- chisq.test(tab)

### merge with DHS vocabulary component ####
## this table has been downloaded from the EpiMap repository
DHS_index <- read.table('~/marenostrum/Data/Hypertension/EpiMap/DHS_Index_and_Vocabulary_hg38_WM20190703.txt', header = F, sep = '\t')
head(DHS_index)
DHS_merged <- merge(DHS, DHS_index[,c('V4','V10','V6')], by = 'V4')
head(DHS_merged)
remove(DHS_index)

#### merge with sentinel info ####
head(sentinel_snps)
sentinel_snps$DHS <- 'No'
sentinel_snps$DHS[sentinel_snps$change %in% DHS_merged$V8] <- 'DHS_EpiMap'

sentinel_snps$Enhancer <- 'No'
sentinel_snps$Enhancer[sentinel_snps$change %in% DHS_merged$V8[DHS_merged$enhancers!='No']] <- 'Enhancer'

sentinel_snps$Promoter <- 'No'
sentinel_snps$Promoter[sentinel_snps$change %in% DHS_merged$V8[DHS_merged$promoters!='No']] <- 'Promoter'

sentinel_snps$Dyadic <- 'No'
sentinel_snps$Dyadic[sentinel_snps$change %in% DHS_merged$V8[DHS_merged$dyadic!='No']] <- 'Dyadic'

#### do plots ####
list_snps <- list(sentinel_snps$snp, sentinel_snps$snp[sentinel_snps$DHS !='No'],
                  sentinel_snps$snp[sentinel_snps$Enhancer!='No'],
                  sentinel_snps$snp[sentinel_snps$Promoter!='No'],
                  sentinel_snps$snp[sentinel_snps$Dyadic!='No'])
lt = list(All = as.vector(sentinel_snps$snp),
          DHS = as.vector(sentinel_snps$snp[sentinel_snps$DHS !='No']),
          Enhancer = as.vector(sentinel_snps$snp[sentinel_snps$Enhancer!='No']),
          Promoter = as.vector(sentinel_snps$snp[sentinel_snps$Promoter!='No']),
          Dyadic = as.vector(sentinel_snps$snp[sentinel_snps$Dyadic!='No']))
names(list_snps) <- c('All', 'DHS', 'Enhancer','Promoter','Dyadic')

library(UpSetR)
m2 = list_to_matrix(list_snps, universal_set = sentinel_snps$snp)
m2_mat <- make_comb_mat(m2)
upset(m2_mat)
UpSet(m2_mat, set_order = c("All", "DHS", "Enhancer","Promoter","Dyadic"),
      comb_order = order(comb_size(m2_mat), decreasing = T),
      right_annotation = upset_right_annotation(m2_mat, gp = gpar(fill = c("#474448","#8D818C",'#A5A299','#B4B8C5','#E9EBF8')),
                                                add_numbers = TRUE),
      top_annotation = upset_top_annotation(m2_mat, add_numbers = TRUE))


ht = draw(UpSet(m2_mat, set_order = c("All", "DHS", "Enhancer","Promoter","Dyadic"),
                comb_order = order(comb_size(m2_mat), decreasing = T),
                right_annotation = upset_right_annotation(m2_mat, gp = gpar(fill = c("#474448","#8D818C",'#A5A299','#B4B8C5','#E9EBF8')),
                                                          add_numbers = TRUE),
                top_annotation = upset_top_annotation(m2_mat, add_numbers = TRUE)))
print(ht)
od = column_order(ht)
cs = comb_size(m2_mat)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
            default.units = "native", just = "bottom", gp = gpar(fontsize = 8))
})


head(sentinel_snps)
unique_dhs <- DHS_merged[,c('V4','V10','enhancers','promoters','dyadic')] %>% distinct()
tissues <- as.data.frame(table(unique_dhs$V10)[table(unique_dhs$V10)>0])

library(RColorBrewer)
# Classic palette BuPu, with 4 colors
coul <- brewer.pal(1, "PuBu")
# Add more colors to this palette :
coul <- colorRampPalette(coul)(17)

ggplot(tissues, aes(y=Freq, x=Var1, fill=Var1)) +
  geom_bar(stat = 'identity')+
  xlab('DHS annotation')+
  coord_flip()+
  scale_fill_manual(values =coul)+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        legend.position = 'None')

### subset by snp status , active ####
results <- lapply(c('DHS','Enhancer','Promoter','Dyadic'), function(value) as.data.frame(table(sentinel_snps[[value]], sentinel_snps$ActiveCM)))
table(sentinel_snps$ActiveVSMC, sentinel_snps$DHS)
names(results) <- c('DHS','Enhancer','Promoter','Dyadic')

results_df <- do.call('rbind.data.frame', results)
results_df <- results_df[results_df$Var2 != 'Not Active',]
results_df <- results_df[results_df$Var1 != 'No',]
results_df$prop <- (results_df$Freq/340)*100

palette =  c("#EFBC9B","#EE92C2",'#9D6A89','#725D68')

results_df$Var1 <- factor(results_df$Var1, levels = c('DHS_EpiMap','Enhancer','Promoter','Dyadic'))

ggplot(results_df, aes(x=Var1, y=prop, fill=Var1))+
  geom_bar(stat = 'identity') +
  scale_x_discrete(limits = rev(levels(results_df$Var1)))+
  xlab('')+
  coord_flip()+
  scale_fill_manual(values = palette)+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        legend.position = 'None')

### subset by snp status , differential ####
results <- lapply(c('DHS','Enhancer','Promoter','Dyadic'), function(value) as.data.frame(table(sentinel_snps[[value]], sentinel_snps$DiffCM)))
table(sentinel_snps$ActiveVSMC, sentinel_snps$DHS)
names(results) <- c('DHS','Enhancer','Promoter','Dyadic')

results_df <- do.call('rbind.data.frame', results)
results_df <- results_df[results_df$Var2 != 'Not Diff',]
results_df <- results_df[results_df$Var1 != 'No',]
results_df$prop <- (results_df$Freq/1788)*100

palette =  c("#EFBC9B","#EE92C2",'#9D6A89','#725D68')

results_df$Var1 <- factor(results_df$Var1, levels = c('DHS_EpiMap','Enhancer','Promoter','Dyadic'))

ggplot(results_df, aes(x=Var1, y=prop, fill=Var1))+
  geom_bar(stat = 'identity') +
  scale_x_discrete(limits = rev(levels(results_df$Var1)))+
  xlab('')+
  coord_flip()+
  scale_fill_manual(values = palette)+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        legend.position = 'None')

### subset by snp status , diffactive ####
sentinel_snps$DiffAct <- 'Not DifAct'
sentinel_snps$DiffAct[sentinel_snps$DiffCM == 'Diff' & sentinel_snps$ActiveCM == 'Active'] <- 'DiffAct'

sentinel_snps$DiffActVSMC <- 'Not DifAct'
sentinel_snps$DiffActVSMC[sentinel_snps$DiffVSMC == 'Diff' & sentinel_snps$ActiveVSMC == 'Active'] <- 'DiffAct'

results <- lapply(c('DHS','Enhancer','Promoter','Dyadic'), function(value) as.data.frame(table(sentinel_snps[[value]], sentinel_snps$DiffAct)))
names(results) <- c('DHS','Enhancer','Promoter','Dyadic')

results_df <- do.call('rbind.data.frame', results)
results_df <- results_df[results_df$Var2 != 'Not DifAct',]
results_df <- results_df[results_df$Var1 != 'No',]
results_df$prop <- (results_df$Freq/152)*100

palette =  c("#EFBC9B","#EE92C2",'#9D6A89','#725D68')

results_df$Var1 <- factor(results_df$Var1, levels = c('DHS_EpiMap','Enhancer','Promoter','Dyadic'))

ggplot(results_df, aes(x=Var1, y=prop, fill=Var1))+
  geom_bar(stat = 'identity') +
  scale_x_discrete(limits = rev(levels(results_df$Var1)))+
  xlab('')+
  coord_flip()+
  scale_fill_manual(values = palette)+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'),
        legend.position = 'None')

### Distance to the nearest enh center ####
nearest <- read.table('../../../analysis/03_EpiMapOverlap/closest_gene_SNPS_overlap.bed', header = F, sep = '\t')
head(nearest)
head(sentinel_snps)

sentinel_snps_merged <- merge(sentinel_snps, nearest[,c('V4','V8','V10')], by.x = 'change', by.y = 'V4')

sentinel_snps_merged$type <- 'other'
sentinel_snps_merged$type[sentinel_snps_merged$DiffCM == 'Diff'] <- 'Diff'
sentinel_snps_merged$type[sentinel_snps_merged$ActiveCM == 'Active'] <- 'Active'
sentinel_snps_merged$type[sentinel_snps_merged$DiffAct == 'DiffAct'] <- 'DiffAct'

sentinel_snps_merged$typeVSMC <- 'other'
sentinel_snps_merged$typeVSMC[sentinel_snps_merged$DiffVSMC == 'Diff'] <- 'Diff'
sentinel_snps_merged$typeVSMC[sentinel_snps_merged$ActiveVSMC == 'Active'] <- 'Active'
sentinel_snps_merged$typeVSMC[sentinel_snps_merged$DiffActVSMC == 'DiffAct'] <- 'DiffAct'

## what is the closes element? ####
sentinel_snps_merged$element <- 'DHS'
sentinel_snps_merged$element[sentinel_snps_merged$V8 %in% new_id_enhancers] <- 'Enhancer'
sentinel_snps_merged$element[sentinel_snps_merged$V8 %in% new_id_promoters] <- 'Promoter'
sentinel_snps_merged$element[sentinel_snps_merged$V8 %in% new_id_dyadic] <- 'Dyadic'

write.table(sentinel_snps_merged, '../../../analysis/03_EpiMapOverlap/sentinel_snps_closest_element.txt',
            sep = '\t', quote = F, col.names = T, row.names = F)

#### Do boxplots with abs(logFC) and MPRAact #####
#### MPRA activity ####

vals_significance_cm <- read.table('../../../data/02_activs/CM_vals.significance.5.txt', sep = '\t', header = T)
head(vals_significance_cm)
vals_significance_cm <- vals_significance_cm %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance_cm)

vals_significance_vsmc <- read.table('../../../data/02_activs/VSMC_vals.significance.4.txt', sep = '\t', header = T)
head(vals_significance_vsmc)
vals_significance_vsmc <- vals_significance_vsmc %>%
  separate(name, c("pos", "snp_info","info"), "__")
head(vals_significance_vsmc)

vals_significance_cm$group <- 'Other'
vals_significance_cm$group[vals_significance_cm$info %in% DHS_merged$V8] <- 'Other DHS'
vals_significance_cm$group[vals_significance_cm$info %in% DHS_merged$V8[DHS_merged$enhancers != 'No']] <- 'Enhancer Core'
vals_significance_cm$group[vals_significance_cm$info %in% DHS_merged$V8[DHS_merged$promoters != 'No']] <- 'Promoter'
vals_significance_cm$group[vals_significance_cm$info %in% DHS_merged$V8[DHS_merged$dyadic != 'No']] <- 'Dyadic'

table(vals_significance_cm$group)
vals_significance_cm <- vals_significance_cm[vals_significance_cm$info %in% sentinel_snps$change,]

vals_significance_vsmc$group <- 'Other'
vals_significance_vsmc$group[vals_significance_vsmc$info %in% DHS_merged$V8] <- 'Other DHS'
vals_significance_vsmc$group[vals_significance_vsmc$info %in% DHS_merged$V8[DHS_merged$enhancers != 'No']] <- 'Enhancer Core'
vals_significance_vsmc$group[vals_significance_vsmc$info %in% DHS_merged$V8[DHS_merged$promoters != 'No']] <- 'Promoter'
vals_significance_vsmc$group[vals_significance_vsmc$info %in% DHS_merged$V8[DHS_merged$dyadic != 'No']] <- 'Dyadic'
vals_significance_vsmc <- vals_significance_vsmc[vals_significance_vsmc$info %in% sentinel_snps$change,]
table(vals_significance_vsmc$group)

vals_significance_cm$group <- factor(vals_significance_cm$group, levels = c('Other','Other DHS','Enhancer Core','Dyadic','Promoter'))
vals_significance_vsmc$group <- factor(vals_significance_vsmc$group, levels = c('Other','Other DHS','Enhancer Core','Dyadic','Promoter'))


foo <- pairwise.wilcox.test(vals_significance_vsmc$VSMC, vals_significance_vsmc$group, p.adjust.method="none")
foo

library(ggpubr)
#c("Other", "Other DHS"), c("Other", "Enhancer Core"), c("Other", "Dyadic"),
my_comparisons <- list(
                        c('Other','Promoter'), c('Promoter','Other DHS'),c('Promoter','Enhancer Core'),
                        c('Promoter','Dyadic'))

ggplot(vals_significance_vsmc, aes(x = group, y = log2(VSMC))) +
  geom_violin(aes(fill=group), notch=F, outlier.size=0.5, outlier.alpha = 0.5) +
  ggtitle('') +
  #stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  ylab('MPRA Activity')+
  xlab('')+
  #ylim(c(-1,20))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_brewer(palette="Reds") + theme_classic() +
  stat_compare_means(comparisons = my_comparisons, label.y = c(10, 12, 14,16))


install.packages("ggridges")
library(ggridges)

ggplot(vals_significance_cm, aes(x = CM_log, y = group, fill=group)) +
  geom_density_ridges(rel_min_height = 0.005, quantile_lines = TRUE, alpha = 0.75,
                      quantiles = 2) +
  xlim(c(-3,2.5)) + theme_classic() +
  theme(legend.position = 'None')+
  scale_fill_brewer(palette="Reds") +
  ylab('')


### index
index <- read.table('../../../data/design/Hypertension__pooled.index.txt', sep = '\t', header = T)
head(index)

### read differential results ####
CM_results <- read.table('../../../analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.CM.new_back.005.txt', sep='\t', header = T)
VSMC_results <- read.table('../../../analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_with_significance.VSMC.new_back.005.txt', sep='\t', header = T)

#### merge results with index ####
index['dupe_info'] <- gsub('\\..*','',index$tile_id)

CM_results = merge(CM_results,index[index$tile_type != 'RANDOM',c("name","dupe_info")], by="dupe_info")
VSMC_results = merge(VSMC_results,index[index$tile_type != 'RANDOM',c("name","dupe_info")], by="dupe_info")

CM_results <- CM_results %>% distinct()
VSMC_results <- VSMC_results %>% distinct()

CM_results$info <- gsub('.*__','',CM_results$name)
VSMC_results$info <- gsub('.*__','',VSMC_results$name)


CM_results$group <- 'Other'
CM_results$group[CM_results$info %in% DHS_merged$V8] <- 'Other DHS'
CM_results$group[CM_results$info %in% DHS_merged$V8[DHS_merged$enhancers != 'No']] <- 'Enhancer Core'
CM_results$group[CM_results$info %in% DHS_merged$V8[DHS_merged$promoters != 'No']] <- 'Promoter'
CM_results$group[CM_results$info %in% DHS_merged$V8[DHS_merged$dyadic != 'No']] <- 'Dyadic'

table(CM_results$group)
CM_results <- CM_results[CM_results$info %in% sentinel_snps$change,]

VSMC_results$group <- 'Other'
VSMC_results$group[VSMC_results$info %in% DHS_merged$V8] <- 'Other DHS'
VSMC_results$group[VSMC_results$info %in% DHS_merged$V8[DHS_merged$enhancers != 'No']] <- 'Enhancer Core'
VSMC_results$group[VSMC_results$info %in% DHS_merged$V8[DHS_merged$promoters != 'No']] <- 'Promoter'
VSMC_results$group[VSMC_results$info %in% DHS_merged$V8[DHS_merged$dyadic != 'No']] <- 'Dyadic'
VSMC_results <- VSMC_results[VSMC_results$info %in% sentinel_snps$change,]
table(VSMC_results$group)

CM_results$group <- factor(CM_results$group, levels = c('Other','Other DHS','Enhancer Core','Dyadic','Promoter'))
VSMC_results$group <- factor(VSMC_results$group, levels = c('Other','Other DHS','Enhancer Core','Dyadic','Promoter'))


foo <- pairwise.wilcox.test(VSMC_results$logFC_comp, VSMC_results$group, p.adjust.method="none")
foo

library(ggpubr)
#c("Other", "Other DHS"), c("Other", "Enhancer Core"), c("Other", "Dyadic"),
my_comparisons <- list(
  c('Other','Promoter'), c('Other','Other DHS'),c('Other','Enhancer Core'),
  c('Other','Dyadic'))

ggplot(CM_results, aes(x = group, y = abs(logFC_comp))) +
  geom_boxplot(aes(fill=group), notch=F, outlier.size=0.5, outlier.alpha = 0.5) +
  ggtitle('') +
  #stat_summary(fun.y=mean, geom="point", shape=23, size=4)+
  ylab('MPRA Activity')+
  xlab('')+
  ylim(c(0,9))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.background = element_rect(fill = 'white', colour = 'grey50'))+
  scale_fill_brewer(palette="Reds") + theme_classic() +
  stat_compare_means(comparisons = my_comparisons, label.y = c(4, 5, 6,7))


install.packages("ggridges")
library(ggridges)

ggplot(VSMC_results, aes(x = abs(logFC_comp), y = group, fill=group)) +
  geom_density_ridges(rel_min_height = 0.005, quantile_lines = TRUE, alpha = 0.75,
                      quantiles = 2) +
  xlim(c(-1,5)) + theme_classic() +
  theme(legend.position = 'None')+
  scale_fill_brewer(palette="Reds") +
  ylab('')


