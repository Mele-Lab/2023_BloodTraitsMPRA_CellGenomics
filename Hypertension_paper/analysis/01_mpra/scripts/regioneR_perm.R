###### more differential in repeats? -------

#### read data -----
sentinel_snps = read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/sentinel_trait_ld_snps.txt', sep = '\t', header = T)
head(sentinel_snps)
sentinel_snps = sentinel_snps[sentinel_snps$Trait %in% c('DBP','PP','SBP'),]

#### activity and differential results -----
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

active_snps_cm <- unique(vals_significance_cm$snp_info[vals_significance_cm$CM_padj<0.05])
active_snps_vsmc <- unique(vals_significance_vsmc$snp_info[vals_significance_vsmc$VSMC_padj<0.05])

CM_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_ordered.CM.005.new_back.txt', sep='\t', header = T)
VSMC_results <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_ordered.VSMC.005.new_back.txt', sep='\t', header = T)

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

### repeat sequences -----
### RepeatMasker ####
repeat_masker <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/snps_repeat_masker.info.txt', sep=' ', header=T)
repeat_masker$query <- gsub('_.*','',repeat_masker$query)
repeat_masker <- repeat_masker %>% distinct()
sentinel = merge(sentinel_snps,repeat_masker, by.x= 'coord', by.y= 'query')

### we need to get info in regions -----
head(sentinel_snps)
sentinel_snps_cord <- as.list(sentinel_snps$coord)
sentinel_snps_gr = lapply(sentinel_snps_cord, function (x) {res=strsplit(x, ':')}) %>%
  unlist %>%
  matrix(ncol=3, byrow=T) %>%
  as.data.frame %>%
  dplyr::select(chrom=V1, start=V2, end=V3) %>%
  makeGRangesFromDataFrame

repeat_cord <- as.list(repeat_masker$query)
repeat_gr = lapply(repeat_cord, function (x) {res=strsplit(x, ':')}) %>%
  unlist %>%
  matrix(ncol=3, byrow=T) %>%
  as.data.frame %>%
  dplyr::select(chrom=V1, start=V2, end=V3) %>%
  makeGRangesFromDataFrame

diff_cord <- as.list(sentinel_snps$coord[sentinel_snps$DiffVSMC == 'Diff'])
diff_gr = lapply(diff_cord, function (x) {res=strsplit(x, ':')}) %>%
  unlist %>%
  matrix(ncol=3, byrow=T) %>%
  as.data.frame %>%
  dplyr::select(chrom=V1, start=V2, end=V3) %>%
  makeGRangesFromDataFrame

##### make test -------
library(regioneR)
pt <- permTest(A=diff_gr, B=repeat_gr, randomize.function=resampleRegions, universe=sentinel_snps_gr,
               evaluate.function=numOverlaps, ntimes=5000,verbose=FALSE)
summary(pt)
pt
plot(pt)
