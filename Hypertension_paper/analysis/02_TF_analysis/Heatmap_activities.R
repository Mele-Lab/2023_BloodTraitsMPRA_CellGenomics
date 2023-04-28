##### heatmap for TF MPRA activity explained --------
library(pheatmap)

#### read data ----
sig_motifs_cm <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/TF_Analysis/sig_motifs.txt', sep = '\t', header = T)
sig_motifs_vsmc <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/TF_Analysis/sig_motifs_VSMC.txt', sep = '\t', header = T)

## merge results ------

all_motifs_results <- merge(sig_motifs_cm, sig_motifs_vsmc, by='short_id', all.x=TRUE, all.y=TRUE)
all_motifs_results <- all_motifs_results %>% 
  group_by(HGNC.symbol.x) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  as.data.frame()  


#### create heatmap -------
rownames(all_motifs_results) <- all_motifs_results$HGNC.symbol.x
all_motifs_results_f <- all_motifs_results[,c('rsq.x','rsq.y')]
colnames(all_motifs_results_f) <- c('CM','VSMC')
all_motifs_results_filt <- all_motifs_results_f[all_motifs_results_f$CM > 0.001 | all_motifs_results_f$VSMC > 0.001,]
all_motifs_results_filt <- all_motifs_results_filt[!is.na(all_motifs_results_filt$VSMC),]

pheatmap(all_motifs_results_filt,main = "pheatmap default",
         cluster_cols = F, color=colorRampPalette(c("white", "dark blue"))(50))
