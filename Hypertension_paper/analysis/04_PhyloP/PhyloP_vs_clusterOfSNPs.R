#### correlate PhyloPscore with cluster of SNPs #####
library(dplyr)
library(plyr)
library(tidyr)
library(ggplot2)

### read data ----
phylop_scores <- read.csv('../../analysis/04_PhyloP/PhyloP_scores_summary_tile_46_placenta_window_a_b.txt',
                            sep = '\t', header = T)

### read differential results ----
CM_results <- read.table('../../analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_ordered.CM.005.new_back.txt', sep='\t', header = T)
VSMC_results <- read.table('../../analysis/01_mpra/06_process_mpranalyze_comp/results_comparisons_ordered.VSMC.005.new_back.txt', sep='\t', header = T)

#### parse data ----
head(CM_results)
diff_snps <- unique(CM_results$snp[CM_results$fdr_comp<0.05])
diff_snps_vsmc <- unique(VSMC_results$snp[VSMC_results$fdr_comp<0.05])

cluster_of_snps <- as.data.frame(table(phylop_scores$sentinel))
colnames(cluster_of_snps) <- c('sentinel','Freq')

#### merge data ----
phylop_scores_freq <- merge(phylop_scores, cluster_of_snps, all.x=T)
phylop_scores_freq$diffCM <- 'No Diff'
phylop_scores_freq$diffCM[phylop_scores_freq$snp %in% diff_snps] <- 'Diff'
phylop_scores_freq$diffVSMC <- 'No Diff'
phylop_scores_freq$diffVSMC[phylop_scores_freq$snp %in% diff_snps_vsmc] <- 'Diff'

#### do some plots ----
## first summarise phylop score per region 

phylop_scores_freq_g <- phylop_scores_freq %>% 
  group_by(sentinel) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  as.data.frame()

ggscatter(phylop_scores_freq_g[phylop_scores_freq_g$Freq < 200,], x = "Freq", y = "mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Clustering", ylab = "mean(PhyloP)")

### summarise differential plots ####
cluster_of_snps <- as.data.frame(table(phylop_scores_freq$sentinel[phylop_scores_freq$diffCM=='Diff']))
colnames(cluster_of_snps) <- c('sentinel','Freq_CM')

phylop_scores_freq <- merge(phylop_scores_freq, cluster_of_snps, all.x=T)

cluster_of_snps <- as.data.frame(table(phylop_scores_freq$sentinel[phylop_scores_freq$diffVSMC=='Diff']))
colnames(cluster_of_snps) <- c('sentinel','Freq_VSMC')

phylop_scores_freq <- merge(phylop_scores_freq, cluster_of_snps, all.x=T)

phylop_scores_freq_g <- phylop_scores_freq %>% 
  group_by(sentinel,diffVSMC) %>% 
  summarise_if(is.numeric, mean, na.rm = TRUE) %>%
  as.data.frame()

library("ggpubr")
ggscatter(phylop_scores_freq_g[phylop_scores_freq_g$diffCM == 'Diff',], x = "Freq_CM", y = "mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Clustering", ylab = "mean(PhyloP)")
ggscatter(phylop_scores_freq_g[phylop_scores_freq_g$diffCM == 'Diff' & phylop_scores_freq_g$Freq_CM<100,], x = "Freq_CM", y = "mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Clustering", ylab = "mean(PhyloP)")

ggscatter(phylop_scores_freq_g[phylop_scores_freq_g$diffVSMC == 'Diff',], x = "Freq_VSMC", y = "mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Clustering", ylab = "mean(PhyloP)")
ggscatter(phylop_scores_freq_g[phylop_scores_freq_g$diffVSMC == 'Diff' & phylop_scores_freq_g$Freq_VSMC<20,], x = "Freq_VSMC", y = "mean", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Clustering", ylab = "mean(PhyloP)")


#### phylop scores differential vs not differential ----
phylop_scores$diffCM <- 'No Diff'
phylop_scores$diffCM[phylop_scores$snp %in% diff_snps] <- 'Diff'
phylop_scores$diffVSMC <- 'No Diff'
phylop_scores$diffVSMC[phylop_scores$snp %in% diff_snps_vsmc] <- 'Diff'

phylop_scores$sig_both <- 'no sig'
phylop_scores$sig_both[phylop_scores$diffCM == 'Diff'] <- 'sig CM'
phylop_scores$sig_both[phylop_scores$diffVSMC == 'Diff'] <- 'sig VSMC'
phylop_scores$sig_both[phylop_scores$diffCM == 'Diff' & phylop_scores$diffVSMC == 'Diff'] <- 'sig both'

library(plyr)
mu <- ddply(phylop_scores, "sig_both", summarise, grp.mean=mean(mean, na.rm = TRUE))
head(mu)

p<-ggplot(phylop_scores, aes(x=mean, color=sig_both)) +
  geom_density()+
  xlim(c(-0.1,1.3))+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=sig_both),
             linetype="dashed") + theme_minimal()
p


library(ggpubr)
p <- ggboxplot(phylop_scores, x = "sig_both", y = "mean",
               color = "sig_both", palette = "jco")
p +   stat_compare_means(method = "anova", label.y = 1.3)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "no sig")

mean(phylop_scores$mean[phylop_scores$sig_both!='no sig'], na.rm = T)
mean(phylop_scores$mean[phylop_scores$sig_both=='no sig'], na.rm = T)

wilcox.test(phylop_scores$mean[phylop_scores$sig_both!='no sig'],phylop_scores$mean[phylop_scores$sig_both=='no sig'],
            alternative = 'less')


