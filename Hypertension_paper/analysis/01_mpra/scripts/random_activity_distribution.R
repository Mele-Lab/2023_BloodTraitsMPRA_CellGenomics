##### analysis of random sequences #####
library(dplyr)
library(tidyr)
library(ggplot2)

#### index and activity 
vals_significance_vsmc <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/VSMC_vals.significance.4.txt', sep = '\t', header = T)
head(vals_significance_vsmc)
vals_significance_vsmc <- vals_significance_vsmc %>%
  separate(name, c("pos", "snp_info","info"), "__")

vals_significance <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/CM_vals.significance.5.txt', sep = '\t', header = T)
head(vals_significance)
vals_significance <- vals_significance %>%
  separate(name, c("pos", "snp_info","info"), "__")

 index <- read.table('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/design/Hypertension__pooled.index.txt', sep = '\t', header = T)
 head(index)

### merge 
index$dupe_info <- gsub('\\..*','',index$tile_id)

### select outliers #### 
random <- vals_significance_vsmc[vals_significance_vsmc$tile_type=='RANDOM',]
sd(random$VSMC)
out <- boxplot.stats(random$VSMC)$out
out_ind <- which(random$VSMC %in% c(out))
out_ind
View(random[out_ind,])

outliers <- random$element[out_ind]
write.table(outliers, '~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/02_activs/random_elements_active_cms.IQR.txt',
            sep = '\n',col.names = F, row.names = F, quote = F)

#random$zscore <- ((random$VSMC-mean(random$VSMC))/sd(random$VSMC))

ggplot(random, aes(x = tile_type, y = VSMC)) +
  #geom_violin(fill='#BADACF')+
  geom_violin(fill='dark grey')+
  geom_boxplot(width=0.25)+
  geom_point(data = random[random$VSMC_padj<0.05,], color='red')+
  ylab("MPRA Activity") +
  xlab("") +
  theme(legend.position="top", legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA)) +
  #labs(color="", fill = "")+
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))+
  ylim(c(0,6))

### correlate number of TF w/ activity ######
## load TF FIMO
Tf <- read.delim('~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/03_fimo/FIMO_results/hum_tf/fimo.txt')
head(Tf)
Tf_random <- Tf[grep('RANDOM',Tf$sequence_name),]
head(Tf_random)

motif_info_dir = "~/marenostrum/Projects/Breast/ANALYSIS/Hypertension/data/03_fimo/Kaia_FIMO/"
motif_map_f = paste0(motif_info_dir, "/curated_motif_map.txt")
motif_info_f = paste0(motif_info_dir,"/motif_info.txt")

motif_map <- read.delim(motif_map_f)
motif_info <- read.delim(motif_info_f)

# only analyze the "best" motifs as determined by lambert et al
best_motifs <- motif_info[!is.na(motif_info$Best.Motif.s....Figure.2A.),]
nrow(best_motifs)

best_motifs$short_id <- gsub('\\..*','',best_motifs$CIS.BP.ID)
head(best_motifs)

Tf_random_best <- Tf_random[Tf_random$X..motif_id %in% best_motifs$short_id,]
Tf_random_number <- as.data.frame(table(Tf_random_best$sequence_name))
Tf_random_number_all <- as.data.frame(table(Tf_random$sequence_name))


### now correlate with activity 
vsmc_activity <- merge(vals_significance_vsmc, index[,c("element","tile_type","dupe_info")]%>%distinct(), by = 'element')
cm_activity <- merge(vals_significance, index[,c("element","tile_type","dupe_info")]%>%distinct(), by = 'element')
random_vsmc <- vsmc_activity[vsmc_activity$tile_type.y=='RANDOM',]
random_cm <- cm_activity[cm_activity$tile_type.y=='RANDOM',]

head(Tf_random_number)
random_cm$Var1 <- paste0('random_sequence_',random_cm$dupe_info.y, '_RANDOM')
active_cm_tf <- merge(random_cm,Tf_random_number, by = 'Var1')

random_vsmc$Var1 <- paste0('random_sequence_',random_vsmc$dupe_info.y, '_RANDOM')
active_vsmc_tf <- merge(random_vsmc,Tf_random_number_all, by = 'Var1')


library("ggpubr")
ggscatter(active_vsmc_tf[active_vsmc_tf$VSMC<6,], x = "VSMC", y = "Freq",  color='#F6552F',
          add = "reg.line", conf.int = F, 
          cor.coef = T, cor.method = "pearson",
          xlab = "MPRA activity", ylab = "Predicted TFBS")

